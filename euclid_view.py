#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["astropy>=6", "matplotlib>=3.8", "numpy"]
# ///
"""View the image of an object from ESA's Euclid Data Release 1 (DR1).

Looks the object up in the DR1 MER catalogues of the public ESA Euclid archive,
downloads a cutout of the mosaic it was detected in, and opens it as a PNG.
"""

import argparse
import csv
import io
import os
import subprocess
import sys
import tempfile
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.visualization import AsinhStretch, ImageNormalize
from astropy.wcs import WCS
from matplotlib.figure import Figure

ARCHIVE = "https://eas.esac.esa.int"
TAP_URL = f"{ARCHIVE}/tap-server/tap/sync"
CUTOUT_URL = f"{ARCHIVE}/sas-cutout/cutout"

# In order of preference when an object_id appears in more than one.
CATALOGUES = ["wide_survey", "deep_survey", "wide_mode", "deep_mode", "deep_visits"]
BAND_ALIASES = {"Y": "NIR_Y", "J": "NIR_J", "H": "NIR_H"}

MER_PIXEL_SCALE = 0.1  # arcsec; all MER mosaics share the VIS pixel grid
MIN_AUTO_SIZE, MAX_AUTO_SIZE = 10.0, 120.0  # arcsec
MAX_SIZE = 3600.0  # arcsec; the cutout service accepts radii up to 0.5 deg


class EuclidError(Exception):
    """A problem worth reporting to the user without a traceback."""


def fetch(url, params, timeout=120):
    try:
        with urllib.request.urlopen(
            f"{url}?{urllib.parse.urlencode(params)}", timeout=timeout
        ) as response:
            return response.read()
    except urllib.error.HTTPError as err:
        raise EuclidError(f"Euclid archive returned HTTP {err.code} for {url}") from err
    except (urllib.error.URLError, TimeoutError) as err:
        raise EuclidError(f"could not reach the Euclid archive: {err}") from err


def query(adql):
    body = fetch(
        TAP_URL,
        {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": adql},
    ).decode()
    if body.lstrip().startswith("<"):  # errors come back as a VOTable
        raise EuclidError(f"archive query failed:\n{body.strip()}")
    return list(csv.DictReader(io.StringIO(body)))


def find_object(object_id):
    """Return every DR1 detection of object_id, preferred catalogue first."""
    rows = query(
        " UNION ALL ".join(
            f"SELECT '{name}' AS catalogue, right_ascension, declination, tile_index,"
            f" basic_download_data_oid, vis_det, kron_radius, patch_id_list"
            f" FROM dr1.mer_catalogue_{name} WHERE object_id={object_id}"
            for name in CATALOGUES
        )
    )
    return sorted(
        rows,
        key=lambda r: (
            CATALOGUES.index(r["catalogue"].strip()),
            int(r["basic_download_data_oid"]),
        ),
    )


def pipeline_run(product_id):
    return product_id.rsplit("-compressed_", 1)[0]


def find_mosaics(detection):
    """Return {band: path} for the mosaics this detection was measured on.

    A tile can hold dozens of mosaics per band (one per visit or patch), so match
    on the pipeline run that produced both the mosaics and the catalogue.
    """
    catalogue = query(
        "SELECT product_id FROM dr1.mer_final_catalogue"
        f" WHERE basic_download_data_oid={int(detection['basic_download_data_oid'])}"
    )
    if not catalogue:
        raise EuclidError("the object's source catalogue is missing from the archive")
    run = pipeline_run(catalogue[0]["product_id"])
    mosaics = query(
        "SELECT filter_name, file_path, file_name, product_id FROM dr1.mosaic_product"
        f" WHERE tile_index={int(detection['tile_index'])}"
    )
    return {
        m["filter_name"]: f"{m['file_path']}/{m['file_name']}"
        for m in mosaics
        if pipeline_run(m["product_id"]) == run
    }


def pick_band(requested, detection, mosaics):
    if not mosaics:
        raise EuclidError(f"no mosaics found for tile {detection['tile_index']}")
    if requested is None:
        # NIR-only detections are invisible in VIS.
        requested = "VIS" if detection["vis_det"] == "1" else "NIR_H"
    available = {name.upper(): name for name in mosaics}
    band = requested.upper()
    band = BAND_ALIASES.get(band, band)
    if band not in available:
        raise EuclidError(
            f"no {requested} image for this object (available: {', '.join(sorted(mosaics))})"
        )
    return available[band]


def auto_size(detection):
    """Pick a field of view (arcsec) that frames the object."""
    try:
        size = 6 * float(detection["kron_radius"]) * MER_PIXEL_SCALE
    except ValueError:
        return MIN_AUTO_SIZE
    if not np.isfinite(size):
        return MIN_AUTO_SIZE
    return min(max(size, MIN_AUTO_SIZE), MAX_AUTO_SIZE)


def download_cutout(path, ra, dec, size):
    body = fetch(
        CUTOUT_URL, {"FILEPATH": path, "POS": f"CIRCLE,{ra},{dec},{size / 2 / 3600}"}
    )
    if not body.startswith(b"SIMPLE"):
        raise EuclidError("the cutout service did not return a FITS image")
    return body


def render(cutout, ra, dec, size, title, subtitle, png_path):
    with fits.open(io.BytesIO(cutout)) as hdus:
        data = np.asarray(hdus[0].data, dtype=float)
        wcs = WCS(hdus[0].header)

    data[data == 0] = np.nan  # exactly zero means no coverage
    pixels = data[np.isfinite(data)]
    if pixels.size == 0:
        raise EuclidError("the image has no data at the object's position")

    x, y = (float(v) for v in wcs.wcs_world2pix(ra, dec, 0))
    ny, nx = data.shape
    scale = MER_PIXEL_SCALE

    # Sky level from the whole frame, peak from around the target so that a bright
    # neighbour does not wash out a faint object.
    sky = np.median(pixels)
    noise = 1.4826 * np.median(np.abs(pixels - sky))
    yy, xx = np.ogrid[:ny, :nx]
    target = data[np.hypot(xx - x, yy - y) <= size / scale / 6]
    target = target[np.isfinite(target)]
    peak = np.percentile(target, 99.5) if target.size else np.percentile(pixels, 99.8)
    norm = ImageNormalize(
        vmin=sky - noise,
        vmax=max(peak, sky + 5 * noise),
        stretch=AsinhStretch(0.1),
        clip=True,
    )

    # Axes as offsets from the object; mosaics are north up, east left.
    extent = [(x + 0.5) * scale, (x - nx + 0.5) * scale, (-0.5 - y) * scale, (ny - 0.5 - y) * scale]

    fig = Figure(figsize=(8, 8.4), dpi=120, facecolor="white", layout="constrained")
    ax = fig.add_subplot()
    ax.set_facecolor("black")
    ax.imshow(data, origin="lower", cmap="gray", norm=norm, extent=extent, interpolation="nearest")
    ax.set_xlim(size / 2, -size / 2)
    ax.set_ylim(-size / 2, size / 2)
    ax.set_xlabel("ΔRA (arcsec, east ←)")
    ax.set_ylabel("ΔDec (arcsec, north ↑)")
    ax.set_title(subtitle, fontsize=9, color="0.35")
    fig.suptitle(title, fontsize=13)
    # Corner ticks mark the target without covering it.
    for (dx, dy) in ((1, 0), (0, 1)):
        ax.plot(
            [dx * size * 0.06, dx * size * 0.12],
            [dy * size * 0.06, dy * size * 0.12],
            color="#ffb000",
            lw=1.5,
        )
    fig.savefig(png_path)


def open_file(path):
    if sys.platform == "win32":
        os.startfile(path)
        return
    opener = "open" if sys.platform == "darwin" else "xdg-open"
    try:
        subprocess.run([opener, str(path)], check=True)
    except (OSError, subprocess.CalledProcessError):
        print(f"could not open a viewer; the image is at {path}", file=sys.stderr)


def positive_size(text):
    value = float(text)
    if not 0 < value <= MAX_SIZE:
        raise argparse.ArgumentTypeError(f"must be between 0 and {MAX_SIZE:.0f} arcsec")
    return value


def parse_args():
    parser = argparse.ArgumentParser(
        description="Open the Euclid DR1 image of an object.",
    )
    parser.add_argument("object_id", type=int, help="DR1 MER object_id (may be negative)")
    parser.add_argument("-b", "--band", help="VIS, Y, J, H, or an external band such as DECAM_g"
                        " (default: VIS, or H for objects only detected in NIR)")
    parser.add_argument("-s", "--size", type=positive_size,
                        help="field of view in arcsec (default: scaled to the object)")
    parser.add_argument("-o", "--output", type=Path, help="where to save the PNG (default: temp dir)")
    parser.add_argument("-m", "--match", type=int, default=1,
                        help="which detection to show if the ID has several (default: 1)")
    parser.add_argument("--no-open", action="store_true", help="save the PNG without opening it")
    return parser.parse_args()


def run(args):
    detections = find_object(args.object_id)
    if not detections:
        raise EuclidError(f"object {args.object_id} not found in Euclid DR1")
    if not 1 <= args.match <= len(detections):
        raise EuclidError(f"--match must be between 1 and {len(detections)} for this object")
    if len(detections) > 1:
        names = ", ".join(
            f"{i}={d['catalogue'].strip()} patch {d['patch_id_list'].strip('() ,')}"
            for i, d in enumerate(detections, 1)
        )
        print(f"{len(detections)} detections of this ID ({names}); showing {args.match}. "
              "Use --match to pick another.", file=sys.stderr)

    detection = detections[args.match - 1]
    ra, dec = float(detection["right_ascension"]), float(detection["declination"])
    mosaics = find_mosaics(detection)
    band = pick_band(args.band, detection, mosaics)
    size = args.size or auto_size(detection)

    output = args.output or Path(tempfile.gettempdir()) / f"euclid_{args.object_id}_{band}.png"
    cutout = download_cutout(mosaics[band], ra, dec, size)
    render(
        cutout, ra, dec, size,
        title=f"Euclid DR1 {args.object_id}",
        subtitle=f"{band} · {detection['catalogue'].strip()} · RA {ra:.5f}° Dec {dec:+.5f}° · {size:.0f}″",
        png_path=output,
    )
    print(f"RA {ra:.6f}  Dec {dec:+.6f}  {band}  {size:.0f}\"\n{output}")
    if not args.no_open:
        open_file(output)


def main():
    try:
        run(parse_args())
    except EuclidError as err:
        sys.exit(f"error: {err}")


if __name__ == "__main__":
    main()
