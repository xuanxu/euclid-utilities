#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["astropy>=6", "matplotlib>=3.8", "numpy", "scipy"]
# ///
"""Search Euclid DR1 for edge-on protoplanetary discs ("butterfly" nebulae).

An edge-on disc around a young star hides the star behind a dark dust lane, with
light scattered off the disc surfaces forming two bright lobes above and below it.
The search runs in two stages:

1. Pre-select resolved, non-spurious sources from a DR1 MER catalogue, favouring
   dusty sightlines (where young stars live) or a sky region you choose.
2. Download a VIS cutout of each and score how butterfly-like it is: two bright
   lobes of similar brightness separated by a dark lane through the centre.

The output is a CSV of every scored source, ranked, and a contact sheet of the best
candidates for inspection by eye. Open one with ./euclid_view.py <object_id>.
"""

import argparse
import csv
import functools
import io
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.visualization import AsinhStretch, ImageNormalize
from astropy.wcs import WCS
from matplotlib.figure import Figure
from scipy import ndimage

from euclid_view import (
    CATALOGUES,
    MER_PIXEL_SCALE,
    EuclidError,
    download_cutout,
    find_mosaics,
    open_file,
    query,
)

MIN_CUTOUT, MAX_CUTOUT = 5.0, 30.0  # arcsec
ANGLES = np.deg2rad(np.arange(0, 180, 7.5))
MIN_LOBE_SEPARATION = 3  # pixels; the VIS PSF is ~1.6 px FWHM
MIN_LOBE_LENGTH = 5  # pixels; longer than a point source
DUPLICATE_RADIUS = 2.0  # arcsec; deblended lobes of one disc show up as two sources

CSV_FIELDS = ["rank", "object_id", "catalogue", "ra", "dec", "score", "dip", "balance",
              "lane_pa", "lobe_separation", "gal_ebv", "kron_radius", "error"]


def select_sources(args):
    """Return catalogue rows for the sources worth looking at."""
    where = [
        "vis_det = 1",
        f"point_like_prob < {args.max_point_like}",
        "spurious_prob < 0.2",
        # Isophotal area in 0.01 arcsec² pixels: roughly 1.5″ to 15″ across.
        f"segmentation_area BETWEEN {args.min_area} AND {args.max_area}",
    ]
    if args.ra is not None:
        where.append(
            "1 = CONTAINS(POINT('ICRS', right_ascension, declination),"
            f" CIRCLE('ICRS', {args.ra}, {args.dec}, {args.radius}))"
        )
        order = "flux_detection_total DESC"
    else:
        # Without a region, look where Planck sees dust: discs live in molecular clouds.
        where.append(f"gal_ebv > {args.min_ebv}")
        order = "gal_ebv DESC"
    return query(
        f"SELECT TOP {args.limit} '{args.catalogue}' AS catalogue, object_id,"
        " right_ascension, declination, tile_index, basic_download_data_oid, vis_det,"
        " kron_radius, gal_ebv"
        f" FROM dr1.mer_catalogue_{args.catalogue}"
        f" WHERE {' AND '.join(where)} ORDER BY {order}"
    )


@functools.cache
def vis_mosaic(oid, tile_index):
    mosaics = find_mosaics({"basic_download_data_oid": oid, "tile_index": tile_index})
    if "VIS" not in mosaics:
        raise EuclidError(f"no VIS mosaic for tile {tile_index}")
    return mosaics["VIS"]


def cutout_size(source):
    try:
        size = 4 * float(source["kron_radius"]) * MER_PIXEL_SCALE
    except ValueError:
        return MIN_CUTOUT
    if not np.isfinite(size):
        return MIN_CUTOUT
    return min(max(size, MIN_CUTOUT), MAX_CUTOUT)


def load_image(cutout, ra, dec):
    with fits.open(io.BytesIO(cutout)) as hdus:
        data = np.asarray(hdus[0].data, dtype=float)
        x, y = WCS(hdus[0].header).wcs_world2pix(ra, dec, 0)
    # Guard against the cutout service handing back another request's image.
    ny, nx = data.shape
    if abs(x - (nx - 1) / 2) > 3 or abs(y - (ny - 1) / 2) > 3:
        raise EuclidError("the cutout service returned an image of the wrong position")
    data[data == 0] = np.nan  # exactly zero means no coverage
    pixels = data[np.isfinite(data)]
    if pixels.size < 100:
        raise EuclidError("the cutout has no data")
    sky = np.median(pixels)
    noise = 1.4826 * np.median(np.abs(pixels - sky))
    return np.nan_to_num(data - sky, nan=0.0), noise


def source_mask(image, noise):
    """Pixels belonging to the source at the centre of the cutout.

    Thresholds low on a smoothed image and closes small gaps, so that both lobes end
    up in one region even when the lane between them dips below the noise.
    """
    smooth = ndimage.gaussian_filter(image, 1.5)
    smooth_noise = noise / (2 * np.sqrt(np.pi) * 1.5)  # noise after the Gaussian filter
    mask = ndimage.binary_closing(smooth > 2 * smooth_noise)
    labels, count = ndimage.label(mask)
    if count == 0:
        return None
    ny, nx = image.shape
    centre = labels[ny // 2, nx // 2]
    if centre == 0:  # the centre lies in the lane: take the region nearest to it
        distances = ndimage.distance_transform_edt(labels == 0, return_indices=True)[1]
        centre = labels[distances[0][ny // 2, nx // 2], distances[1][ny // 2, nx // 2]]
    return labels == centre


def lobe_length(row):
    """Pixels along the lane at which a lobe is above half its peak."""
    return np.count_nonzero(row > 0.5 * row.max()) if row.size else 0


def butterfly_score(image, noise):
    """Measure how much the source looks like two lobes split by a dark lane.

    For each orientation of a candidate lane through the source's centroid, collapse
    the source onto the axis across the lane and look for a bright peak on each side
    with a dip between them. Returns the best orientation's measurements.
    """
    mask = source_mask(image, noise)
    result = {"score": 0.0, "dip": 0.0, "balance": 0.0, "lane_pa": np.nan,
              "lobe_separation": np.nan}
    if mask is None or mask.sum() < 30:
        return result

    smooth = ndimage.gaussian_filter(image, 1.0)
    weights = np.clip(smooth, 0, None) * mask
    yy, xx = np.indices(image.shape)
    total = weights.sum()
    if total <= 0:
        return result
    cx, cy = (weights * xx).sum() / total, (weights * yy).sum() / total
    dx, dy = xx[mask] - cx, yy[mask] - cy
    values = smooth[mask]
    brightest = np.argmax(values)

    for angle in ANGLES:
        # u runs along the lane, v across it.
        u = dx * np.cos(angle) + dy * np.sin(angle)
        v = -dx * np.sin(angle) + dy * np.cos(angle)
        half_length = np.percentile(np.abs(u), 80)
        strip = np.abs(u) <= half_length
        if strip.sum() < 20:
            continue
        bins = np.round(v[strip]).astype(int)
        offset = -bins.min()
        sums = np.bincount(bins + offset, weights=values[strip])
        counts = np.bincount(bins + offset)
        profile = np.where(counts > 0, sums / np.maximum(counts, 1), 0.0)
        centre = offset
        if centre < 2 or centre > profile.size - 3:
            continue

        below, above = profile[:centre], profile[centre + 1:]
        i_low = int(np.argmax(below))
        i_high = centre + 1 + int(np.argmax(above))
        low, high = profile[i_low], profile[i_high]
        separation = i_high - i_low
        if separation < MIN_LOBE_SEPARATION or min(low, high) < 3 * noise:
            continue
        # Two neighbouring blobs also give a dip; a disc's lobes are spread along the
        # lane, so they are longer than their separation.
        if 2 * half_length < separation:
            continue
        i_floor = i_low + int(np.argmin(profile[i_low:i_high + 1]))
        if abs(i_floor - centre) > max(1, separation / 4):  # lane off-centre
            continue
        # The star is hidden: a bright core in the lane means a bulge, not a disc.
        if abs(v[brightest] - (i_floor - offset)) < max(1, separation / 4):
            continue
        # Each lobe is a band along the lane, not a point source.
        rows = bins + offset
        if min(lobe_length(values[strip][rows == i]) for i in (i_low, i_high)) < max(
            MIN_LOBE_LENGTH, 0.7 * separation
        ):
            continue
        floor = profile[i_floor]
        dip = 1 - max(floor, 0) / min(low, high)
        balance = min(low, high) / max(low, high)
        score = dip * np.sqrt(balance)
        if score > result["score"]:
            result = {
                "score": score,
                "dip": dip,
                "balance": balance,
                # Position angle of the lane, east of north (images are north up, east left).
                "lane_pa": (90 + np.rad2deg(angle)) % 180,
                "lobe_separation": separation * MER_PIXEL_SCALE,
                "centre": (cx, cy),
                "angle": angle,
                "half_length": half_length,
            }
    return result


def score_source(source):
    """Download a source's VIS cutout and score it. Never raises."""
    ra, dec = float(source["right_ascension"]), float(source["declination"])
    row = {
        "object_id": source["object_id"],
        "catalogue": source["catalogue"].strip(),
        "ra": ra,
        "dec": dec,
        "gal_ebv": source["gal_ebv"],
        "kron_radius": source["kron_radius"],
        "error": "",
    }
    try:
        path = vis_mosaic(int(source["basic_download_data_oid"]), int(source["tile_index"]))
        cutout = download_cutout(path, ra, dec, cutout_size(source))
        image, noise = load_image(cutout, ra, dec)
        row.update(butterfly_score(image, noise))
        row["image"] = image
        row["noise"] = noise
    except (EuclidError, OSError, ValueError) as err:
        row.update(score=0.0, dip=0.0, balance=0.0, lane_pa=np.nan,
                   lobe_separation=np.nan, error=str(err))
    return row


def drop_duplicates(rows):
    """Keep only the best-scoring source within DUPLICATE_RADIUS of each position."""
    kept = []
    for row in rows:
        cos_dec = np.cos(np.deg2rad(row["dec"]))
        if all(
            np.hypot((row["ra"] - k["ra"]) * cos_dec, row["dec"] - k["dec"]) * 3600
            > DUPLICATE_RADIUS
            for k in kept
        ):
            kept.append(row)
    return kept


def write_csv(rows, path):
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_FIELDS, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({
                key: f"{value:.{7 if key in ('ra', 'dec') else 4}f}"
                if isinstance(value, float) else value
                for key, value in row.items()
            })


def contact_sheet(rows, path):
    columns = min(5, len(rows))
    lines = -(-len(rows) // columns)
    fig = Figure(figsize=(2.4 * columns, 2.6 * lines), dpi=120, facecolor="white",
                 layout="constrained")
    fig.suptitle("Euclid DR1 edge-on disc candidates (VIS)", fontsize=12)
    for i, row in enumerate(rows):
        ax = fig.add_subplot(lines, columns, i + 1)
        image, noise = row["image"], row["noise"]
        ny, nx = image.shape
        peak = np.percentile(image[ny // 4:3 * ny // 4, nx // 4:3 * nx // 4], 99.5)
        norm = ImageNormalize(vmin=-noise, vmax=max(peak, 5 * noise),
                              stretch=AsinhStretch(0.1), clip=True)
        ax.imshow(image, origin="lower", cmap="gray", norm=norm, interpolation="nearest")
        # Mark the lane just outside the source, so the lane itself stays visible.
        (cx, cy), angle, length = row["centre"], row["angle"], row["half_length"]
        for side in (-1, 1):
            r = np.array([1.3, 1.8]) * length * side
            ax.plot(cx + r * np.cos(angle), cy + r * np.sin(angle), color="#ffb000", lw=1)
        ax.set_xlim(-0.5, nx - 0.5)
        ax.set_ylim(-0.5, ny - 0.5)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(f"#{row['rank']}  score {row['score']:.2f}\n{row['object_id']}",
                     fontsize=7)
    fig.savefig(path)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Search Euclid DR1 for edge-on protoplanetary discs (butterfly-like"
        " images). Without a region, searches the dustiest sightlines in the catalogue.",
    )
    region = parser.add_argument_group("region (optional)")
    region.add_argument("--ra", type=float, help="centre RA in degrees")
    region.add_argument("--dec", type=float, help="centre Dec in degrees")
    region.add_argument("-r", "--radius", type=float, default=0.2,
                        help="search radius in degrees (default: 0.2)")
    parser.add_argument("-c", "--catalogue", choices=CATALOGUES, default="wide_survey",
                        help="DR1 MER catalogue to search (default: wide_survey)")
    parser.add_argument("--min-ebv", type=float, default=0.2,
                        help="minimum Planck E(B-V) when no region is given (default: 0.2)")
    parser.add_argument("--max-point-like", type=float, default=0.1,
                        help="maximum point_like_prob (default: 0.1)")
    parser.add_argument("--min-area", type=int, default=200,
                        help="minimum isophotal area in VIS pixels (default: 200)")
    parser.add_argument("--max-area", type=int, default=20000,
                        help="maximum isophotal area in VIS pixels (default: 20000)")
    parser.add_argument("-n", "--limit", type=int, default=300,
                        help="maximum number of sources to download and score (default: 300)")
    parser.add_argument("-t", "--top", type=int, default=25,
                        help="number of candidates on the contact sheet (default: 25)")
    parser.add_argument("-o", "--output", type=Path, default=Path("discs"),
                        help="prefix for the .csv and .png outputs (default: discs)")
    parser.add_argument("--no-open", action="store_true",
                        help="save the contact sheet without opening it")
    args = parser.parse_args()
    if (args.ra is None) != (args.dec is None):
        parser.error("--ra and --dec must be given together")
    return args


def run(args):
    sources = select_sources(args)
    if not sources:
        raise EuclidError("no sources match the selection")
    print(f"scoring {len(sources)} sources from dr1.mer_catalogue_{args.catalogue}",
          file=sys.stderr)

    # One at a time: the cutout service mixes up the images of concurrent requests.
    rows = []
    for done, source in enumerate(sources, 1):
        rows.append(score_source(source))
        print(f"\r  {done}/{len(sources)}", end="", file=sys.stderr)
    print(file=sys.stderr)

    failed = sum(1 for r in rows if r["error"])
    if failed:
        print(f"{failed} sources could not be scored; see the error column", file=sys.stderr)

    rows = drop_duplicates(sorted(rows, key=lambda r: r["score"], reverse=True))
    for rank, row in enumerate(rows, 1):
        row["rank"] = rank

    csv_path = args.output.with_suffix(".csv")
    write_csv(rows, csv_path)
    best = [r for r in rows[:args.top] if not r["error"] and r["score"] > 0]
    if not best:
        print(f"no butterfly-like sources found\n{csv_path}")
        return

    png_path = args.output.with_suffix(".png")
    contact_sheet(best, png_path)
    for row in best[:10]:
        print(f"{row['rank']:3d}  {row['object_id']:>20}  score {row['score']:.2f}"
              f"  RA {row['ra']:.5f}  Dec {row['dec']:+.5f}")
    print(f"{csv_path}\n{png_path}")
    if not args.no_open:
        open_file(png_path)


def main():
    try:
        run(parse_args())
    except EuclidError as err:
        sys.exit(f"error: {err}")


if __name__ == "__main__":
    main()
