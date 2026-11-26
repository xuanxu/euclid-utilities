# euclid_view

Open the image of any object in ESA's Euclid Data Release 1 (DR1) from its `object_id`.

The tool looks the object up in the DR1 MER catalogues on the public
[ESA Euclid archive](https://eas.esac.esa.int/), downloads a cutout of the mosaic it was
detected in, renders it as a PNG and opens it in your default image viewer.
No archive account is needed.

## Requirements

- [uv](https://docs.astral.sh/uv/). The script declares its dependencies inline,
  and `uv` installs them on first run.
- Network access to `eas.esac.esa.int`.

Without `uv`, install the dependencies yourself and run the script with Python 3.11+:

```sh
pip install astropy matplotlib numpy
python euclid_view.py <object_id>
```

## Usage

```sh
./euclid_view.py <object_id> [options]
```

| Option | Description |
|---|---|
| `-b`, `--band BAND` | `VIS`, `Y`, `J`, `H`, or an external band such as `DECAM_g`. Default: `VIS`, or `H` for objects only detected in the near-infrared. |
| `-s`, `--size SIZE` | Field of view in arcsec, up to 3600. Default: scaled to the object (10–120″). |
| `-o`, `--output PATH` | Where to save the PNG. Default: the system temp directory. |
| `-m`, `--match N` | Which detection to show when the ID appears more than once. Default: 1. |
| `--no-open` | Save the PNG without opening it. |

The image is shown north up and east left, with axes in arcsec from the object.
Orange ticks point to the object's catalogue position.

## Examples

Open an object (DR1 IDs are often negative; pass them as-is):

```sh
./euclid_view.py -683484103657468118
```

Zoom in on the Einstein ring around NGC 6505:

```sh
./euclid_view.py 2677813145655310071 -s 12
```

View the same object in the NISP H band, or in a ground-based band:

```sh
./euclid_view.py -683484103657468118 --band H
./euclid_view.py -705289463376508502 --band DECAM_g
```

Pick another detection of an ID that appears more than once (common in the deep fields):

```sh
./euclid_view.py 2631515946654747504 --match 2
```

Save the image without opening it:

```sh
./euclid_view.py -683484103657468118 --no-open -o galaxy.png
```

The command prints the object's position, the band and field of view, and where the PNG was saved:

```
RA 68.348410  Dec -65.746812  VIS  57"
/var/folders/.../euclid_-683484103657468118_VIS.png
```

## How it works

1. **Find the object.** It searches the five DR1 MER catalogues, in this order:
   `wide_survey`, `deep_survey`, `wide_mode`, `deep_mode` and `deep_visits`.
   If the ID appears more than once, the tool lists every detection with its
   catalogue and patch, and `--match` picks one.
2. **Find its mosaic.** A sky tile can hold dozens of mosaics per band, one per
   visit or patch. The tool picks those from the same pipeline run as the
   catalogue entry, so the image is the one the object was measured on.
3. **Cut out and render.** The archive's cutout service returns a small FITS
   image, which is drawn with an asinh stretch.

## Notes

- **Only DR1 IDs work.** IDs from the earlier Q1 release are not found, because
  DR1 assigned new IDs. Search DR1 near the same position instead: an ID encodes
  its coordinates, so `2677813028655307424` is near RA 267.7813028°, Dec +65.5307424°.
- **Contrast is per image.** The brightness range is set from the object itself,
  so don't compare brightness between images by eye.
- **Automatic framing can be too wide.** It uses the catalogue's Kron radius,
  which can be inflated for faint or blended sources. Set `--size` if needed.
- Tested on macOS. On Linux the image opens with `xdg-open`.

# euclid_discs

Search DR1 for edge-on protoplanetary discs: young stars hidden behind a dark dust
lane, with two bright scattered-light lobes on either side (a "butterfly").

```sh
./euclid_discs.py                                  # the dustiest sightlines in the catalogue
./euclid_discs.py --ra 68.348 --dec -65.747 -r 0.2 # a region of your choice
```

It writes `discs.csv`, every scored source ranked, and `discs.png`, a contact sheet of
the best candidates with orange ticks along each candidate's dark lane. Look at one more
closely with `./euclid_view.py <object_id>`.

| Option | Description |
|---|---|
| `--ra`, `--dec`, `-r`, `--radius` | Search a cone (degrees). Default: no region; select sightlines by dust instead. |
| `--min-ebv E` | Without a region, keep sources with Planck E(B-V) above this. Default: 0.2. |
| `-c`, `--catalogue NAME` | DR1 MER catalogue to search. Default: `wide_survey`. |
| `-n`, `--limit N` | Sources to download and score. Default: 300. |
| `-t`, `--top N` | Candidates on the contact sheet. Default: 25. |
| `--max-point-like P` | Maximum `point_like_prob`. Default: 0.1. |
| `--min-area`, `--max-area` | Isophotal area in VIS pixels (0.01 arcsec²). Default: 200–20000. |
| `-o`, `--output PREFIX` | Prefix for the `.csv` and `.png`. Default: `discs`. |
| `--no-open` | Save the contact sheet without opening it. |

## How the search works

1. **Pre-select.** Resolved (`point_like_prob` low), non-spurious, VIS-detected sources
   1.5–15″ across, from a cone or from the dustiest sightlines: discs sit in molecular
   clouds, which show up as high E(B-V).
2. **Score.** For each source, a VIS cutout is downloaded and the source is collapsed
   onto the axis across every candidate lane orientation. A butterfly shows a bright peak
   on each side with a dip between them. The score is `dip × √balance`, where `dip` is
   how dark the lane is relative to the fainter lobe (0–1) and `balance` is the
   brightness ratio of the lobes (0–1). Orientations are rejected when the lane is
   off-centre, the brightest pixel lies in the lane (a galaxy bulge), or a lobe is
   point-like or shorter than the lobe separation (a pair of neighbouring sources).
3. **Rank.** Sources within 2″ of a better one are dropped, because a disc's two lobes
   are often catalogued as separate objects.

A well-defined butterfly scores above about 0.5. Normal galaxies and stars score near
0, but close galaxy pairs and clumpy galaxies can still score 0.2–0.4, so check the
contact sheet by eye.

Cutouts are downloaded one at a time: the archive's cutout service returns the wrong
images when it gets concurrent requests. A search of 300 sources takes about 3 minutes.
