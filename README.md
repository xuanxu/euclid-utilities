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
