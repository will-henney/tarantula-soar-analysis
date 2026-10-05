"""
Apply various preprocessing fixes to SOAR Ha spectral cube:

* WCS astrometric corrections
* Convert from wavelength to LSRK velocities
* Subtract continuum from cube
* Save 2D maps of continuum, summed line, and equivalent width
"""

from astropy.io import fits
from astropy.wcs import WCS
from pathlib import Path
import numpy as np
import re

# Location of data cube (assumed to be sibling folder to this project)
PATH = Path(__file__).resolve().parent.parent.parent / "tarantula-soar-cubo"

# Prefix for all new files
PREFIX = "soar-30dor-ha"

# WCS keywords that we might want to nuke
WCS_PATTERNS = [
    r"WCSDIM",
    r"WCSAXES",
    r"CTYPE\d+",
    r"CUNIT\d+",
    r"CRVAL\d+",
    r"CRPIX\d+",
    r"CDELT\d+",
    r"CD\d+_\d+",
    r"C\d+_\d+",
    r"PC\d+_\d+",
    r"CROTA\d+",
    r"LTV\d+",
    r"LTM\d+_\d+",
    r"WAT\d+_\d+",
]

# Load original cube
hdu = fits.open(PATH / "new_cube30Dor_SOAR_invertedY_NEWWAVELENGTH.fits")[0]

# Save the wavelength pixel size in Angstrom
dwav = hdu.header["CD3_3"]

# Get the by-hand corrections to the WCS (astrometry fix and conversion to VLSRK)
with open(PATH / "soar-vlsrk.wcs") as f:
    new_header = fits.Header.fromstring(f.read(), sep="\n")

# Remove old WCS from header
for key in list(hdu.header.keys()):
    if any(re.fullmatch(pattern, key) for pattern in WCS_PATTERNS):
        del hdu.header[key]
# Apply the new WCS to header
hdu.header.update(new_header)

#
# Make maps of line and continuum
#

# Limits of the line along velocity axis (0-based indices)
k1, k2 = 3, 37
# Take median of the far-blue and far-red edges to get the continuum
cont_map = np.median(
    np.concatenate(
        (hdu.data[:k1, ...], hdu.data[k2:, ...]),
        axis=0,
    ),
    axis=0,
)

# Continuum-subtracted cube
line_cube = hdu.data - cont_map[None, ...]

# Integrate the continuum-subtracted cube to get the line emission
line_map = np.trapezoid(line_cube[k1:k2, ...], dx=dwav, axis=0)

# Equivalent width of emission line
ew_map = line_map / cont_map

# Save the maps
wcs_2d = WCS(hdu.header).celestial
fits.PrimaryHDU(
    data=cont_map,
    header=wcs_2d.to_header(),
).writeto(PATH / f"{PREFIX}-cont.fits", overwrite=True)
fits.PrimaryHDU(
    data=line_map,
    header=wcs_2d.to_header(),
).writeto(PATH / f"{PREFIX}-sum.fits", overwrite=True)
fits.PrimaryHDU(
    data=ew_map,
    header=wcs_2d.to_header(),
).writeto(PATH / f"{PREFIX}-ew.fits", overwrite=True)

# Save the continuum-subtracted cube
fits.PrimaryHDU(
    data=line_cube,
    header=hdu.header,
).writeto(PATH / f"{PREFIX}-cube-lsrk-csub.fits", overwrite=True)
