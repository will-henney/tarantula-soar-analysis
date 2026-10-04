"""
Use the tetrabloks algorithm to make spatially binned versions of the
SOAR Ha cube: 2x2, 4x4, 8x8, 16x16 binning of the native 0.2 arcsec
pixels

The cubes are saved at their smaller grid sizes rather than being
resampled back to the fine grid.  This makes them suitable for
accelerating the fitting of Gaussians

Will Henney 2026-10-04
"""

from astropy.io import fits
from pathlib import Path
import numpy as np
from tetrabloks import rebin_utils

# Folder that contains this script
script_dir = Path(__file__).resolve().parent
# Assume we are a top-level subfolder of the project
project_dir = script_dir.parent

# Location of data cube: assume that it is in sibling folder of this project
cube_dir = project_dir.parent / "tarantula-soar-cubo"

# This is the continuum-subtracted cube in LSRK velocities
prefix = "soar-30dor-ha"

# Load full-resolution cube
cube_hdu = fits.open(cube_dir / f"{prefix}-cube-lsrk-csub.fits")[0]

# List of resampling factors and minimum good pixels per level
nlist = [2, 4, 8, 16]
mingoods = [2, 2, 2, 2]
nmax = max(nlist)

nk, nj, ni = cube_hdu.data.shape
print("Original shape:", cube_hdu.data.shape)
# tetrabloks is written for 2d images, so we will just loop over velocity slices
cube = np.stack([rebin_utils.pad_array(_, nmax) for _ in cube_hdu.data])
print("Padded shape:", cube.shape)

ew_map = rebin_utils.pad_array(fits.open(cube_dir / f"{prefix}-ew.fits")[0].data, nmax)
# Mask out pixels with EW(Ha) < 10 AA
m = ew_map > 10.0

# Use uniform weights for now
w = np.where(m, 1.0, 0.0)


def downsample_header(header: fits.Header) -> fits.Header:
    """Adjust WCS keywords for rebinning of image by factor of 2"""
    new_header = header.copy()
    new_header["CRPIX1"] /= 2.0
    new_header["CRPIX2"] /= 2.0
    new_header["CDELT1"] *= 2.0
    new_header["CDELT2"] *= 2.0
    return new_header


# rebin_utils.downsample() accepts a list of images, so we can do the entrire cube at once
imlist = [_ for _ in cube]
header = cube_hdu.header
for n, mingood in zip(nlist, mingoods):
    imlist, m, w = rebin_utils.downsample(imlist, m, weights=w, mingood=mingood)
    header = downsample_header(header)
    outfile = cube_dir / f"{prefix}-cube-lsrk-csub-bin{n:03d}.fits"
    print("Saving", outfile)
    fits.HDUList(
        [
            fits.PrimaryHDU(),
            fits.ImageHDU(data=np.stack(imlist), header=header, name="SCI"),
            fits.ImageHDU(data=np.stack([w] * nk), header=header, name="WHT"),
        ]
    ).writeto(outfile, overwrite=True)
