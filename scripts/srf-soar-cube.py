"""
Investigate spectral reponse function (SRF) of the SOAR spectra

It seems to have exponential wings, which we want to include in our fitting profile

Use the opportunity to learn how to use the spectral-cube and pyspeckit libraries

"""

from pathlib import Path
from spectral_cube import SpectralCube
import astropy.units as u
import pyspeckit
from pyspeckit.spectrum.units import SpectroscopicAxis


# Folder that contains this script
script_dir = Path(__file__).resolve().parent
# Assume we are a top-level subfolder of the project
project_dir = script_dir.parent

# Location of data cube: assume that it is in sibling folder of this project
cube_dir = project_dir.parent / "tarantula-soar-cubo"

# This is the continuum-subtracted cube in LSRK velocities
prefix = "soar-30dor-ha"

cube = SpectralCube.read(cube_dir / f"{prefix}-cube-lsrk-csub.fits")

# print(cube.spectral_axis.to(u.km / u.s))

# Extract spectrum to fit, which has no obvious high-velocity components
region_strs = [
    'icrs; box(84.6877883,-69.0854684,2.0",2.0",0.0)',
    'icrs; box(84.6991375,-69.0866617,2.0",2.0",0.0)',
    'icrs; box(84.7239776,-69.0922928,12.0",5.0",55)',
    'icrs; box(84.6607049,-69.1081513,4.0",2.0",55)',
    'icrs; box(84.6638244,-69.0943192,4.0",5.5",0)',
]
region_str = region_strs[0]

sub_cube = cube.subcube_from_ds9region(region_str)
# Convert to 1D spectrum
spectrum = sub_cube.mean(axis=(1, 2))

#
# Transfer spectrum to pyspeckit
#
xarr = SpectroscopicAxis(
    spectrum.spectral_axis.to(u.km / u.s),
    velocity_convention="optical",
    refX=6562.78 * u.AA,
)
sp = pyspeckit.Spectrum(
    xarr=xarr,
    data=spectrum.value,
    header={},
    unit=spectrum.unit.to_string(),
)

sp.plotter()

sp.specfit(
    fittype="gaussian",
    guesses=[spectrum.max().value, 275.0, 20.0],
)

sp.specfit.plot_fit()
sp.plotter.axis.set_yscale("log")
sp.plotter.axis.set_ylim(1.0, None)

sp.plotter.savefig(script_dir / "srf-gauss-fit.jpg")


import numpy as np
from tarantula_soar.soar_profiles import soar_moffat_fitter, soar_moffat_bfix_fitter

# Use a non-uniform weight so that the wings get taken into account in the fit
y = spectrum.value
y_floor = 20.0  # for example: perhaps 3–5 times the wing noise
alpha = 0.5
scale = np.maximum(np.abs(y), y_floor)
error_eff = y_floor * (scale / y_floor) ** alpha

spp = pyspeckit.Spectrum(
    xarr=xarr,
    data=y,
    error=error_eff,
    header={},
    unit=spectrum.unit.to_string(),
)
spp.plotter()
spp.Registry.add_fitter("soar-moffat-bfix", soar_moffat_bfix_fitter(), 3)

spp.specfit(
    fittype="soar-moffat-bfix",
    guesses=[spectrum.max().value, 275.0, 20.0, 1.0],
)
spp.specfit.plot_fit()

spp.plotter.savefig(script_dir / "srf-custom-fit.jpg")

spp.plotter.axis.set_yscale("log")
spp.plotter.axis.set_ylim(1.0, None)

spp.plotter.savefig(script_dir / "srf-custom-fit-semilog.jpg")
