import numpy as np
import pyspeckit


def moffat_wings(
    x,
    amplitude,
    shift,
    width,
    beta,
    single_component_flux=False,
    return_components=False,
    return_hyperfine_components=False,
):
    x = np.asarray(x, dtype=float)
    dx = x - shift
    sigma = width

    return amplitude * (1.0 + dx**2 / (2.0 * beta * sigma**2)) ** (-beta)


BETA_FROZEN = 1.8  # From fitting narrow bright profiles with no high-velocity emission


def moffat_wings_frozen(
    x,
    amplitude,
    shift,
    width,
    single_component_flux=False,
    return_components=False,
    return_hyperfine_components=False,
):
    return moffat_wings(x, amplitude, shift, width, BETA_FROZEN)


def gaussian_exp_wings(
    x,
    amplitude,
    shift,
    width,
    wing_scale,
    single_component_flux=False,
    return_components=False,
    return_hyperfine_components=False,
):
    x = np.asarray(x, dtype=float)
    dx = x - shift
    sigma = width
    b = wing_scale

    phi = (sigma / b) ** 2 * (np.sqrt(1.0 + (b * dx / sigma**2) ** 2) - 1.0)
    return amplitude * np.exp(-phi)


def soar_gew_fitter():
    """
    Generator for SOAR SRF fitter class with Gausian+Exponential-Wings model
    """
    myclass = pyspeckit.models.model.SpectralModel(
        gaussian_exp_wings,
        4,
        parnames=["amplitude", "shift", "width", "wing_scale"],
        parlimited=[(False, False), (False, False), (True, False), (True, False)],
        parlimits=[(0, 0), (0, 0), (0, 0), (0, 0)],
        shortvarnames=("A", r"\Delta x", r"\sigma", "b"),
        centroid_par="shift",
        fitunit="km/s",
    )
    myclass.__name__ = "soar-gew"
    return myclass


def soar_moffat_fitter():
    """
    Generator for SOAR SRF fitter class with Moffat profile
    """
    myclass = pyspeckit.models.model.SpectralModel(
        moffat_wings,
        4,
        parnames=["amplitude", "shift", "width", "beta"],
        parlimited=[(False, False), (False, False), (True, False), (True, False)],
        parlimits=[(0, 0), (0, 0), (0, 0), (0, 0)],
        shortvarnames=("A", r"\Delta x", r"\sigma", r"\beta"),
        centroid_par="shift",
        fitunit="km/s",
    )
    myclass.__name__ = "soar-moffat"
    return myclass


def soar_moffat_bfix_fitter():
    """
    Generator for SOAR SRF fitter class with Moffat profile at fixed beta parameter
    """
    myclass = pyspeckit.models.model.SpectralModel(
        moffat_wings_frozen,
        3,
        parnames=["amplitude", "shift", "width"],
        parlimited=[(True, False), (False, False), (True, False)],
        parlimits=[(0, 0), (0, 0), (0, 0)],
        shortvarnames=("A", r"\Delta x", r"\sigma"),
        centroid_par="shift",
        fitunit="km/s",
    )
    myclass.__name__ = "soar-moffat-bfix"
    return myclass
