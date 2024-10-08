"""Core module of the package."""

import logging
import importlib

import astropy.units as u
from gammapy.stats import WStatCountsStatistic
import numpy as np
from scipy import interpolate
from scipy.integrate import quad

from . import (
    LOW_ZENITH_PERFORMANCE,
    MID_ZENITH_PERFORMANCE,
)
from .io import load_ebl
from .spectral import crab_nebula_spectrum

__all__ = [
    "setup_logging",
    "check_input_configuration",
    "initialize_model",
    "observed_flux",
    "get_sed",
    "significance_li_ma",
    "prepare_data",
    "source_detection",
    "calculate",
    "get_horizon_stereo_profile",
]

logger = logging.getLogger(__name__)


def setup_logging(log_level, source_name):
    """
    Create a logger.

    The logger will have a console and file handler,
    saving to a file with the name of the required source.

    Parameters
    ----------
    log_level : `str`
        Logging level to use for the console handler.
        The file handler always uses ``DEBUG``
        (all calls are saved).
    source_name : `str`
        Name of the source provided by the user via the command
        line interface. A log file with the source name as
        prefix will be created at the output path.

    Returns
    -------
    logger : `~logging.Logger`
        Logger instance.
    """
    numeric_level = getattr(logging, log_level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError("Invalid log level: %s" % log_level)

    logger = logging.getLogger("iact_estimator")
    logger.setLevel(numeric_level)

    console_handler = logging.StreamHandler()
    console_handler.setLevel(numeric_level)

    file_handler = logging.FileHandler(
        f"{source_name}_iact_estimator.log",
        mode="w",
        encoding=None,
        delay=False,
        errors=None,
    )

    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    console_handler.setFormatter(formatter)
    file_handler.setFormatter(formatter)

    logger.addHandler(console_handler)
    logger.addHandler(file_handler)
    return logger


def get_horizon_stereo_profile(M1_data, M2_data):
    """Merge the horizon profiles into a stereo one.

    Parameters
    ----------
    M1_data : :class:`~astropy.table.QTable`
        Horizone profile seen from M1.
        The table must have columns
        'azimuth' and 'zenith', both in units
        of degrees.
    M2_data : :class:`~astropy.table.QTable`
        Horizone profile seen from M2.
        The table must have columns
        'azimuth' and 'zenith', both in units
        of degrees.

    Returns
    -------
    az : array-like
        Azimuth values within the maximum
        range of the separate measurements.
    zd : array-like
        Zenith angle values for the merged horizon
        profile.
    """
    m2_az = M2_data["azimuth"].value
    m2_zd = M2_data["zenith"].value
    m1_az = M1_data["azimuth"].value
    m1_zd = M1_data["zenith"].value

    if np.max(m2_az) > 360:
        m2_az = m2_az - 360.0
    if np.max(m1_az) > 360:
        m1_az = m1_az - 360.0

    min_az = np.min(np.concatenate((m2_az, m1_az)))
    max_az = np.max(np.concatenate((m2_az, m1_az)))

    az = np.linspace(min_az, max_az, 100)

    interp_M1_zd = np.interp(az, m1_az, m1_zd)
    interp_M2_zd = np.interp(az, m2_az, m2_zd)

    zd = np.minimum(interp_M1_zd, interp_M2_zd)

    return az * u.deg, zd * u.deg


def check_input_configuration(config):
    """
    Import and initialize a spectral model.

    Parameters
    ----------
    config : `dict`
        Configuration data in form of a Python dictionary.

    Returns
    -------
    valid : `bool`
        Initialized instance of a spectral model.
    """

    # Assume a valid configuration
    is_valid = True

    extension = u.Quantity(config["extension"]).to_value("deg")

    if extension > 1:
        logger.warning(
            "Extension comparable to the size of the MAGIC camera cannot be simulated"
        )
        is_valid = False
    if (config["n_off_regions"] <= 0) or (config["n_off_regions"] > 7):
        logger.warning("Number of OFF estimation regions must be in the range 1-7")
        is_valid = False
    if (extension > 0.5) and (config["n_off_regions"] > 1):
        logger.warning(
            "For large source extensions 1 OFF estimation"
            " region (n_off_regions) should be used."
        )
        is_valid = False
    if config["sum_trigger"] and (config["zenith_performance"] == "mid"):
        logger.warning("This functionality has not been yet implemented.")
        # raise NotImplementedError("This functionality has not been yet implemented.")
        is_valid = False
    if config["offset_degradation_factor"] > 1.00001:
        logger.warning(
            "No cheating! The performance degradation (%f) should not be larger than 1",
            config["offset_degradation_factor"],
        )
        is_valid = False
    if config["pulsar_mode"]["enable"]:
        if config["pulsar_on_range"] <= 0 or config["pulsar_on_range"] >= 1:
            logger.warning(
                "Pulsar mode ON phase range is %f, and it should be in range (0,1)",
                config["pulsar_on_range"],
            )
            is_valid = False
        if config["pulsar_off_range"] <= 0 or config["pulsar_off_range"] >= 1:
            logger.warning(
                "Pulsar mode OFF phase range is %f, and it should be in range (0,1)",
                config["pulsar_off_range"],
            )
            is_valid = False
        if config["redshift"] > 0:
            logger.warning(
                "Do you really want to observe a pulsar at redshift of %f ??",
                config["redshift"],
            )
            is_valid = False
    return is_valid


def initialize_model(config):
    """
    Import and initialize a spectral model.

    Parameters
    ----------
    config : `dict`
        Configuration data in form of a Python dictionary.

    Returns
    -------
    initialized_model : `~gammapy.modeling.models.SpectralModel`
        Initialized instance of a spectral model.
    """
    model_name = config["assumed_model"]["name"]
    module_name = ".".join(model_name.split(".")[:-1])
    class_name = model_name.split(".")[-1]
    module = importlib.import_module(module_name)
    model = getattr(module, class_name)
    model_parameters = config["assumed_model"]["parameters"]

    if (
        class_name == "LogParabolaSpectralModel"
        and config["assumed_model"]["from_log10"]
    ):
        initialized_model = model.from_log10(**model_parameters)
    else:
        initialized_model = model(**model_parameters)
    return initialized_model


@u.quantity_input(energy=u.TeV)
def observed_flux(energy, redshift, flux_int, ebl_file_path):
    """
    Get the attenuated flux of the source.

    Parameters
    ----------
    energy : `~astropy.units.Quantity`
        Array of energy values.
    redshift : `float`
        Redshift of the source.
    flux_int : `~astropy.units.Quantity`
        TBD.

    Returns
    -------
    valid : `bool`
        Initialized instance of a spectral model.
    """
    ebl_redshift, ebl_energy, ebl_taus = load_ebl(ebl_file_path)
    ftau = interpolate.interp2d(
        ebl_redshift, ebl_energy.to_value(energy.unit), ebl_taus, kind="cubic"
    )
    tau = ftau(redshift, energy.value)
    atten = np.exp(-tau)
    return flux_int * atten


# should be able to do it using gammapy
@u.quantity_input(energy=u.TeV, flux=u.Unit("TeV^-1 cm^-2 s^-1"))
def get_sed(energy, flux):
    """
    Compute the Spectral Energy Distribution (SED) from the source flux.

    Parameters
    ----------
    energy : `~astropy.units.Quantity`
        Energy values.
    flux : `~astropy.units.Quantity`
        Source energy flux values.

    Returns
    -------
    sed : `~astropy.units.Quantity`
        Spectral Energy Distribution.
    """
    sed = energy * energy * flux
    return sed


def significance_li_ma(n_on, n_off, alpha, mu_sig=None):
    """
    Get the Li & Ma significance.

    This is equivalent to eq.17 of [1]_.

    Parameters
    ----------
    n_on : `int`
        Measured counts in ON region.
    n_off : `int`
        Measured counts in OFF region.
    alpha : `float`
        Acceptance ratio of ON and OFF measurements.
    mu_sig : `float`
        Expected signal counts in ON region.

    Returns
    -------
    sqrt_ts : `float``
        Significance as the square root of the Test Statistic.

    Notes
    -----
    The implementation uses `gammapy.stats.WStatCountsStatistic`
    and takes the square root of the Test Statistic.

    References
    ----------
    .. [1] Li, T.-P. & Ma, Y.-Q., ApJ, 1983, 272, 317, 10.1086/161295.
    """
    statistics = WStatCountsStatistic(n_on, n_off, alpha, mu_sig)
    sqrt_ts = statistics.sqrt_ts
    return sqrt_ts


def prepare_data(config):
    """
    Extract the performance data.

    Parameters
    ----------
    config : `dict`
        Configuration data in form of a Python dictionary.

    Returns
    -------
    energy_bins : `~astropy.units.Quantity`
        Values of the energy bin edges.
    gamma_rate : `~astropy.units.Quantity`
        Rate of gamma-ray events from performance data.
    background_rate : `~astropy.units.Quantity`
        Rate of background events from performance data.
    """

    performance_data = {
        "low": LOW_ZENITH_PERFORMANCE,
        "mid": MID_ZENITH_PERFORMANCE,
    }

    if config["sum_trigger"] and config["zenith_performance"] == "low":
        message = "Low zenith performance with the SUM trigger"
        " is not currently available."
        logger.critical(message)
        raise NotImplementedError(message)

    min_energy = performance_data[config["zenith_performance"]]["min_energy"]
    max_energy = performance_data[config["zenith_performance"]]["max_energy"]
    energy_bins = np.append(min_energy.value, max_energy[-1].value) * min_energy.unit
    gamma_rate = performance_data[config["zenith_performance"]]["gamma_rate"]
    background_rate = performance_data[config["zenith_performance"]]["background_rate"]

    return energy_bins, gamma_rate, background_rate


@u.quantity_input(observation_time=u.h)
def source_detection(sigmas, observation_time):
    """
    Determine if the source can be detected.

    Parameters
    ----------
    sigmas : `list` of `float`
        Values of the significance for each energy bin.
    observation_time : `float`
        Observation time.

    Returns
    -------
    combined_significance : `float`
        Combined significance.
    """

    time = observation_time.to("h")

    if len(sigmas) > 0:
        combined_significance = sum(sigmas) / np.sqrt(len(sigmas))

        print(
            f"Combined significance (using the {len(sigmas):d} data points"
            f" shown in the SED) = {combined_significance:.1f}"
        )
    else:
        print(f"The source will not be detected in {time}.")

    if combined_significance < 4:
        print(f"The source probably will not be detected in {time}.")
    elif combined_significance < 6:
        print(f"The source probably might be detected in {time}.")
    else:
        print(f"The source probably will be detected in {time}.")

    return combined_significance


def calculate(energy_bins, gamma_rate, background_rate, config, assumed_spectrum):
    """
    Produce the necessary information to build an SED.

    Parameters
    ----------
    energy_bins : `~astropy.units.Quantity`
        Values of the energy bin edges.
    gamma_rate : `~astropy.units.Quantity`
        Rate of gamma-ray events from performance data.
    background_rate : `~astropy.units.Quantity`
        Rate of background events from performance data.
    config : `dict`
        Configuration data as a Python disctionary.
    assumed_spectrum : `~gammapy.modeling.models.SpectralModel`
        Assumed spectral model of the source.

    Returns
    -------
    en : `list` of `float`
        Energy values.
    sed : `list` of `float`
        Values for the Spectral energy distribution.
    dsed : `list` of `float`
        Error on Spectral energy distribution.
    sigmas : `list` of `float`
        Significance for each energy bin.
    detected : `list` of `bool`
        True if source is detected in energy bin,
        False otherwise.
    """

    n_off_regions = config["n_off_regions"]
    redshift = config["redshift"]
    observation_time_min = u.Quantity(config["observation_time"]).to("min")
    pulsar_mode_config = config["pulsar_mode"]
    pulsar_mode = pulsar_mode_config["enable"]
    pulsar_on_range = pulsar_mode_config["pulsar_on_range"]
    pulsar_off_range = pulsar_mode_config["pulsar_off_range"]
    min_num_events = config["min_number_events"]
    min_SBR = config["min_ratio_signal_background"]
    min_error = config["plotting_options"]["min_error"]
    psf = u.Quantity(config["PSF"])
    extension = u.Quantity(config["extension"])

    en = []
    sed = []
    dsed = []
    sigmas = []
    detected = []
    print(f"There are {len(energy_bins)} energy bins")
    for i, e1, e2 in zip(range(len(energy_bins)), energy_bins, energy_bins[1:]):
        # integral_crab_spectrum, _ = quad(
        #     lambda x: crab_nebula_spectrum()(x * e1.unit).value, e1.value, e2.value
        # )
        integral_crab_spectrum = crab_nebula_spectrum().integral(e1, e2)
        if redshift > 0:
            integral_assumed_spectrum, _ = quad(
                lambda x: observed_flux(x, redshift, assumed_spectrum(x)), e1, e2
            )
        else:
            # integral_assumed_spectrum, _ = quad(
            #     lambda x: assumed_spectrum(x * e1.unit).value, e1.value, e2.value
            # )
            integral_assumed_spectrum = assumed_spectrum.integral(e1, e2)
        n_off_events = background_rate[i].to("1/min") * observation_time_min
        n_off_events *= (psf**2 + extension**2) / (
            psf**2
        )  # larger integration cut due to extension
        n_off_events_error = np.sqrt(n_off_events / n_off_regions)
        n_excess = (
            gamma_rate[i].to("1/min")
            * observation_time_min
            * integral_assumed_spectrum
            / integral_crab_spectrum
        )

        error_n_excess = np.sqrt(n_excess + n_off_events + n_off_events_error**2)

        n_off_events_on_phase = 0
        if pulsar_mode:
            n_off_events_on_phase = n_off_events * pulsar_on_range
            n_off_events *= pulsar_off_range  # number of bgd events in OFF phase
            n_off_events_error = (
                np.sqrt(n_off_events) * pulsar_on_range / pulsar_off_range
            )  # ignoring numoff for pulsars and scaling for the phase difference
            error_n_excess = np.sqrt(
                n_excess + n_off_events_on_phase + n_off_events_error**2
            )

        # for tiny excesses (1.e-7 events) the function below can
        # have numerical problems, and either way sigma should be 0 then

        sigma = 0
        if n_excess > 0.01:
            if pulsar_mode:
                sigma = significance_li_ma(
                    n_excess + n_off_events_on_phase,
                    n_off_events,
                    pulsar_on_range / pulsar_off_range,
                )
                n_off_events = n_off_events_on_phase  # needed later for SBR
            else:
                sigma = significance_li_ma(
                    n_excess + n_off_events,
                    n_off_events * n_off_regions,
                    1.0 / n_off_regions,
                )

        detect = False  # assume no detection

        if pulsar_mode:
            if (sigma >= 5.0) and (n_excess > min_num_events):
                detect = True
        else:
            if (
                sigma >= 5.0
                and n_excess / n_off_events > min_SBR
                and n_excess > min_num_events
            ):
                detect = True
        print(
            "{0:.1f}-{1:.1f} GeV: exc. = {2:.1f}+-{3:.1f} ev., SBR={4:.2f}%, sigma = {5:.1f}".format(
                energy_bins[i],
                energy_bins[i + 1],
                n_excess,
                error_n_excess,
                100.0 * n_excess / n_off_events,
                sigma,
            ),
            " DETECTION" if detect else "",
        )
        logger.debug(
            "n_excess = %f, min_error = %f, error_n_excess = %f",
            n_excess,
            min_error,
            error_n_excess,
        )
        if n_excess > min_error * error_n_excess:
            tmpen = np.sqrt(energy_bins[i] * energy_bins[i + 1])
            en.append(tmpen.value)
            if redshift > 0:
                tmpsed = get_sed(
                    tmpen, observed_flux(tmpen, redshift, assumed_spectrum(tmpen))
                )
            else:
                tmpsed = get_sed(tmpen, assumed_spectrum(tmpen))

            sed.append(float(tmpsed.to_value(tmpsed.unit)))
            dsed.append(float(tmpsed.to_value(tmpsed.unit) * error_n_excess / n_excess))
            sigmas.append(sigma)
            detected.append(detect)
    return (
        en * energy_bins.unit,
        sed * tmpsed.unit,
        dsed * tmpsed.unit,
        sigmas,
        detected,
    )
