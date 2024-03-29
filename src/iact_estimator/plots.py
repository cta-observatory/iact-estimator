"""Plotting functions."""

import logging
from pathlib import Path

import astropy.units as u
from astropy.visualization import quantity_support
from astroplan import FixedTarget
from astroplan.plots import plot_sky_24hr, plot_altitude
import matplotlib.pyplot as plt
import numpy as np

from .core import observed_flux, get_horizon_stereo_profile
from .spectral import crab_nebula_spectrum
from iact_estimator import HORIZON_PROFILE_M1, HORIZON_PROFILE_M2

__all__ = [
    "plot_spectrum",
    "plot_sed",
    "plot_transit",
    "plot_altitude_airmass",
    "plot_exposure",
    "plot_rates",
]

logger = logging.getLogger(__name__)


def plot_spectrum(
    config,
    energy_bounds,
    model,
    source_name,
    plotting_options,
    savefig=True,
    output_path=None,
    **kwargs,
):
    """
    Plot a spectrum from a model.

    Parameters
    ----------
    energy_bounds : `~astropy.units.Quantity`
        Plot energy bounds.
    model : `~gammapy.modeling.models.SpectralModel`
        Spectral model to plot.
    output_path : `str` or `pathlib.Path`
        Path to the output directory where to save the plot.
    source_name : `str`
        Name of the source.
    plotting_options : `dict`
        Dictionary of options related to plotting
        from the configuration file.
    **kwargs :
        Keyword arguments for `~matplotlib.pyplot.plot`.
    """
    fig, ax = plt.subplots(figsize=plotting_options["figure_size"])
    model.plot(energy_bounds, **kwargs)
    crab_nebula_spectrum().plot(
        energy_bounds, label="Crab Aleksic 2016", linewidth=5, alpha=0.5
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.grid(which="both", axis="both", visible=True)
    ax.set_xlabel(rf"Energy ({energy_bounds[0].unit})")
    ax.set_ylabel(
        rf"$\dfrac{{dN}}{{dE dA dt}}$ [ {model.amplitude.unit.to_string('latex', fraction=False)} ]"
    )
    plt.legend()

    if savefig:
        output_path = output_path if output_path is not None else Path.cwd()
        fig.savefig(
            output_path
            / f"{source_name}_spectrum.{config['plotting_options']['file_format']}",
            bbox_inches=plotting_options["bbox_inches"],
        )
        logger.debug("Plot has been successfully saved at %s", output_path)


def plot_sed(
    config,
    sigmas,
    sig_comb,
    source_name,
    assumed_spectrum,
    en,
    sed,
    dsed,
    detected,
    savefig=True,
    output_path=None,
    annotation_options={
        "rotation": 45,
        "xytext": (10, 10),
        "size": 15,
        "horizontalalignment": "left",
        "verticalalignment": "bottom",
    },
):
    """
    Plot the Spectral Energy distribution with significances.

    Parameters
    ----------
    config : dict
        Loaded configuration
    sigmas : array-like
    sig_comb : float
        Combined significance
    source_name : str
    assumed_spectrum : `~gammapy.modeling.models.SpectralModel`
    en :
    sed :
    dsed :
    detected :
    savefig : bool
    output_path : str or `~pathlib.Path`
    annotation_options : dict
        Options for `matplotlib.axes.Axes.annotate`.

    Notes
    -----
    Spectral points following the assumed spectrum are shown.
    Their error bars reflect the performance of the instrument for such a source.

    Significance for each bin is given. Bins without detection
    have gray numbers.

    At the top of the plot a simple number using
    information from all the shown bins is given to evaluate
    if the source can be detected, roughly if

    .. math::

        \dfrac{ \sum  sigma_{i} } { \sqrt{N} } \gtrsim 5
    """
    fig, ax = plt.subplots(figsize=config["plotting_options"]["figure_size"])
    ax.set_xscale("log")
    ax.set_yscale("log")

    min_energy = u.Quantity(config["plotting_options"]["min_energy"])
    max_energy = u.Quantity(config["plotting_options"]["max_energy"])
    energy_unit = u.Unit(config["plotting_options"]["energy_unit"])
    energy_flux_unit = u.Unit(config["plotting_options"]["energy_flux_unit"])
    redshift = config["redshift"]

    ax.set(
        xlabel=rf"Energy [ {energy_unit} ]",
        ylabel=rf"$E^{2}$ $\dfrac{{dN}}{{dE dA dt}}$ [ {energy_flux_unit.to_string('latex', fraction=False)} ]",
        title=rf"$\Sigma\sigma_i / \sqrt{{{len(sigmas)}}} = {{{round(sig_comb.value, 1)}}}$",
    )
    energy = np.logspace(
        np.log10(min_energy.to_value(energy_unit)),
        np.log10(max_energy.to_value(energy_unit)),
        50,
    ) * u.Unit(energy_unit)
    labeltext = rf"Expected SED ($T_{{obs}}$ = {config['observation_time']})"
    plt.plot(
        energy,
        energy * energy * crab_nebula_spectrum()(energy),
        linewidth=10,
        alpha=0.5,
        label="Crab (Aleksic et al 2016)",
    )
    if redshift > 0:
        attenuated_flux = [
            float(observed_flux(ee, redshift, assumed_spectrum(ee))) for ee in energy
        ]
        plt.plot(
            energy,
            energy * energy * attenuated_flux,
            label=f"{source_name} (Assumed, z={redshift:.2f})",
        )
    else:
        plt.plot(
            energy,
            energy * energy * assumed_spectrum(energy),
            label=f"{source_name} (Assumed)",
        )
    if len(en) > 0:
        ax.errorbar(en, sed, yerr=dsed, label=labeltext, color="0", fmt="o")
    ax.legend(loc="upper right")
    ax.grid(True, which="both", axis="both", ls="--", color="0.95")

    if config["plotting_options"]["draw_sigma"]:
        for i in range(len(sigmas)):
            col = "0" if detected[i] else "0.75"
            ax.annotate(
                f"{sigmas[i]:.1f}$\sigma$",
                (en[i], sed[i]),
                color=col,
                xycoords="data",
                textcoords="offset points",
                **annotation_options,
            )
    if savefig:
        output_path = output_path if output_path is not None else Path.cwd()
        fig.savefig(
            output_path
            / f"SED_{source_name}.{config['plotting_options']['file_format']}",
        )


def plot_transit(
    config,
    source_name,
    target_source,
    observer,
    time,
    merge_profiles=True,
    plot_crab=True,
    style_kwargs=None,
    savefig=True,
    output_path=None,
):
    """
    Plot the Spectral Energy distribution with significances.

    Parameters
    ----------
    config : dict
        Loaded configuration
    output_path : str or `~pathlib.Path`
    source_name : str
    target_source : `astroplan.Target`
    observer : `astroplan.Observer`
    time : `~astropy.time.Time`
        Datetime of the planned observation
    merge_profiles : bool, default=True
        If True plot the combined horizon profile
        from both telescopes.
    crab : bool, default=True
        If True plot the Crab together with
        the target source for comparison.
    style_kwargs : dict, default=None
        Dictionary of keywords passed into
        `~matplotlib.pyplot.scatter`
        to set plotting styles.
    """
    fig = plt.figure(figsize=config["plotting_options"]["figure_size"])
    ax = fig.add_subplot(projection="polar")
    ax.tick_params(pad=10)

    plot_sky_24hr(
        target_source,
        observer,
        time,
        delta=1 * u.h,
        ax=ax,
        style_kwargs=style_kwargs,
        north_to_east_ccw=True,
        grid=True,
        az_label_offset=0 * u.deg,
        center_time_style_kwargs=None,
    )
    if plot_crab:
        crab_nebula = FixedTarget.from_name("Crab Nebula")
        plot_sky_24hr(
            crab_nebula,
            observer,
            time,
            delta=1 * u.h,
            ax=ax,
            style_kwargs={"label": "Crab Nebula"},
            north_to_east_ccw=True,
            grid=True,
            az_label_offset=0 * u.deg,
            center_time_style_kwargs=None,
        )

    if merge_profiles:
        az, zd = get_horizon_stereo_profile(HORIZON_PROFILE_M1, HORIZON_PROFILE_M2)
        ax.plot(
            az.to_value("rad"),
            zd.to_value("deg"),
            linestyle="-",
            label="Horizon profile",
            alpha=0.5,
            color="#4daf4a",
        )
    else:
        ax.plot(
            HORIZON_PROFILE_M2["azimuth"].to("rad").value,
            HORIZON_PROFILE_M2["zenith"].value,
            linestyle="-",
            label="Horizon profile M2",
            alpha=0.5,
            color="#4daf4a",
        )
        ax.plot(
            HORIZON_PROFILE_M1["azimuth"].to("rad").value,
            HORIZON_PROFILE_M1["zenith"].value,
            linestyle="-",
            label="Horizon profile M1",
            alpha=0.5,
            color="#f781bf",
        )

    ax.legend(loc="center", bbox_to_anchor=(0.5, 0.8))
    if savefig:
        output_path = output_path if output_path is not None else Path.cwd()
        fig.savefig(
            output_path
            / f"{source_name}_skyplot.{config['plotting_options']['file_format']}",
            bbox_inches=config["plotting_options"]["bbox_inches"],
        )
        logger.debug("Plot has been successfully saved at %s", output_path)
    return ax


def plot_altitude_airmass(
    config,
    source_name,
    target_source,
    observer,
    time,
    ax=None,
    savefig=True,
    output_path=None,
    **kwargs,
):
    fig, ax = plt.subplots(figsize=config["plotting_options"]["figure_size"])

    plot_altitude(target_source, observer, time, ax=ax, **kwargs)
    ax.grid(False)

    if savefig:
        output_path = output_path if output_path is not None else Path.cwd()
        fig.savefig(
            output_path
            / f"{source_name}_altitude_airmass.{config['plotting_options']['file_format']}",
            bbox_inches=config["plotting_options"]["bbox_inches"],
        )
    return ax


def plot_exposure(data):
    _, ax = plt.subplots()
    reco_energy_edges = data.to_numpy()[1] * u.GeV
    reco_energy = 0.5 * (reco_energy_edges[1:] + reco_energy_edges[:-1])
    exposure = data.to_numpy()[0] * u.cm * u.s
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel(r"Exposure [cm$^{2}$ s]")
    ax.set_xlabel("Reconstructed energy [GeV]")
    ax.step(reco_energy, exposure, where="mid")
    ax.grid(which="both", axis="both")


def plot_rates(performance_data, title=None, ax=None):
    ax = plt.gca() if ax is None else ax
    if title:
        ax.set_title(title)
    with quantity_support():
        energy_edges = np.concatenate(
            (
                performance_data["min_energy"].value,
                np.array([performance_data["max_energy"][-1].value]),
            )
        )
        energy = 0.5 * (energy_edges[1:] + energy_edges[:-1])
        ax.step(
            energy * performance_data["max_energy"].unit,
            performance_data["gamma_rate"],
            label="gamma rate",
        )
        ax.step(
            energy * performance_data["max_energy"].unit,
            performance_data["background_rate"],
            label="bkg rate",
        )
        ax.legend()
        ax.set_xscale("log")
        ax.set_yscale("log")
    return ax
