"""Plotting functions related to source observability from a location."""

import logging
from pathlib import Path

import astropy.units as u
from astroplan import FixedTarget
from astroplan.plots import plot_sky_24hr, plot_altitude
import matplotlib.pyplot as plt

from ..core import get_horizon_stereo_profile
from iact_estimator import HORIZON_PROFILE_M1, HORIZON_PROFILE_M2

logger = logging.getLogger(__name__)

__all__ = ["plot_transit", "plot_altitude_airmass"]


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
