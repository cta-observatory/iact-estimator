"""Plotting functions related to source observability from a location."""

import logging
import warnings

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from datetime import timedelta
from pathlib import Path

from astroplan import FixedTarget
from astroplan.plots import plot_sky_24hr, plot_altitude
from astroplan.utils import time_grid_from_range
from astropy.coordinates.errors import NonRotationTransformationWarning
from astropy.time import Time

from ..observability import get_total_available_time

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


def plot_observability_constraints_grid(
    source_name,
    config,
    observer,
    target,
    start_time,
    end_time,
    time_resolution,
    constraints,
    ax=None,
    savefig=True,
    output_path=None,
):
    """
    Plot a grid representing the observability of a target based on a set of constraints
    over a specified time range.

    Parameters
    ----------
    source_name : str
        The name of the source being observed.
    config : dict
        Configuration dictionary containing plotting options, such as file format.
    observer : `~astroplan.Observer`
        The observer object that provides the position and time of observation.
    target : `~astroplanb.FixedTarget`
        The target object representing the source of interest.
    start_time : datetime
        The starting time of the observation period.
    end_time : datetime
        The ending time of the observation period.
    time_resolution : `~astropy.units.Quantity`
        The time resolution (step size) for the grid.
    constraints : `~astroplan.constraints.Constraint` or list(~astroplan.constraints.Constraint)`
        A list of constraint functions that evaluate the observability of the target
        at given times for the observer and target.
    ax : `~matplotlib.axes.Axes`, optional
        The axis on which to plot the grid. If None, a new figure and axis are created.
    savefig : bool, optional
        Whether to save the generated plot as a file. Default is True.
    output_path : str or Path, optional
        The directory where the plot should be saved. If None, saves to the current
        working directory.

    Returns
    -------
    ax : `~matplotlib.axes.Axes`
        The axis with the generated observability grid plot.

    Notes
    -----
    The observability grid is created by evaluating each constraint function
    at every time step within the specified time range. The grid is plotted as
    an image, with time on the x-axis and constraints on the y-axis.
    """
    fig, ax = plt.subplots(layout="constrained")
    ax = plt.cga() if ax is None else ax

    # Create grid of times from ``start_time`` to ``end_time``
    # with resolution ``time_resolution``
    time_grid = time_grid_from_range(
        [start_time, end_time], time_resolution=time_resolution
    )

    observability_grid = np.zeros((len(constraints), len(time_grid)))

    constraint_labels = {}
    for i, constraint in enumerate(constraints):
        # Evaluate each constraint
        observability_grid[i, :] = constraint(observer, target, times=time_grid)
        constraint_labels[i] = f"{constraint.__class__.__name__.split('Constraint')[0]}"

    # Create plot showing observability of the target

    numcols = len(time_grid)
    numrows = len(constraint_labels)
    extent = [-0.5, numcols - 0.5, numrows - 0.5, -0.5]

    ax.imshow(observability_grid, extent=extent)

    ax.set_yticks(range(0, 3))
    ax.set_yticks(list(constraint_labels.keys()), constraint_labels.values())

    ax.set_xticks(range(len(time_grid)))
    ax.set_xticklabels([t.datetime.strftime("%H:%M") for t in time_grid], fontsize=7)

    ax.set_xticks(np.arange(extent[0], extent[1]), minor=True)
    ax.set_yticks(np.arange(extent[2], extent[3]), minor=True)

    ax.grid(which="minor", color="w", linestyle="-", linewidth=2)
    ax.tick_params(axis="x", which="minor", bottom="off")
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right")

    ax.tick_params(axis="y", which="minor", left="off")

    if savefig:
        output_path = output_path if output_path is not None else Path.cwd()
        fig.savefig(
            output_path
            / f"observability_grid_{source_name}.{config['plotting_options']['file_format']}",
        )

    return ax


def create_observability_heatmap(
    target_source,
    observer,
    constraints,
    start_date,
    end_date,
    time_resolution=1 * u.hour,
    cmap="YlGnBu",
    sns_plotting_context="paper",
    sns_axes_style="whitegrid",
    savefig=True,
    output_path=None,
    save_format="png",
):
    """Plot an annotated heatmap showing the amount of available hours per day.

    Parameters
    ==========
    target_source: `~astroplan.FixedTarget`
    observer: `~astroplan.Observer`
    constraints: `~astroplan.Constraint` or `list(~astroplan.Constraint)`
    start_date: `~astropy.time.Time`
    end_date: `~astropy.time.Time`
    time_resolution: `u.Quantity`, default=1h
    cmap: str, default="YlGnBu"
    sns_plotting_context: str, default="paper"
    sns_axes_style: str, default="darkgrid"
    savefig: bool, default=True
    output_path: `pathlib.Path`, default=None
        If unspecified the figure is saved in the current working directory.
    save_format: str, default="png"
    """
    date_range = pd.date_range(
        start=start_date.datetime, end=end_date.datetime, freq="D"
    )

    data = []

    for date in date_range:
        # Define time range for this particular day
        time_range = Time([date, date + timedelta(days=1)])

        # Get available observation hours for the day
        # see also https://github.com/astropy/astroplan/issues/598
        warnings.filterwarnings("ignore", category=NonRotationTransformationWarning)
        # see also https://github.com/astropy/astroplan/issues/132#issuecomment-225471962
        warnings.filterwarnings(
            "ignore", module="astropy.coordinates.builtin_frames.utils"
        )
        available_hours = get_total_available_time(
            target_source, observer, constraints, time_range, time_resolution
        ).value

        # Append data with month, day, and available hours
        data.append(
            {
                "year": date.year,
                "month": date.month,
                "day": date.day,
                "available_hours": available_hours,
            }
        )

    # Create DataFrame
    df = pd.DataFrame(data)

    if df["year"].nunique() > 1:
        df["month_label"] = (
            df["year"].astype(str) + "-" + df["month"].astype(str).str.zfill(2)
        )
    else:
        df["month_label"] = df["month"]

    heatmap_data = df.pivot(
        index="month_label", columns="day", values="available_hours"
    )

    with (
        sns.axes_style(sns_axes_style),
        sns.plotting_context(sns_plotting_context),
        sns.color_palette("deep"),
    ):
        # Dynamically set the figure size based on data dimensions
        fig_width = 1 + heatmap_data.shape[1] * 0.3  # 0.3 inch per day
        fig_height = 1 + heatmap_data.shape[0] * 0.5  # 0.5 inch per month
        fig, ax = plt.subplots(figsize=(fig_width, fig_height), layout="constrained")
        fmt = ".0f" if time_resolution == 1 * u.h else ".1f"
        sns.heatmap(
            heatmap_data,
            ax=ax,
            annot=True,
            fmt=fmt,
            cmap=cmap,
            annot_kws={"size": 12, "fontweight": 10},
            cbar_kws={"label": "Available Observation Hours"},
        )
        ax.set_xlabel("Day of the Month")
        ax.set_ylabel("Month")

    if savefig:
        output_path = output_path if output_path is not None else Path.cwd()
        fig.savefig(
            output_path / f"observability_heatmap_{target_source.name}.{save_format}",
        )
