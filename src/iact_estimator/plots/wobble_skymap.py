"""Plotting API for skymap with planned wobbles."""

from pathlib import Path

import astropy.units as u
import matplotlib.colors as matplotlib_colors
import matplotlib.pyplot as plt
import numpy as np
from starplot.styles import PolygonStyle
from starplot import MapPlot, Projection, Star, DSO
from starplot.styles import PlotStyle, extensions
from starplot.data.stars import STAR_NAMES
from starplot.data.bayer import hip

__all__ = ["calculate_square_coordinates", "load_wobbles", "plot_skymap_with_wobbles"]


def calculate_square_coordinates(ra0, dec0, delta_dec):
    """
    Calculate the required Right Ascension (RA) range to maintain a square plot.

    Parameters:
    ra0 (float): Central Right Ascension in degrees.
    dec0 (float): Central Declination in degrees.
    delta_dec (float): Size of the square in Declination degrees.

    Returns:
    tuple: (RA_min, RA_max, Dec_min, Dec_max)
    """

    dec0_rad = np.radians(dec0)

    delta_ra = delta_dec / np.cos(dec0_rad)

    ra_min = ra0 - delta_ra
    ra_max = ra0 + delta_ra
    dec_min = dec0 - delta_dec
    dec_max = dec0 + delta_dec

    return ra_min, ra_max, dec_min, dec_max


def load_wobbles(wobbles_config):
    """Load estimated wobble information from configuration.

    Parameters
    ----------
    wobbles_config : dict
        Dictionary with 'fov_offsets' and 'position_angles'
        as keys and a list of floats as values.

    Returns
    -------
    offsets : array-like
        Array of FoV offset angles in degrees.
    angles : array-like
        Array or position angles in degrees.
    """
    cfg_angles = (
        [wobbles_config["position_angles"]]
        if not isinstance(wobbles_config["position_angles"], list)
        else wobbles_config["position_angles"]
    )
    angles = (
        np.asarray(cfg_angles) * u.deg if len(cfg_angles) > 1 else cfg_angles * u.deg
    )

    cfg_offsets = (
        [wobbles_config["fov_offsets"]]
        if not isinstance(wobbles_config["fov_offsets"], list)
        else wobbles_config["fov_offsets"]
    )
    offsets = (
        (np.repeat(cfg_offsets, len(cfg_angles)) * u.deg)
        if len(cfg_offsets) == 1
        else cfg_offsets * u.deg
    )

    return offsets, angles


def plot_skymap_with_wobbles(
    target_source,
    observer,
    instrument_field_of_view,
    wobble_angles,
    wobble_offsets,
    config,
    savefig=True,
    output_path=None,
):
    target_coordinates = target_source.coord

    plotting_options = config["wobble_skymap_plot_options"]

    style = PlotStyle().extend(
        getattr(extensions, plotting_options["map_color_scheme"]),
        extensions.MAP,
        {
            "legend": plotting_options["legend"],
        },
    )

    ra_min, ra_max, dec_min, dec_max = calculate_square_coordinates(
        target_coordinates.ra, target_coordinates.dec, 1.5 * instrument_field_of_view
    )

    if -70 < target_coordinates.dec.deg < 70:
        projection = Projection.MERCATOR
    else:
        if observer.latitude < 0 and target_coordinates.dec.deg < -70:
            projection = Projection.STEREO_SOUTH
        if observer.latitude > 0 and target_coordinates.dec.deg > 70:
            projection = Projection.STEREO_NORTH

    p = MapPlot(
        projection=projection,  # specify a non-perspective projection
        ra_min=ra_min.to_value("hourangle"),
        ra_max=ra_max.to_value("hourangle"),
        dec_min=dec_min.to_value("deg"),
        dec_max=dec_max.to_value("deg"),
        style=style,
        resolution=1500,
    )

    p.gridlines()

    DANGER_MAGNITUDE = plotting_options["magnitude"]["danger"]
    MAGNITUDE_LIMIT = plotting_options["magnitude"]["max"]

    # first we define the callable:
    def color_by_mag(star: Star) -> str:
        if star.magnitude <= DANGER_MAGNITUDE:
            return "red"
        else:
            return "black"

    def bright_star_labels() -> dict:
        return {
            k: f"{v} {Star.get(name=v).magnitude}"
            if Star.get(name=v).magnitude <= DANGER_MAGNITUDE
            else v
            for k, v in STAR_NAMES.items()
        }

    def star_label_with_magnitude(star: Star) -> str:
        try:
            star_label = hip[star.hip]

        except KeyError:
            star_label = star.tyc if star.tyc is not None else star.hip

        label = f"{star_label} ({star.magnitude})"

        return label

    p.stars(
        where=[
            (Star.magnitude > DANGER_MAGNITUDE) & (Star.magnitude < MAGNITUDE_LIMIT)
        ],
    )

    p.stars(
        where=[Star.magnitude < DANGER_MAGNITUDE],
        where_labels=[Star.magnitude < DANGER_MAGNITUDE],
        bayer_labels=True,
        label_fn=star_label_with_magnitude,
        color_fn=color_by_mag,
    )

    p.galaxies(
        where=[DSO.magnitude < MAGNITUDE_LIMIT],
    )

    p.nebula(
        where=[
            DSO.magnitude < MAGNITUDE_LIMIT,
        ],
        true_size=False,
        label_fn=lambda d: d.ic,
    )
    p.open_clusters(
        where=[
            DSO.magnitude < MAGNITUDE_LIMIT,
        ],
    )
    p.milky_way()
    p.globular_clusters()

    p.marker(
        ra=target_coordinates.ra.hour,
        dec=target_coordinates.dec.deg,
        label=None,
        legend_label="target source",
        style=plotting_options["target_source"],
    )

    color_map = plt.get_cmap(plotting_options["wobbles_colormap"])

    for angle, offset, color in zip(wobble_angles, wobble_offsets, color_map.colors):
        color = matplotlib_colors.to_hex(color)

        # Transform IACT convention to astropy
        position_angle = (90 * u.deg - angle) % (360 * u.deg)

        wobble_center = target_coordinates.directional_offset_by(position_angle, offset)
        wobble_ra = wobble_center.ra.to_value("hourangle")
        wobble_dec = wobble_center.dec.to_value("deg")

        p.circle(
            center=(wobble_ra, wobble_dec),
            style=PolygonStyle(
                fill_color=None,
                edge_width=1,
                edge_color=color,
                alpha=1,
            ),
            radius_degrees=instrument_field_of_view.to_value("deg") / 2.0,
        )

        p.marker(
            ra=wobble_ra,
            dec=wobble_dec,
            label=None,
            legend_label=f"W{offset}+{angle}",
            style={
                "marker": {
                    "size": 15,
                    "symbol": "point",
                    "fill": "full",
                    "color": color,
                    "alpha": 0.4,
                },
                "label": {
                    "font_size": 15,
                    "font_weight": "bold",
                    "font_color": color,
                    "font_alpha": 0.8,
                },
            },
        )

    p.legend()

    if savefig:
        output_path = output_path if output_path is not None else Path.cwd()
        p.export(
            f"{target_source.name}_skymap_with_wobbles.{config['wobble_skymap_plot_options']['export']['format']}",
            **config["wobble_skymap_plot_options"]["export"],
        )
