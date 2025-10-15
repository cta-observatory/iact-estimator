"""Plotting functions related to other wavelengths."""

import astropy.units as u
import matplotlib.pyplot as plt

from astroplan.plots import plot_finder_image


@u.quantity_input(fov_radius=u.arcmin)
def plot_from_skyview_survey(
    target_source,
    survey_name="DSS",
    fov_radius=10 * u.arcmin,
    log=False,
    ax=None,
    reticle=False,
    style_kwargs=None,
    reticle_style_kwargs=None,
):
    """Plot a survey image from the SkyView service centered on ``target``.

    Parameters
    ----------
    target : `~astroplan.FixedTarget`, `~astropy.coordinates.SkyCoord`
        Coordinates of celestial object

    survey : string
        Name of survey to retrieve image from. For dictionary of
        available surveys, use
        ``from astroquery.skyview import SkyView; SkyView.list_surveys()``.
        Defaults to ``'DSS'``, the Digital Sky Survey.

    fov_radius : `~astropy.units.Quantity`
        Radius of field of view of retrieved image. Defaults to 10 arcmin.

    log : bool, optional
        Take the natural logarithm of the FITS image if `True`.
        False by default.

    ax : `~matplotlib.axes.Axes` or None, optional.
        The `~matplotlib.axes.Axes` object to be drawn on.
        If None, uses the current `~matplotlib.axes.Axes`.

    reticle : bool, optional
        Draw reticle on the center of the FOV if `True`. Default is `False`.

    style_kwargs : dict or `None`, optional.
        A dictionary of keywords passed into `~matplotlib.pyplot.imshow`
        to set plotting styles.

    reticle_style_kwargs : dict or `None`, optional
        A dictionary of keywords passed into `~matplotlib.pyplot.axvline` and
        `~matplotlib.pyplot.axhline` to set reticle style.

    Returns
    -------
    ax : `~matplotlib.axes.Axes`
        Matplotlib axes with survey image centered on ``target``

    hdu : `~astropy.io.fits.PrimaryHDU`
        FITS HDU of the retrieved image

    Notes
    -----
    This is wrapper function around `astroplan.plots.plot_finder_image()`.

    """

    ax = plt.gca() if ax is None else ax

    ax = plot_finder_image(
        target_source,
        survey=survey_name,
        fov_radius=fov_radius,
        log=log,
        ax=ax,
        reticle=reticle,
        style_kwargs=style_kwargs,
        reticle_style_kwargs=reticle_style_kwargs,
    )

    return ax
