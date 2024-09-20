"""Module that stores all functions related to observability estimation from the given location."""

from astroplan import AtNightConstraint, MoonSeparationConstraint, AltitudeConstraint
import astropy.units as u

from astroplan import is_observable, months_observable


__all__ = ["define_constraints"]


def define_constraints(config):
    constraints_config = config["observation"]["constraints"]

    moon = constraints_config["moon_separation"]
    zenith = constraints_config["zenith"]

    constraints = [
        AtNightConstraint(u.Quantity(constraints_config["max_solar_altitude"])),
        MoonSeparationConstraint(
            min=u.Quantity(moon["min"]),
            max=u.Quantity(moon["max"]),
            ephemeris=moon["ephemeris"],
        ),
        AltitudeConstraint(
            min=90 * u.deg - u.Quantity(zenith["max"]),
            max=90 * u.deg - u.Quantity(zenith["min"]),
        ),
    ]

    return constraints


def check_observability(
    constraints, observer, target_source, time_range, time_grid_resolution
):
    # Are targets *ever* observable in the time range?
    ever_observable = is_observable(
        constraints,
        observer,
        target_source,
        time_range=time_range,
        time_grid_resolution=time_grid_resolution,
    )

    # During what months are the targets ever observable?
    best_months = months_observable(
        constraints,
        observer,
        target_source,
        time_range,
        time_grid_resolution=time_grid_resolution,
    )

    return ever_observable, best_months
