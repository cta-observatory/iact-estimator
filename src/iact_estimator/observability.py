"""Module that stores all functions related to observability estimation from the given location."""

from datetime import datetime

from astropy.time import Time

from astroplan import (
    AtNightConstraint,
    MoonSeparationConstraint,
    AltitudeConstraint,
    MoonIlluminationConstraint,
)
import astropy.units as u

from astroplan import (
    is_observable,
    months_observable,
    is_event_observable,
    time_grid_from_range,
)


__all__ = ["define_constraints"]


def get_days_in_this_year():
    """Get the number of days in the current year.

    Accounts for leap years."""
    this_year = datetime.today().year
    start_of_year = Time(f"{this_year}-01-01")
    end_of_year = Time(f"{this_year}-12-31")

    # Add 1 to include the last day
    num_days = int((end_of_year - start_of_year).to_value("day") + 1)
    return num_days * u.day


def define_constraints(config):
    constraints_config = config["observation"]["constraints"]

    moon_separation = constraints_config["moon_separation"]
    moon_illumination = constraints_config["moon_illumination_fraction"]
    zenith = constraints_config["zenith"]

    constraints = [
        AtNightConstraint(u.Quantity(constraints_config["max_solar_altitude"])),
        MoonSeparationConstraint(
            min=u.Quantity(moon_separation["min"]),
            max=u.Quantity(moon_separation["max"]) if moon_separation["max"] else None,
            ephemeris=moon_separation["ephemeris"],
        ),
        MoonIlluminationConstraint(
            min=u.Quantity(moon_illumination["min"])
            if moon_illumination["min"]
            else None,
            max=u.Quantity(moon_illumination["max"]),
            ephemeris=moon_separation["ephemeris"],
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
    """Check is the target source is observable and which months are best.

    Parameters
    ----------
    constraints : `~astroplan.constraints.Constraint` or list(~astroplan.constraints.Constraint)`
        One or more observational constraints.
    observer : `~astroplan.Observer`
        Instance of an _astroplan_ observer.
    target_source : `~astroplan.FixedTarget`
        Instance of an _astroplan_ fixed target source,
    time_range : list(~astropy.time.Time), shape=(2,)
        Start and stop observation time.
    time_grid_resolution : `~astropy.units.Quantity`
        Unit of time to test observability along the provided time range.

    Returns
    -------
    ever_observable : bool
        True if the target source is observable within at least one time unit.
    best_months : set
        Set of months where the target source is more suited for observation.
    """
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


def get_total_available_time(
    target_source, observer, constraints, time_range, time_resolution
):
    """Get available observation time.

    Args:
        target_source (astroplan.FixedTarget)
        observer (astroplan.Observer)
        constraints (astroplan.Constraint or list(astroplan.Constraint)):
            Set of constraining limits for the observation.
        time_range (list(astropy.time.Time)): [start, stop] time stamps
        time_resolution (u.Quantity): time step over which to evaluate the constraints.

    Returns:
        total_available_time (u.Quantity): in hours.
    """
    time_grid = time_grid_from_range(time_range, time_resolution=time_resolution)

    is_observable_at_time = is_event_observable(
        constraints, observer, target_source, times=time_grid
    )

    total_observable_time = is_observable_at_time.sum()

    total_available_time = total_observable_time * time_resolution.to("h")

    return total_available_time
