"""Tests for the observability module."""

import numpy as np
import astropy.units as u
from astropy.time import Time
from astroplan import AtNightConstraint, MoonSeparationConstraint, AltitudeConstraint


def test_get_days_in_this_year():
    """Test getting the number of days in the current year."""
    from ..observability import get_days_in_this_year

    num_days = get_days_in_this_year()

    # Should be either 365 or 366 days
    assert num_days in [365 * u.day, 366 * u.day]
    assert num_days.unit == u.day


def test_define_constraints(sample_config):
    """Test defining observational constraints from config."""
    from ..observability import define_constraints

    constraints = define_constraints(sample_config)

    # Check that we have the expected number of constraints
    assert len(constraints) == 4

    # Check constraint types
    assert any(isinstance(c, AtNightConstraint) for c in constraints)
    assert any(isinstance(c, MoonSeparationConstraint) for c in constraints)
    assert any(isinstance(c, AltitudeConstraint) for c in constraints)


def test_define_constraints_custom_values(sample_config):
    """Test defining constraints with custom values."""
    from ..observability import define_constraints

    # Modify config
    sample_config["observation"]["constraints"]["zenith"]["max"] = "45 deg"
    sample_config["observation"]["constraints"]["moon_separation"]["min"] = "45 deg"

    constraints = define_constraints(sample_config)

    # Find altitude constraint and check values
    altitude_constraints = [c for c in constraints if isinstance(c, AltitudeConstraint)]
    assert len(altitude_constraints) == 1

    # Altitude = 90 - zenith, so zenith max of 45 deg should give altitude min of 45 deg
    altitude_constraint = altitude_constraints[0]
    assert altitude_constraint.min == 45 * u.deg


def test_define_constraints_no_moon_max():
    """Test defining constraints with no maximum moon separation."""
    from ..observability import define_constraints

    config = {
        "observation": {
            "constraints": {
                "max_solar_altitude": "-18 deg",
                "moon_separation": {
                    "min": "30 deg",
                    "max": None,  # No maximum
                    "ephemeris": "builtin",
                },
                "moon_illumination_fraction": {
                    "min": None,
                    "max": 1.0,
                },
                "zenith": {
                    "min": "0 deg",
                    "max": "30 deg",
                },
            },
        },
    }

    constraints = define_constraints(config)
    assert len(constraints) == 4

    # Find moon separation constraint
    moon_constraints = [
        c for c in constraints if isinstance(c, MoonSeparationConstraint)
    ]
    assert len(moon_constraints) == 1


def test_check_observability(sample_target_source, sample_observer, sample_config):
    """Test checking observability of a target."""
    from ..observability import check_observability, define_constraints

    constraints = define_constraints(sample_config)

    time_range = [Time("2024-06-01"), Time("2024-08-31")]
    time_grid_resolution = 1 * u.hour

    ever_observable, best_months = check_observability(
        constraints,
        sample_observer,
        [sample_target_source],
        time_range,
        time_grid_resolution,
    )

    # Check return types - is_observable returns array for list of targets
    assert isinstance(ever_observable[0], (bool, np.bool_))
    # months_observable returns a list of sets when given a list of targets
    assert isinstance(best_months, list)
    assert isinstance(best_months[0], set)


def test_get_total_available_time(
    sample_target_source, sample_observer, sample_config, sample_time_range
):
    """Test computing total available observation time."""
    from ..observability import get_total_available_time, define_constraints

    constraints = define_constraints(sample_config)
    time_resolution = 30 * u.min

    total_time = get_total_available_time(
        sample_target_source,
        sample_observer,
        constraints,
        sample_time_range,
        time_resolution,
    )

    # Check that we get a valid time quantity
    assert total_time.unit.is_equivalent(u.hour)
    assert total_time.value >= 0

    # For a year, should have some observable time
    # (exact value depends on constraints and target)
    assert total_time < 365 * 24 * u.hour  # Can't exceed total hours in a year


def test_get_total_available_time_different_resolutions(
    sample_target_source, sample_observer, sample_config
):
    """Test that different time resolutions give consistent results."""
    from ..observability import get_total_available_time, define_constraints

    constraints = define_constraints(sample_config)

    # Short time range for faster testing
    time_range = [Time("2024-06-01"), Time("2024-06-30")]

    # Test with different resolutions
    time_coarse = get_total_available_time(
        sample_target_source, sample_observer, constraints, time_range, 60 * u.min
    )

    time_fine = get_total_available_time(
        sample_target_source, sample_observer, constraints, time_range, 30 * u.min
    )

    # Results should be similar (within reasonable tolerance due to discretization)
    # Fine resolution should be more accurate
    assert abs(time_coarse - time_fine) < 10 * u.hour


def test_get_total_available_time_unreachable_target(sample_observer, sample_config):
    """Test total available time for an unreachable target."""
    from astroplan import FixedTarget
    from astropy.coordinates import SkyCoord
    from ..observability import get_total_available_time, define_constraints

    # Create a target in the southern hemisphere (unreachable from La Palma)
    south_pole_coords = SkyCoord(ra=0 * u.deg, dec=-89 * u.deg, frame="icrs")
    south_pole_target = FixedTarget(coord=south_pole_coords, name="Near South Pole")

    constraints = define_constraints(sample_config)
    time_range = [Time("2024-06-01"), Time("2024-06-30")]

    total_time = get_total_available_time(
        south_pole_target, sample_observer, constraints, time_range, 60 * u.min
    )

    # Should have very little or no observable time
    assert total_time.value < 1.0  # Less than 1 hour


def test_check_observability_always_observable_target(sample_observer, sample_config):
    """Test observability for a target that should be observable during the year."""
    from astroplan import FixedTarget
    from astropy.coordinates import SkyCoord
    from ..observability import check_observability, define_constraints

    # Use Crab Nebula which is known to be observable from La Palma
    # and is a common target for IACTs
    crab_coords = SkyCoord(ra=83.63 * u.deg, dec=22.01 * u.deg, frame="icrs")
    crab_target = FixedTarget(coord=crab_coords, name="Crab Nebula")

    constraints = define_constraints(sample_config)

    time_range = [Time("2024-01-01"), Time("2024-12-31")]
    time_grid_resolution = 24 * u.hour  # Coarse resolution for speed

    ever_observable, best_months = check_observability(
        constraints, sample_observer, [crab_target], time_range, time_grid_resolution
    )

    # Crab should be observable at some point during the year
    # ever_observable[0] is numpy.bool_, not Python bool, so use truthiness
    assert ever_observable[0]
    assert len(best_months) > 0
