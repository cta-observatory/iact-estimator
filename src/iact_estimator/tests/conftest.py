"""Shared pytest fixtures for the test suite."""

import numpy as np
import pytest
from astropy.table import QTable
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroplan import FixedTarget, Observer
from astropy.time import Time


@pytest.fixture
def sample_config():
    """Sample configuration dictionary for testing."""
    return {
        "target_source": {
            "name": "Test Source",
            "coordinates": {
                "ra_l": "05h34m31.94s",
                "dec_b": "+22d00m52.2s",
                "frame": "icrs",
                "force": False,
            },
            "assumed_model": {
                "name": "gammapy.modeling.models.LogParabolaSpectralModel",
                "from_log10": True,
                "parameters": {
                    "amplitude": "3e-11 TeV-1 cm-2 s-1",
                    "reference": "1 TeV",
                    "alpha": 2.5,
                    "beta": 0.2,
                },
            },
        },
        "extension": "0.0 deg",
        "n_off_regions": 3,
        "sum_trigger": False,
        "magic_lst1": False,
        "zenith_range": "low",
        "offset_degradation_factor": 1.0,
        "pulsar_mode": {
            "enable": False,
            "pulsar_on_range": 0.1,
            "pulsar_off_range": 0.9,
        },
        "pulsar_on_range": 0.1,
        "pulsar_off_range": 0.9,
        "redshift": 0.0,
        "observation": {
            "time": "50 h",
            "constraints": {
                "max_solar_altitude": "-18 deg",
                "moon_separation": {
                    "min": "30 deg",
                    "max": None,
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
        "min_number_events": 10,
        "min_ratio_signal_background": 0.05,
        "PSF": "0.1 deg",
        "plotting_options": {
            "min_error": 3,
        },
    }


@pytest.fixture
def sample_crab_config(sample_config):
    """Configuration for Crab Nebula observations."""
    config = sample_config.copy()
    config["target_source"]["name"] = "Crab Nebula"
    return config


@pytest.fixture
def sample_target_source():
    """Sample FixedTarget for testing."""
    coordinates = SkyCoord(ra="05h34m31.94s", dec="+22d00m52.2s", frame="icrs")
    return FixedTarget(coord=coordinates, name="Test Source")


@pytest.fixture
def sample_observer():
    """Sample Observer for testing (Roque de los Muchachos Observatory)."""
    return Observer.at_site("Roque de los Muchachos")


@pytest.fixture
def sample_time_range():
    """Sample time range for testing."""
    start_time = Time("2024-01-01 00:00:00")
    end_time = Time("2024-12-31 23:59:59")
    return [start_time, end_time]


@pytest.fixture
def sample_performance_data():
    """Sample performance data table for testing."""
    n_bins = 5
    min_energy = np.logspace(-1, 1, n_bins) * u.TeV
    max_energy = np.logspace(-0.8, 1.2, n_bins) * u.TeV
    gamma_rate = np.array([100, 50, 25, 10, 5]) * u.Unit("1/min")
    background_rate = np.array([1000, 500, 250, 100, 50]) * u.Unit("1/min")

    table = QTable(
        {
            "min_energy": min_energy,
            "max_energy": max_energy,
            "gamma_rate": gamma_rate,
            "background_rate": background_rate,
        }
    )
    return table


@pytest.fixture
def sample_horizon_profile():
    """Sample horizon profile data for testing."""
    azimuth = np.linspace(0, 360, 100) * u.deg
    zenith = (85 + 5 * np.sin(np.linspace(0, 2 * np.pi, 100))) * u.deg
    return QTable({"azimuth": azimuth, "zenith": zenith})


@pytest.fixture
def sample_ebl_file(tmp_path):
    """Create a sample EBL file for testing."""
    ebl_file = tmp_path / "test_ebl.txt"
    content = """0.0 0.01 0.03 0.05 0.1 0.2 0.5 1.0
10.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
100.0 0.01 0.02 0.03 0.05 0.08 0.15 0.30
1000.0 0.1 0.2 0.3 0.5 0.8 1.5 3.0
10000.0 1.0 2.0 3.0 5.0 8.0 15.0 30.0
"""
    ebl_file.write_text(content)
    return ebl_file


@pytest.fixture
def sample_yaml_config(tmp_path, sample_config):
    """Create a sample YAML configuration file for testing."""
    import yaml

    yaml_file = tmp_path / "test_config.yml"
    with open(yaml_file, "w") as f:
        yaml.safe_dump(sample_config, f)
    return yaml_file


@pytest.fixture
def sample_energy_bins():
    """Sample energy bins for testing."""
    return np.logspace(-1, 2, 11) * u.TeV


@pytest.fixture
def sample_spectrum():
    """Sample spectral model for testing."""
    from gammapy.modeling.models import PowerLawSpectralModel

    return PowerLawSpectralModel(
        amplitude="1e-11 TeV-1 cm-2 s-1",
        reference="1 TeV",
        index=2.5,
    )
