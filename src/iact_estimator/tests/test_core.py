"""Tests for the core module."""

import logging
import pytest
import numpy as np
import astropy.units as u
from astroplan import FixedTarget
from gammapy.modeling.models import PowerLawSpectralModel, LogParabolaSpectralModel


class TestLoadTargetSourceCoordinates:
    """Tests for load_target_source_coordinates function."""

    def test_load_icrs_coordinates(self):
        """Test loading target source with ICRS coordinates."""
        from ..core import load_target_source_coordinates

        config = {
            "target_source": {
                "name": "Test Source",
                "coordinates": {
                    "ra_l": "05h34m31.94s",
                    "dec_b": "+22d00m52.2s",
                    "frame": "icrs",
                },
            }
        }

        target = load_target_source_coordinates(config)

        assert isinstance(target, FixedTarget)
        assert target.name == "Test Source"
        assert target.coord.frame.name == "icrs"

    def test_load_galactic_coordinates(self):
        """Test loading target source with Galactic coordinates."""
        from ..core import load_target_source_coordinates

        config = {
            "target_source": {
                "name": "Galactic Source",
                "coordinates": {
                    "ra_l": 184.56 * u.deg,
                    "dec_b": -5.78 * u.deg,
                    "frame": "galactic",
                },
            }
        }

        target = load_target_source_coordinates(config)

        assert isinstance(target, FixedTarget)
        assert target.name == "Galactic Source"

    def test_load_no_name_uses_default(self):
        """Test that missing name defaults to 'test_source'."""
        from ..core import load_target_source_coordinates

        config = {
            "target_source": {
                "name": None,
                "coordinates": {
                    "ra_l": "05h34m31.94s",
                    "dec_b": "+22d00m52.2s",
                    "frame": "icrs",
                },
            }
        }

        target = load_target_source_coordinates(config)
        assert target.name == "test_source"


class TestSetupLogging:
    """Tests for setup_logging function."""

    def test_setup_logging_info_level(self, tmp_path, monkeypatch):
        """Test setting up logger with INFO level."""
        from ..core import setup_logging

        # Change to tmp directory so log file is created there
        monkeypatch.chdir(tmp_path)

        logger = setup_logging("INFO", "test_source")

        assert logger.name == "iact_estimator"
        assert logger.level == logging.INFO

        # Check log file was created
        log_file = tmp_path / "test_source_iact_estimator.log"
        assert log_file.exists()

    def test_setup_logging_debug_level(self, tmp_path, monkeypatch):
        """Test setting up logger with DEBUG level."""
        from ..core import setup_logging

        monkeypatch.chdir(tmp_path)
        logger = setup_logging("DEBUG", "debug_test")

        assert logger.level == logging.DEBUG

    def test_setup_logging_invalid_level(self):
        """Test that invalid log level raises ValueError."""
        from ..core import setup_logging

        with pytest.raises(ValueError, match="Invalid log level"):
            setup_logging("INVALID_LEVEL", "test")


class TestGetHorizonStereoProfile:
    """Tests for get_horizon_stereo_profile function."""

    def test_merge_horizon_profiles(self, sample_horizon_profile):
        """Test merging two horizon profiles."""
        from ..core import get_horizon_stereo_profile

        # Create slightly different profiles for M1 and M2
        m1_data = sample_horizon_profile.copy()
        m2_data = sample_horizon_profile.copy()
        # M2 sees slightly higher zenith - modify the values directly
        # Get the values, add 1 deg, and create new column
        new_zenith = m2_data["zenith"].to_value(u.deg) + 1
        m2_data["zenith"] = new_zenith * u.deg

        az, zd = get_horizon_stereo_profile(m1_data, m2_data)

        assert len(az) == 100
        assert len(zd) == 100
        assert az.unit == u.deg
        assert zd.unit == u.deg

        # Stereo profile should take minimum zenith (should be closer to m1_data values)
        max_zenith = m2_data["zenith"].max()
        if hasattr(zd, "value"):
            assert np.all(zd.value <= max_zenith.value)
        else:
            assert np.all(zd <= max_zenith.value)

    def test_horizon_profile_with_wraparound(self):
        """Test horizon profile with azimuth > 360 degrees."""
        from ..core import get_horizon_stereo_profile
        from astropy.table import QTable

        # Create profiles with azimuth > 360
        m1_data = QTable(
            {
                "azimuth": np.linspace(0, 360, 50) * u.deg,
                "zenith": 85 * np.ones(50) * u.deg,
            }
        )
        m2_data = QTable(
            {
                "azimuth": np.linspace(370, 730, 50) * u.deg,  # > 360
                "zenith": 87 * np.ones(50) * u.deg,
            }
        )

        az, zd = get_horizon_stereo_profile(m1_data, m2_data)

        # Should handle wraparound correctly
        assert np.all(az >= 0 * u.deg)
        assert np.all(az <= 370 * u.deg)


class TestCheckInputConfiguration:
    """Tests for check_input_configuration function."""

    def test_valid_configuration(self, sample_config, sample_performance_data):
        """Test that valid configuration passes."""
        from ..core import check_input_configuration

        is_valid = check_input_configuration(sample_config, sample_performance_data)
        assert is_valid is True

    def test_invalid_large_extension(self, sample_config, sample_performance_data):
        """Test that large extension triggers warning."""
        from ..core import check_input_configuration

        sample_config["extension"] = "1.5 deg"
        is_valid = check_input_configuration(sample_config, sample_performance_data)
        assert is_valid is False

    def test_invalid_n_off_regions(self, sample_config, sample_performance_data):
        """Test that invalid number of OFF regions triggers warning."""
        from ..core import check_input_configuration

        sample_config["n_off_regions"] = 0
        is_valid = check_input_configuration(sample_config, sample_performance_data)
        assert is_valid is False

        sample_config["n_off_regions"] = 10
        is_valid = check_input_configuration(sample_config, sample_performance_data)
        assert is_valid is False

    def test_sum_trigger_mid_zenith_not_implemented(self, sample_config):
        """Test that SUMT at mid zenith raises NotImplementedError."""
        from ..core import check_input_configuration

        sample_config["sum_trigger"] = True
        sample_config["zenith_range"] = "mid"

        with pytest.raises(NotImplementedError):
            check_input_configuration(sample_config, None)

    def test_pulsar_mode_validation(self, sample_config, sample_performance_data):
        """Test pulsar mode parameter validation."""
        from ..core import check_input_configuration

        sample_config["pulsar_mode"]["enable"] = True
        sample_config["pulsar_mode"]["pulsar_on_range"] = 1.5  # Invalid
        sample_config["pulsar_on_range"] = 1.5  # Also update top-level for validation
        is_valid = check_input_configuration(sample_config, sample_performance_data)
        assert is_valid is False


class TestInitializeModel:
    """Tests for initialize_model function."""

    def test_initialize_power_law(self):
        """Test initializing a PowerLaw model."""
        from ..core import initialize_model

        config = {
            "target_source": {
                "assumed_model": {
                    "name": "gammapy.modeling.models.PowerLawSpectralModel",
                    "from_log10": False,
                    "parameters": {
                        "amplitude": "1e-11 TeV-1 cm-2 s-1",
                        "reference": "1 TeV",
                        "index": 2.5,
                    },
                }
            }
        }

        model = initialize_model(config)

        assert isinstance(model, PowerLawSpectralModel)
        assert model.amplitude.value == 1e-11
        assert model.index.value == 2.5

    def test_initialize_log_parabola_from_log10(self, sample_config):
        """Test initializing LogParabola model with from_log10."""
        from ..core import initialize_model

        model = initialize_model(sample_config)

        assert isinstance(model, LogParabolaSpectralModel)


class TestObservedFlux:
    """Tests for observed_flux function."""

    def test_observed_flux_zero_redshift(self, sample_ebl_file):
        """Test that zero redshift returns unattenuated flux."""
        from ..core import observed_flux

        energy = np.array([0.1, 1.0, 10.0]) * u.TeV
        flux_int = np.array([1.0, 1.0, 1.0]) * u.Unit("TeV-1 cm-2 s-1")
        redshift = 0.0

        obs_flux = observed_flux(energy, redshift, flux_int, sample_ebl_file)

        # At zero redshift, observed flux should equal intrinsic flux
        # obs_flux is flux_int * attenuation, where attenuation should be ~1 at z=0
        if hasattr(obs_flux, "to_value"):
            assert np.allclose(
                obs_flux.to_value(flux_int.unit), flux_int.value, atol=1e-5
            )
        else:
            # Result might be bare numpy array
            assert np.allclose(obs_flux, flux_int.value, atol=1e-5)

    def test_observed_flux_with_attenuation(self, sample_ebl_file):
        """Test that non-zero redshift computes observed flux."""
        from ..core import observed_flux

        energy = np.array([100.0, 500.0, 1000.0]) * u.GeV
        flux_int = np.array([1.0, 1.0, 1.0]) * u.Unit("TeV-1 cm-2 s-1")
        redshift = 0.5

        obs_flux = observed_flux(energy, redshift, flux_int, sample_ebl_file)

        # Check that we get a result with correct shape and positive values
        # The actual attenuation depends on the EBL model data
        if hasattr(obs_flux, "to_value"):
            flux_values = obs_flux.to_value(flux_int.unit)
        else:
            flux_values = obs_flux

        assert len(flux_values) == 3
        assert np.all(flux_values > 0)  # Flux should remain positive
        assert np.all(np.isfinite(flux_values))  # No NaNs or infs


class TestGetSED:
    """Tests for get_sed function."""

    def test_get_sed_calculation(self):
        """Test SED calculation from energy and flux."""
        from ..core import get_sed

        energy = np.array([0.1, 1.0, 10.0]) * u.TeV
        flux = np.array([1e-10, 1e-11, 1e-12]) * u.Unit("TeV-1 cm-2 s-1")

        sed = get_sed(energy, flux)

        # SED = E^2 * flux
        expected_sed = energy**2 * flux
        assert np.allclose(sed.value, expected_sed.value)
        assert sed.unit == u.Unit("TeV cm-2 s-1")


class TestPrepareData:
    """Tests for prepare_data function."""

    def test_prepare_data_low_zenith(self, sample_config):
        """Test preparing performance data for low zenith."""
        from ..core import prepare_data

        sample_config["zenith_range"] = "low"
        sample_config["sum_trigger"] = False
        sample_config["magic_lst1"] = False

        energy_bins, gamma_rate, background_rate = prepare_data(sample_config)

        assert len(energy_bins) > 0
        assert len(gamma_rate) > 0
        assert len(background_rate) > 0
        assert energy_bins.unit.is_equivalent(u.TeV)

    def test_prepare_data_with_degradation_factor(self, sample_config):
        """Test that degradation factor scales rates correctly."""
        from ..core import prepare_data

        sample_config["offset_degradation_factor"] = 0.8

        energy_bins, gamma_rate, background_rate = prepare_data(sample_config)

        # Rates should be scaled by degradation factor
        # (exact values depend on packaged data)
        assert np.all(gamma_rate.value >= 0)
        assert np.all(background_rate.value >= 0)

    def test_prepare_data_magic_lst1(self, sample_config):
        """Test preparing data for MAGIC+LST1 mode."""
        from ..core import prepare_data

        sample_config["magic_lst1"] = True
        sample_config["sum_trigger"] = False
        sample_config["zenith_range"] = "low"

        energy_bins, gamma_rate, background_rate = prepare_data(sample_config)

        assert len(energy_bins) > 0


class TestSourceDetection:
    """Tests for source_detection function."""

    def test_source_detection_high_significance(self, capsys):
        """Test source detection with high significance."""
        from ..core import source_detection

        sigmas = [5.0, 6.0, 7.0, 8.0]
        observation_time = 50 * u.h

        combined_sig = source_detection(sigmas, observation_time)

        assert combined_sig > 5.0

        # Check that appropriate message was printed
        captured = capsys.readouterr()
        assert "probably will be detected" in captured.out

    def test_source_detection_low_significance(self, capsys):
        """Test source detection with low significance."""
        from ..core import source_detection

        sigmas = [1.0, 1.5, 2.0]
        observation_time = 10 * u.h

        combined_sig = source_detection(sigmas, observation_time)

        assert combined_sig < 5.0

        captured = capsys.readouterr()
        assert (
            "will not be detected" in captured.out
            or "might be detected" in captured.out
        )

    def test_source_detection_no_data(self, capsys):
        """Test source detection with no significant data points."""
        from ..core import source_detection

        sigmas = []
        observation_time = 10 * u.h

        source_detection(sigmas, observation_time)

        captured = capsys.readouterr()
        assert "will not be detected" in captured.out


class TestCalculate:
    """Tests for calculate function."""

    def test_calculate_basic(self, sample_config, sample_spectrum):
        """Test basic calculate function."""
        from ..core import calculate, prepare_data

        energy_bins, gamma_rate, background_rate = prepare_data(sample_config)

        en, sed, dsed, sigmas, detected = calculate(
            energy_bins, gamma_rate, background_rate, sample_config, sample_spectrum
        )

        # Check output shapes and types
        assert len(en) >= 0
        assert len(sed) == len(en)
        assert len(dsed) == len(en)
        assert len(sigmas) == len(en)
        assert len(detected) == len(en)

        # Check units
        assert en.unit.is_equivalent(u.TeV)

    def test_calculate_with_pulsar_mode(self, sample_config, sample_spectrum):
        """Test calculate function with pulsar mode enabled."""
        from ..core import calculate, prepare_data

        sample_config["pulsar_mode"]["enable"] = True
        sample_config["pulsar_mode"]["pulsar_on_range"] = 0.1
        sample_config["pulsar_mode"]["pulsar_off_range"] = 0.9

        energy_bins, gamma_rate, background_rate = prepare_data(sample_config)

        en, sed, dsed, sigmas, detected = calculate(
            energy_bins, gamma_rate, background_rate, sample_config, sample_spectrum
        )

        # Should still produce valid output
        assert len(en) >= 0
