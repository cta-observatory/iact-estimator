"""Tests for the statistics module."""

import numpy as np
from ..statistics import (
    sigma_to_probability,
    probability_to_sigma,
    significance_li_ma,
)


class TestSigmaToProbability:
    """Tests for sigma_to_probability function."""

    def test_sigma_1(self):
        """Test conversion of 1-sigma significance."""
        probability = sigma_to_probability(1.0)
        assert np.isclose(probability, 0.1586553, atol=1e-6)

    def test_sigma_2(self):
        """Test conversion of 2-sigma significance."""
        probability = sigma_to_probability(2.0)
        assert np.isclose(probability, 0.0227501, atol=1e-6)

    def test_sigma_3(self):
        """Test conversion of 3-sigma significance."""
        probability = sigma_to_probability(3.0)
        assert np.isclose(probability, 0.0013499, atol=1e-6)

    def test_sigma_5(self):
        """Test conversion of 5-sigma significance."""
        probability = sigma_to_probability(5.0)
        assert np.isclose(probability, 2.866516e-07, atol=1e-10)

    def test_zero_sigma(self):
        """Test conversion of 0-sigma significance."""
        probability = sigma_to_probability(0.0)
        assert np.isclose(probability, 0.5)

    def test_negative_sigma(self):
        """Test conversion of negative sigma."""
        probability = sigma_to_probability(-1.0)
        assert probability > 0.5

    def test_array_input(self):
        """Test conversion with array input."""
        sigmas = np.array([1.0, 2.0, 3.0])
        probabilities = sigma_to_probability(sigmas)
        assert len(probabilities) == 3
        assert all(probabilities > 0)
        assert all(probabilities < 1)


class TestProbabilityToSigma:
    """Tests for probability_to_sigma function."""

    def test_inversion_1_sigma(self):
        """Test that probability_to_sigma inverts sigma_to_probability."""
        sigma_original = 1.0
        probability = sigma_to_probability(sigma_original)
        sigma_recovered = probability_to_sigma(probability)
        assert np.isclose(sigma_original, sigma_recovered)

    def test_inversion_5_sigma(self):
        """Test inversion for 5-sigma significance."""
        sigma_original = 5.0
        probability = sigma_to_probability(sigma_original)
        sigma_recovered = probability_to_sigma(probability)
        assert np.isclose(sigma_original, sigma_recovered, rtol=1e-6)

    def test_probability_half(self):
        """Test conversion of probability=0.5 (0-sigma)."""
        sigma = probability_to_sigma(0.5)
        assert np.isclose(sigma, 0.0)

    def test_small_probability(self):
        """Test conversion of small probability (high significance)."""
        probability = 1e-6
        sigma = probability_to_sigma(probability)
        assert sigma > 4.0

    def test_array_input(self):
        """Test conversion with array input."""
        probabilities = np.array([0.1586553, 0.0227501, 0.0013499])
        sigmas = probability_to_sigma(probabilities)
        expected_sigmas = np.array([1.0, 2.0, 3.0])
        assert np.allclose(sigmas, expected_sigmas, rtol=1e-4)


class TestSignificanceLiMa:
    """Tests for significance_li_ma function."""

    def test_no_excess(self):
        """Test with no excess events."""
        n_on = 100
        n_off = 100
        alpha = 1.0
        sigma = significance_li_ma(n_on, n_off, alpha)
        assert np.isclose(sigma, 0.0, atol=1e-2)

    def test_positive_excess(self):
        """Test with positive excess."""
        n_on = 150
        n_off = 100
        alpha = 1.0
        sigma = significance_li_ma(n_on, n_off, alpha)
        assert sigma > 0

    def test_off_region_scaling(self):
        """Test with different alpha (OFF region scaling)."""
        n_on = 200  # Signal + background
        n_off = 300  # More OFF events
        alpha = 0.5  # ON region half the size of OFF
        # Expected background in ON region: 0.5 * 300 = 150
        # So excess is 200 - 150 = 50 events
        sigma = significance_li_ma(n_on, n_off, alpha)
        assert sigma > 0

    def test_high_significance(self):
        """Test case with high significance."""
        n_on = 1000
        n_off = 100
        alpha = 1.0
        sigma = significance_li_ma(n_on, n_off, alpha)
        assert sigma > 20  # Should be high significance

    def test_typical_iact_case(self):
        """Test typical IACT observation case."""
        # Typical case: 100 excess, 300 background, 3 OFF regions
        n_on = 400  # signal + background
        n_off = 900  # 3 OFF regions with 300 each
        alpha = 1.0 / 3.0
        sigma = significance_li_ma(n_on, n_off, alpha)
        assert 4.0 < sigma < 6.0  # Expected range

    def test_with_expected_signal(self):
        """Test with expected signal parameter."""
        n_on = 150  # Observed ON events
        n_off = 100  # Observed OFF events
        alpha = 1.0
        mu_sig = 40  # Expected signal (less than observed excess of 50)
        sigma = significance_li_ma(n_on, n_off, alpha, mu_sig)
        # When observed excess > expected signal, should have positive significance
        assert sigma > 0

    def test_zero_counts(self):
        """Test with zero counts."""
        n_on = 0
        n_off = 0
        alpha = 1.0
        sigma = significance_li_ma(n_on, n_off, alpha)
        assert np.isclose(sigma, 0.0, atol=1e-10)

    def test_small_alpha(self):
        """Test with small alpha value (many OFF regions)."""
        n_on = 100
        n_off = 700  # 7 OFF regions
        alpha = 1.0 / 7.0
        sigma = significance_li_ma(n_on, n_off, alpha)
        # Small excess should give small significance
        assert sigma < 1.0


class TestStatisticsRoundTrip:
    """Tests for round-trip conversions."""

    def test_sigma_probability_roundtrip(self):
        """Test round-trip conversion for various sigma values."""
        test_sigmas = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0]
        for sigma_orig in test_sigmas:
            prob = sigma_to_probability(sigma_orig)
            sigma_recovered = probability_to_sigma(prob)
            assert np.isclose(sigma_orig, sigma_recovered, rtol=1e-6)

    def test_probability_sigma_roundtrip(self):
        """Test round-trip conversion for various probability values."""
        test_probs = [0.4, 0.3, 0.1, 0.01, 0.001]
        for prob_orig in test_probs:
            sigma = probability_to_sigma(prob_orig)
            prob_recovered = sigma_to_probability(sigma)
            assert np.isclose(prob_orig, prob_recovered, rtol=1e-6)
