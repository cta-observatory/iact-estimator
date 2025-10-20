"""Functions that implement some statistical formulas."""

from gammapy.stats import WStatCountsStatistic
import numpy as np
from scipy.special import erfc, erfcinv

__all__ = ["sigma_to_probability", "probability_to_sigma", "significance_li_ma"]


def sigma_to_probability(sigma):
    """Convert significance to one side of the two-sided probability.

    Parameters
    ----------
    sigma : float
        Significance value.

    Returns
    -------
    probability : float
        Probability value."""
    probability = erfc(sigma / np.sqrt(2.0)) / 2
    return probability


#
def probability_to_sigma(probability):
    """Inversion function of `.sigma_to_probability()`.

    Parameters
    ----------
    probability : float
        Probability value.


    Returns
    -------
    sigma : float
        Significance value.
    """
    sigma = erfcinv(probability * 2) * np.sqrt(2)
    return sigma


def significance_li_ma(n_on, n_off, alpha, mu_sig=None):
    """
    Get the Li & Ma significance.

    This is equivalent to eq.17 of [1]_.

    Parameters
    ----------
    n_on : `int`
        Measured counts in ON region.
    n_off : `int`
        Measured counts in OFF region.
    alpha : `float`
        Acceptance ratio of ON and OFF measurements.
    mu_sig : `float`
        Expected signal counts in ON region.

    Returns
    -------
    sqrt_ts : `float``
        Significance as the square root of the Test Statistic.

    Notes
    -----
    The implementation uses `gammapy.stats.WStatCountsStatistic`
    and takes the square root of the Test Statistic.

    References
    ----------
    .. [1] Li, T.-P. & Ma, Y.-Q., ApJ, 1983, 272, 317, 10.1086/161295.
    """
    statistics = WStatCountsStatistic(n_on, n_off, alpha, mu_sig)
    sqrt_ts = statistics.sqrt_ts
    return sqrt_ts
