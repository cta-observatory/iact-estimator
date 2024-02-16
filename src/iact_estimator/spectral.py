"""Spectral references and operations."""

import astropy.units as u

from gammapy.modeling.models import LogParabolaSpectralModel

__all__ = ["crab_nebula_spectrum"]


def crab_nebula_spectrum():
    """
    Spectrum of the Crab Nebula.

    Reference data comes from [1]_.

    Parameters
    ----------
    energy : `~astropy.units.Quantity`
      Energy values where to evaluate the spectrum.

    Returns
    -------
    spectrum : `~astropy.units.Quantity`
      Energy spectrum evaluated at ``energy`` values.

    References
    ----------

    .. [1] Aleksić, J. et al. The major upgrade of the MAGIC telescopes, Part II:
          A performance study using observations of the Crab Nebula.
          Astroparticle Physics 72, 76-94 (2015).
    """
    alpha = 2.51
    beta = 0.21
    amplitude = 3.39e-11 * u.Unit("TeV^-1 s^-1 cm^-2")
    reference = 1 * u.TeV

    spectrum = LogParabolaSpectralModel.from_log10(amplitude, reference, alpha, beta)

    return spectrum
