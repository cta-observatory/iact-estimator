.. iact-estimator documentation master file, created by
   sphinx-quickstart on Fri Jun 23 09:39:53 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Documentation of ``iact-estimator``
===================================

|ci| |docs|

.. |ci| image:: https://github.com/cta-observatory/iact-estimator/actions/workflows/ci.yml/badge.svg?branch=main
    :alt: CI
    :target: https://github.com/cta-observatory/iact-estimator/actions/workflows/ci.yml

.. |docs| image:: https://readthedocs.org/projects/iact-estimator/badge/?version=latest
    :target: https://iact-estimator.readthedocs.io/latest/?badge=latest
    :alt: Documentation Status

``iact-estimator`` is a Python3 package which allows
to evaluate the ability of an IACT telescope system
to detect a gamma-ray source.

It started from legacy scripts used by the
`MAGIC telescopes <https://magic.mpp.mpg.de/>`_
collaboration to aid users in the development
of their observation proposals.


About
-----

The output is an  estimate of what kind of signal can be observed by the
telescope system given a spectral shape.
The signal significances of each spectral point are computed according to
Eq. 17 definition from [2]_.

The currently available performance data publicly shipped with the package is
summarized by this table:

=================  ===============  ============  ==============  ==========
Instrument         Zenith range ID  Zenith range  Energy range    References
=================  ===============  ============  ==============  ==========
MAGIC              low              0 to 30 deg   40 GeV-16 TeV   [1]_
MAGIC              mid              30 to 45 deg  40 GeV-16 TeV   [1]_
MAGIC              high             ~60 deg       250 GeV-16 TeV  [4]_
MAGIC+LST1         low              0 to 30 deg   40 GeV-16 TeV   [5]_
MAGIC+LST1         mid              30 to 45 deg  80 GeV-16 TeV   [3]_
MAGIC+LST1         high             ~60 deg       250 GeV-16 TeV  [5]_
=================  ===============  ============  ==============  ==========

.. note::
   MAGIC+LST1 high zenith performance is based on Monte Carlo simulations with
   conservative background estimates. Please refer to :ref:`userguide` for
   detailed caveats when using this dataset.

Current caveats
---------------

.. warning::
   Please refer to the :ref:`userguide` for comprehensive information about
   caveats and limitations. Key points include:

- The tool operates on **estimated energy** by comparing with differential rates
  from Crab Nebula observations; for sources with soft spectra, energy migration
  effects may result in different performance

- The treatment of **extended sources** is approximate, accounting only for
  increased background (without energy-dependent PSF); for sources with
  extension > 0.4°, camera offset effects will further reduce performance

- Significances are given for **differential energy bins**; in practice, broader
  energy cuts and optimized analysis yield better integral sensitivity than
  simple integration of differential bins

- The **combined significance** is a crude approximation of detection capability

- **High zenith** observations have higher energy thresholds and specific
  limitations (see userguide)

References
----------

.. [1] Aleksić, J., et al., 2016, Astroparticle Physics, 72, 76
.. [2] Li, T.-P. & Ma, Y.-Q., 1983, ApJ, 272, 317
.. [3] Abe, H., et al., 2023, A&A, 670, A145
.. [4] MAGIC Crab observations at 55-62° zenith (2016-2018), Dataset: `Dr. J. van Scherpenberg PhD Thesis <https://nbn-resolving.org/urn:nbn:de:bvb:91-diss-20250606-1775491-0-4>`_
.. [5] MAGIC+LST1 Monte Carlo simulations (conservative estimates for high zenith)

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   :hidden:

   installation
   userguide
   example
   api/index
   contribute
   changelog
   authors

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
