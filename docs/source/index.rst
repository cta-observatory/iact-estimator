.. iact-estimator documentation master file, created by
   sphinx-quickstart on Fri Jun 23 09:39:53 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Documentation of ``iact-estimator``
===================================

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

The currently available performance data shipped with the package is:

- MAGIC *low* zenith angle (0 to 30 degrees) in the range 40GeV-16TeV from [1]_,
- MAGIC *mid* zenith angle (30 to 45 degrees) in the range 40GeV-16TeV from [1]_.

Current caveats
---------------

- the tool is operating on estimated energy doing simple
  comparisons with differential in estimated energy rates seen for the
  Crab Nebula; for softer sources the differences in energy migration
  will result in different performance than the ones produced by this tool

- the treatment of extended sources is very approximate, only the
  increase in background is taken into account (without energy
  dependence of Point Spread Function); for extended sources with extension >~ 0.4 deg the
  dependence on the offset from the centre of the camera will further
  worsen the performance w.r.t. one produced by this tool

- significances are given for each differential energy bin
  separately, but to detect a source one normally applies a cut that
  keeps a broad range of energies inside resulting in better integral
  sensitivity than differential one; also, optimization of cuts for a
  broad energy range usually results in somewhat
  better sensitivity than what one can get by simply integrating
  the used here signal in
  differential energy bins; as a very crude approximation for detection
  capability we calculate here also a sum of significances of all the
  points in Spectral Energy Distribution
  divided by the sqrt of the number of those points

References
----------

.. [1] Aleksic, J., et al., 2016, Astroparticle Physics, 72, 76
.. [2] Li & Ma 1983, ApJ, 272, 317

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
