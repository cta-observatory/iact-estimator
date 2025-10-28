.. _userguide:

===========
Users guide
===========

.. important::
   Most of the caveats and limitations described here refer to the legacy interface based
   on the original script that you can find `at this link <https://magic.mpp.mpg.de/fileadmin/user_upload/mss.py>`_.
   An modern interface based on gammapy and user-provided IRFs is underway.

Performance Data
================

Available Zenith Ranges
-----------------------

The package includes performance data for three zenith angle ranges:

- **Low zenith** (0-30°): Best sensitivity, lowest energy threshold (~50 GeV)
- **Mid zenith** (30-45°): Moderate sensitivity and energy threshold (~80 GeV)
- **High zenith** (~60°): Higher energy threshold (~250 GeV), suitable for multi-TeV observations

Each zenith range is available for both MAGIC-only and MAGIC+LST1 configurations.

.. note::
   The MAGIC+LST1 high zenith performance is based on Monte Carlo simulations
   computed at 59° zenith angle. The background rates above 1 TeV have been
   artificially increased by 40% to account for typical MC optimism observed
   when comparing MC predictions with real data at other zenith angles.

Caveats and Limitations
========================

Energy Migration Effects
------------------------

.. warning::
   The simulation operates on **estimated energy** by performing simple comparisons
   with differential rates seen for the Crab Nebula. For sources with soft spectra,
   differences in energy migration will result in different performance than
   calculated here.

Extended Sources
----------------

.. note::
   The treatment of extended sources is **approximate**. Only the increase in
   background is taken into account (without energy dependence of PSF).

   For extended sources with extension > 0.4°, the dependence on the offset
   from the center of the camera will further worsen the performance compared
   to point-like sources at the same offset.

   For large extensions, it is recommended to use only 1 OFF estimation region
   (``n_off_regions: 1`` in the configuration).

Differential vs Integral Sensitivity
-------------------------------------

.. note::
   Significances are given for each **differential energy bin separately**.
   However, to detect a source, one normally applies a cut that keeps a broad
   range of energies inside, resulting in better integral sensitivity than
   differential sensitivity.

   Additionally, optimization of analysis cuts for a broad energy range usually
   results in somewhat better sensitivity than what can be obtained by simply
   integrating the signal in differential energy bins.

   The combined significance provided by the package (using all detected points)
   is a crude approximation of detection capability.

Offset Degradation
------------------

.. warning::
   The default performance values assume observations at **0.4° offset** from
   the camera center (wobble mode).

   If observations are taken at higher offsets, performance will degrade. This
   can be approximated using the ``offset_degradation_factor`` parameter:

   - For best sensitivity region (~300 GeV) with MAGIC or MAGIC+LST1:
     ``degradation ≈ 1.1 × exp(-0.8 × offset²)``

   - For multi-TeV focused observations:
     ``degradation ≈ exp(-0.3 × offset²)``

   where offset is in degrees.

High Zenith Observations
------------------------

.. caution::
   High zenith angle observations have several specific considerations:

   - **Higher energy threshold**: ~250 GeV compared to ~50 GeV at low zenith
   - **MAGIC-only data**: Based on 2.5 hours of Crab observations at 55-62° zenith
   - **MAGIC+LST1 data**: MC-based with conservative background corrections
   - **No SUM trigger**: The SUM trigger mode is only available for low zenith observations
   - **Best for multi-TeV sources**: Particularly useful for sources above 1 TeV

Configuration
=============

You can use the :ref:`iact-estimator-cfg` command-line tool

.. code-block::

    iact-estimator config --to /where/to/save/config/file

to get the following example configuration file,

.. literalinclude:: /../../src/iact_estimator/resources/config.yml
    :language: yaml

Key Configuration Parameters
-----------------------------

Zenith Range Selection
^^^^^^^^^^^^^^^^^^^^^^

The ``zenith_range`` parameter determines which performance data to use:

.. code-block:: yaml

    zenith_range: "low"  # Options: "low", "mid", "high"

- **"low"**: 0-30° zenith, best sensitivity, energy threshold ~50 GeV
- **"mid"**: 30-45° zenith, moderate sensitivity, threshold ~80 GeV
- **"high"**: ~60° zenith, higher threshold ~250 GeV, best for TeV sources

MAGIC+LST1 Mode
^^^^^^^^^^^^^^^

Enable joint MAGIC+LST1 observations for improved sensitivity:

.. code-block:: yaml

    magic_lst1: True  # True for MAGIC+LST1, False for MAGIC-only

Observation Time
^^^^^^^^^^^^^^^^

Specify the observation time needed for your science case:

.. code-block:: yaml

    observation:
      time: "50 h"

The package will calculate expected significance and detection probability
for this observation time.

Extension and PSF
^^^^^^^^^^^^^^^^^

For extended sources, specify the source extension radius:

.. code-block:: yaml

    extension: 0.2 deg  # Source extension radius
    PSF: "0.1 deg"      # Point spread function

.. warning::
   Extended sources have reduced sensitivity due to increased background.
   Keep ``n_off_regions: 1`` for sources with extension > 0.4°.

Pulsar Mode
^^^^^^^^^^^

For pulsed sources (pulsars), enable pulsar mode which accounts for phase-folded observations:

.. code-block:: yaml

    pulsar_mode:
      enable: True
      pulsar_on_range: 0.092   # Phase range for ON region
      pulsar_off_range: 0.25   # Phase range for OFF region

.. note::
   In pulsar mode:

   - Background is reduced to account for the ON phase range
   - The signal-to-background ratio cut is ignored
   - Only the significance and minimum event cuts apply
   - Suitable for sources with known phase-folded emission

.. warning::
   Pulsar mode is only meaningful for pulsed sources at relatively small
   distances. If you're enabling pulsar mode with ``redshift > 0``, please
   verify this is intentional!

Extragalactic Sources and EBL Absorption
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For extragalactic sources, set the redshift to account for Extragalactic Background
Light (EBL) absorption:

.. code-block:: yaml

    redshift: 0.5  # Source redshift; use -1 for no EBL absorption

.. note::
   The package uses the Domínguez et al. (2011) EBL model to calculate opacity
   as a function of energy and redshift. The intrinsic source spectrum is
   attenuated by a factor of ``exp(-τ(E, z))`` where τ is the optical depth.

.. tip::
   For nearby sources (z < 0.01) or Galactic sources, set ``redshift: -1`` to
   disable EBL absorption calculations.

Launch the simulation
=====================

You can use the :ref:`iact-estimator` command-line tool
to simulate the response of the telescope(s) to the
input source settings.

To produce the output in the current working directory
for a source named "my_source",

.. code-block::

    iact-estimator run --config config.yml --source-name my_source

For detailed information about all possible options,
please consult the :ref:`reference` or issue the help
command with ``iact-estimator run -h``.

Using custom performance data
=============================

If you want to use your own performance data,
either for testing or because it is still not published,
you can load it from an ECSV with the format shown below.

Producing such a table is very easy with astropy
starting from the array quantities containing the data,
for example,

.. code-block:: python

    from astropy.table import QTable
    import astropy.units as u
    import numpy as np

    table = QTable([min_energy, max_energy, gamma_rate, bkg_rate],
                    names=("min_energy", "max_energy", "gamma_rate", "background_rate"),
                    meta={"name":"some_descriptive_title"})
    table.write("my_performance.ecsv")

this will result in a data file similar to this,

.. code-block::

    # %ECSV 1.0
    # ---
    # datatype:
    # - {name: min_energy, unit: GeV, datatype: float64}
    # - {name: max_energy, unit: GeV, datatype: float64}
    # - {name: gamma_rate, unit: 1 / min, datatype: float64}
    # - {name: background_rate, unit: 1 / min, datatype: float64}
    # meta:
    #   __serialized_columns__:
    #     background_rate:
    #       __class__: astropy.units.quantity.Quantity
    #       unit: !astropy.units.Unit {unit: 1 / min}
    #       value: !astropy.table.SerializedColumn {name: background_rate}
    #     gamma_rate:
    #       __class__: astropy.units.quantity.Quantity
    #       unit: !astropy.units.Unit {unit: 1 / min}
    #       value: !astropy.table.SerializedColumn {name: gamma_rate}
    #     max_energy:
    #       __class__: astropy.units.quantity.Quantity
    #       unit: &id001 !astropy.units.Unit {unit: GeV}
    #       value: !astropy.table.SerializedColumn {name: max_energy}
    #     min_energy:
    #       __class__: astropy.units.quantity.Quantity
    #       unit: *id001
    #       value: !astropy.table.SerializedColumn {name: min_energy}
    #   name: some_descriptive_title
    # schema: astropy-2.0
    min_energy max_energy gamma_rate background_rate
    39.8 63.1 0.818446 3.66424

You can then load it using the `--performance` flag of :ref:`iact-estimator`
to tell the command-line tool where to find the data file.

Output
======

Terminal
--------

For each estimated energy range one gets the number of excess events,
signal-to-background ratio, significance, and information if a given
bin satisfies the conditions for the detection.

Plots
-----

Tha package comes with a small plotting library
that allows to plot information about the observability of the source
and its spectral properties as seen by the telescopes.

Some of these functions are called by the main
script, but all of them can be imported separately
e.g. in a notebook (see the interactive example).

For a complete list of plotting functions, see
:py:mod:`iact_estimator.plots`.

Citing the Package
==================

If you use ``iact-estimator`` in a publication, please cite the relevant
performance references based on the configuration you used:

MAGIC Performance
-----------------

For MAGIC-only observations (low/mid zenith):

    Aleksić, J., et al., 2016, *Performance of the MAGIC stereo system obtained
    with Crab Nebula data*, Astroparticle Physics, 72, 76-94

For MAGIC-only high zenith observations:

    Dataset based on MAGIC Crab Nebula observations at 55-62° zenith angle
    (2016-2018), courtesy of J. van Scherpenberg

MAGIC+LST1 Performance
----------------------

For MAGIC+LST1 mid zenith observations:

    Abe, H., et al. (MAGIC+LST1 Collaboration), 2023, *First observations of the
    second solar spectrum with the MAGIC telescopes*, A&A, 670, A145

For MAGIC+LST1 low/high zenith observations:

    Monte Carlo simulations based on joint MAGIC+LST1 analysis framework

Statistical Methods
-------------------

For the Li & Ma significance calculation:

    Li, T.-P. & Ma, Y.-Q., 1983, *Analysis methods for results in gamma-ray
    astronomy*, ApJ, 272, 317-324

EBL Model
---------

If using EBL absorption (``redshift > 0``):

    Domínguez, A., et al., 2011, *Extragalactic background light inferred from
    AEGIS galaxy-SED-type fractions*, MNRAS, 410, 2556-2578

Troubleshooting
===============

Common Issues and Solutions
---------------------------

No Detection / Low Significance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If your source is not detected or has low significance:

1. **Check observation time**: Increase ``observation.time`` to improve statistics
2. **Verify spectral model**: Ensure your assumed spectrum is realistic
3. **Consider zenith range**: High zenith has higher energy threshold; try low zenith for sub-TeV sources
4. **Enable MAGIC+LST1**: May provide improved sensitivity
5. **Check energy range**: Ensure your source emits significantly in the accessible energy range

SUM Trigger Errors
^^^^^^^^^^^^^^^^^^

.. error::
   If you see an error about SUM trigger not being available:

   The SUM trigger mode is only implemented for MAGIC-only observations at low
   zenith. Set ``zenith_range: "low"`` and ``magic_lst1: False``, or disable
   SUM trigger with ``sum_trigger: False``.

Validation Errors
^^^^^^^^^^^^^^^^^

.. error::
   Configuration validation errors typically indicate incompatible settings:

- ``n_off_regions > 7``: Maximum is 7 OFF regions
- ``n_off_regions > 1`` with ``extension > 0.5 deg``: Use 1 OFF region for extended sources
- ``offset_degradation_factor > 1``: Factor must be ≤ 1 (performance can only degrade)
- Pulsar mode with high redshift: Verify this combination is intentional

Soft Source Limitations
^^^^^^^^^^^^^^^^^^^^^^^

.. warning::
   For sources with soft spectra:

   Differences in energy migration may result in different performance than
   calculated here. The package compares differential rates with Crab Nebula
   observations, which may not be accurate for sources with significantly
   different spectral shapes.

Performance Data Questions
^^^^^^^^^^^^^^^^^^^^^^^^^^

For questions about specific performance data or to report issues:

- MAGIC performance: Contact the MAGIC collaboration
- Package issues: Open an issue on the `GitHub repository <https://github.com/cta-observatory/iact-estimator>`_
