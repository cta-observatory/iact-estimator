.. _reference:

=============
Reference API
=============

Command-line tools
==================

.. argparse::
   :module: iact_estimator.scripts.main
   :func: parser
   :prog: iact-estimator
   :nosubcommands:

.. _iact-estimator-cfg:

Get the default config file
---------------------------

.. argparse::
   :module: iact_estimator.scripts.main
   :func: parser
   :prog: iact-estimator
   :path: config

.. _iact-estimator:

Run the estimation process
--------------------------

.. argparse::
   :module: iact_estimator.scripts.main
   :func: parser
   :prog: iact-estimator
   :path: run

Package resources
=================

The package comes shipped with a set of package resources which are available
after installation.

Configuration File
------------------

A default configuration file ``config.yml`` can be obtained with the
:ref:`iact-estimator-cfg` command.

Performance Data Tables
-----------------------

All performance data are stored as `~astropy.table.QTable` objects and can be
imported directly:

MAGIC-only Performance
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    from iact_estimator import (
        LOW_ZENITH_PERFORMANCE,      # 0-30° zenith
        MID_ZENITH_PERFORMANCE,      # 30-45° zenith
        HIGH_ZENITH_PERFORMANCE,     # ~60° zenith
    )

MAGIC+LST1 Performance
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    from iact_estimator import (
        MAGIC_LST1_LOW_ZENITH_PERFORMANCE,   # 0-30° zenith
        MAGIC_LST1_MID_ZENITH_PERFORMANCE,   # 30-45° zenith
        MAGIC_LST1_HIGH_ZENITH_PERFORMANCE,  # ~60° zenith (MC-based)
    )

Each table contains the following columns:

- ``min_energy``: Lower energy bin edge (GeV)
- ``max_energy``: Upper energy bin edge (GeV)
- ``gamma_rate``: Rate of gamma-ray events (1/min)
- ``background_rate``: Rate of background events (1/min)

Data Sources and References
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table:: Performance Data Provenance
   :header-rows: 1
   :widths: 30 20 50

   * - Dataset
     - Zenith
     - Source
   * - MAGIC Low/Mid
     - 0-45°
     - Aleksić et al. 2016, Astroparticle Physics, 72, 76
   * - MAGIC High
     - 55-62°
     - Crab observations (2016-2018), 2.5 hours
   * - MAGIC+LST1 Low
     - 0-30°
     - Monte Carlo simulations, validated with data
   * - MAGIC+LST1 Mid
     - 30-45°
     - Abe et al. 2023, A&A, 670, A145
   * - MAGIC+LST1 High
     - ~59°
     - Monte Carlo simulations (conservative background)

.. note::
   When using MAGIC+LST1 high zenith performance data, please note that it is
   based on Monte Carlo simulations. Background rates above 1 TeV have been
   increased by 40% to provide conservative estimates.

API description
===============

This section contains auto-generated API reference documentation [#f1]_.

.. toctree::
   :glob:

   **

.. [#f1] Created with `sphinx-autoapi <https://github.com/readthedocs/sphinx-autoapi>`_
