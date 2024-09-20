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
after installation,

- a configuration file ``config.yml``,
  which can be obtained with the :ref:`iact-estimator-cfg` command
- a `~astropy.table.QTable` containing the low zenith (0 to 30 degrees)
  performance (``from iact_estimator import LOW_ZENITH_PERFORMANCE``)
- a `~astropy.table.QTable` zenith (0 to 30 degrees)
  performance (``from iact_estimator import MID_ZENITH_PERFORMANCE``)

API description
===============

This section contains auto-generated API reference documentation [#f1]_.

.. toctree::
   :glob:

   **

.. [#f1] Created with `sphinx-autoapi <https://github.com/readthedocs/sphinx-autoapi>`_
