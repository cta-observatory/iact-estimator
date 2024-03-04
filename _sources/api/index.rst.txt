.. _reference:

=============
Reference API
=============

Command-line tools
==================

.. _iact-estimator-cfg:

iact-estimator-cfg
-------------------

.. argparse::
   :module: iact_estimator.scripts.get_config
   :func: parser
   :prog: iact-estimator-cfg

.. _iact-estimator:

iact-estimator
---------------

.. argparse::
   :module: iact_estimator.scripts.estimator
   :func: parser
   :prog: iact-estimator

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
