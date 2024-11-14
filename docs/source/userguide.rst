.. _userguide:

===========
Users guide
===========

Configuration
=============

You can use the :ref:`iact-estimator-cfg` command-line tool

.. code-block::

    iact-estimator config --to /where/to/save/config/file

to get the example configuration files:

- a global configuration file

.. literalinclude:: /../../src/iact_estimator/resources/config.yml
    :language: yaml

- a configuration files to define models to be used with the :ref:`gammapy_interface`,`

.. literalinclude:: /../../src/iact_estimator/resources/models.yml
    :language: yaml

Launch the simulation
=====================

You can use the :ref:`iact-estimator` command-line tool
to simulate the response of the telescope(s) to the
input source settings.

To produce the output in the current working directory
for a source named "my_source",

.. code-block::

    iact-estimator --config config.yml --source-name my_source

Gammapy interface
-----------------

In addition to the legacy scripts (which can still be used
in case all you have are gamma and background rates as a function of reconstructed energy)
_iact-estimator_ provides also an interface based on
`gammapy <https://gammapy.org/>`_.

In order to use it you only need to add the `--use-gammapy` flag
to the `run` subcommand.

This interface defines all sky models via the `models.yml` file,
ignoring the `assumed_model` section in the global configuration file.

.. note::

    The configuration schema will soon change,
    unifying the sky model specifications.

At the moment only the 1D analysis case without EBL absorption is supported.

Output
======

Terminal
--------

For each estimated energy range one gets the number of excess events,
signal-to-background ratio, significance, and information if a given
bin satisfies the conditions for the detection.

Plots
-----

The package comes with a plotting library
that allows to plot information about the observability of the source
and its spectral properties as seen by the telescopes.

Some of these functions are called by the main
script, but all of them can be imported separately
e.g. in a notebook (see the interactive example).

For a complete list of plotting functions, see
:py:mod:`iact_estimator.plots`.
