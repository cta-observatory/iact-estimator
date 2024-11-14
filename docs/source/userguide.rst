.. _userguide:

===========
Users guide
===========

Configuration
=============

You can use the :ref:`iact-estimator-cfg` command-line tool

.. code-block::

    iact-estimator config --to /where/to/save/config/file

to get the following example configuration file,

.. literalinclude:: /../../src/iact_estimator/resources/config.yml
    :language: yaml

Launch the simulation
=====================

You can use the :ref:`iact-estimator` command-line tool
to simulate the response of the telescope(s) to the
input source settings.

To produce the output in the current working directory
for a source named "my_source",

.. code-block::

    iact-estimator run --config config.yml --source-name my_source

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
