.. _userguide:

===========
Users guide
===========

Configuration
=============

You can use the :ref:`iact-estimator-cfg` command-line tool

.. code-block::

    iact-estimator-cfg --output-path /where/to/save/config/file

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

    iact-estimator --config config.yml --source-name my_source

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
