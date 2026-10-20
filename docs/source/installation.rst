.. _installation:

Installation
============

.. note::

   The project is under development and not yet released as a package on PyPI or any other package managers.
   The instructions refer to the development version from Github only.

.. code::

   git clone https://github.com/cta-observatory/iact-estimator.git
   cd iact-estimator


.. tab-set::

    .. tab-item:: conda/mamba

        .. code-block:: bash

            conda env create -f environment.yml
            conda activate iact-estimator
            pip install .

    .. tab-item:: pixi

        .. code-block:: bash

            pixi install
            pixi workspace register --name iact-estimator --path $PWD

            # then from wherever you want
            pixi shell -w iact-estimator
