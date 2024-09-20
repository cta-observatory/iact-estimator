.. _installation:

Installation
============

Requirements
------------

- Python >= 3.11
- ``pip``, the Python package manager

It is always recommended to work from a virtual environment.

Users
-----

.. important::

    This type of installation is not yet available.
    A Python package on `PyPI <https://pypi.org/>`_ and on the
    `conda-forge <https://conda-forge.org/>`_ channel are planned.

.. _installation-dev:

Developers
----------

1. clone the project's repository and enter it
2. create a virtual environment
3. install the package in development mode

    .. code-block:: shell

        pip install -e '.[dev]'

.. tip::

    Any Python virtual environment can be used.

    If you like to use ``conda`` or ``mamba``,
    (the `miniforge <https://github.com/conda-forge/miniforge?tab=readme-ov-file#install>`_
    distribution is recommended)
    you can create a virtual environment
    with the following environment file, which you can find also in the root
    of the cloned repository.

    .. code-block:: shell

        mamba env create -f environment.yml

    .. literalinclude:: /../../environment.yml
       :language: yaml
