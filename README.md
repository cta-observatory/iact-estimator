# IACT Telescopes estimator

[![CI](https://github.com/cta-observatory/iact-estimator/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/cta-observatory/iact-estimator/actions/workflows/ci.yml)

[![docs](https://readthedocs.org/projects/iact-estimator/badge/?version=latest)](https://iact-estimator.readthedocs.io/latest/?badge=latest)

``iact-estimator`` is a Python3-based package that allows
the estimation of the performance of an IACT telescope system
to detect a gamma-ray source.

It started from legacy scripts which were developed
by Dr. Julian Sitarek for the [MAGIC telescopes](https://magic.mpp.mpg.de/)
to aid users in preparing their observational proposals.

## Installation

For the moment there is no public packaged release.

1. clone the repository
2. (recommended) create a virtual environment (a `conda` environment file is provided)
3. `pip install -e .` (if you are a developer add `[dev]`)

## Usage

You can use the package using the available command line entry points.

Assuming you are working from a directory called `project`

1. get the default configuration file with `iact-estimator config`
2. launch the estimation with `iact-estimator run --config config.yml --source-name "Crab Nebula"`

> [!NOTE]
> Currently the source name must be a valid identifier from the SIMBAD, NED or VizieR databases.

For more details use the `-h` option to show the help menu.

The package is also a library, so you can use it to make your
scripts or in a notebook (see example notebook in the documentation).

The packaged data can be extended to other IACT systems,
provided the data has been published.

The support to load your own data files will be provided.

## Documentation

Documentation can be compiled locally after a
successful installation,

1. `make -C docs html`
2. open `docs/build/html/index.html` with your favourite browser

# How to contribute

This project is public and everyone can contribute to it.

Instructions on how to contribute to the project
are provided in the documentation.
