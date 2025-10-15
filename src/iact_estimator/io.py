"""Functions for input / output operations."""

import logging
from pathlib import Path

from numpy import loadtxt
from yaml import safe_load

__all__ = ["read_yaml", "load_ebl"]

logger = logging.getLogger(__name__)


def read_yaml(input_file_path):
    """
    Load data from a YAML file as a dictionary.

    Parameters
    ----------
    input_file_path : `str` or `pathlib.Path`
        Path to the input YAML file.

    Returns
    -------
    data : `dict`
        Contents of the YAML file in form
        of a Python dictionary.
    """

    input_file_path = Path(input_file_path).resolve()
    if not input_file_path.is_file():
        raise ValueError(f"Configuration file not found at {input_file_path}")

    with open(input_file_path, "r") as input_file:
        data = safe_load(input_file)

    return data


def load_ebl(ebl_file_path):
    """
    Load Extragalactic Background Light (EBL) data.

    Parameters
    ----------
    ebl_file_path : `str` or `~pathlib.Path`
        Path to an EBL data file.

    Returns
    -------
    zz : `np.array`
        TBD.
    energies : `~astopy.units.Quantity`
        Energy values.
    taus : `np.array`
        TBD.

    Notes
    -----
    This is a legacy function which works
    only with the default data file packages with
    *iact-estimator*.
    """
    if not Path(ebl_file_path).exists():
        raise ValueError("EBL file not found at", Path(ebl_file_path).absolute())

    with open(ebl_file_path, "r") as file:
        line = file.readline()
        firstlineEBL = list(map(float, line.split()))
        zz = firstlineEBL[1:]
    ebl_body = loadtxt(ebl_file_path, delimiter=" ", skiprows=1)
    energies = ebl_body[:, 0] * 1.0e3  # in GeV  #en*=1.e3; // GeV
    taus = ebl_body[:, 1:]
    if len(taus > 0):
        return zz, energies, taus
