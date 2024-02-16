"""Estimate the performance of an IACT telescopes system to an observation."""

from importlib.resources import files
from pathlib import Path

from astropy.table import QTable

from .version import __version__

__all__ = [
    "__version__",
    "RESOURCES_PATH",
    "LOW_ZENITH_PERFORMANCE",
    "MID_ZENITH_PERFORMANCE",
]


RESOURCES_PATH = files("iact_estimator") / "resources"
"""Installation path of package resource files."""

HORIZON_PROFILE_M1 = QTable.read(RESOURCES_PATH / Path("horizon_profile_M1.ecsv"))
"""`~astropy.table.QTable` with energy zenith azimuth for
the horizon profile has seeen from M1."""

HORIZON_PROFILE_M2 = QTable.read(RESOURCES_PATH / Path("horizon_profile_M2.ecsv"))
"""`~astropy.table.QTable` with energy zenith azimuth for
the horizon profile has seeen from M1."""

LOW_ZENITH_PERFORMANCE = QTable.read(
    RESOURCES_PATH / Path("low_zenith_performance.ecsv")
)
"""`~astropy.table.QTable` with energy bins, rate of gammas and background for
low zenith performance (0 to 30 degrees)."""

MID_ZENITH_PERFORMANCE = QTable.read(
    RESOURCES_PATH / Path("mid_zenith_performance.ecsv")
)
"""`~astropy.table.QTable` with energy bins, rate of gammas and background for
mid zenith performance (30 to 45 degrees)."""
