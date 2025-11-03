"""Tests for the io module."""

import pytest
import numpy as np
from astropy.table import QTable
from astropy.io.fits import PrimaryHDU
import astropy.units as u


def test_read_yaml(tmp_path):
    """Test reading YAML configuration file."""
    from ..io import read_yaml

    test_file = tmp_path / "test_file.yml"
    test_file.write_text(
        """
        a: 1
        b:
            c: 3
            d: 4
        """
    )

    data = read_yaml(test_file)

    assert data["a"] == 1
    assert data["b"]["c"] == 3
    assert data["b"]["d"] == 4


def test_read_yaml_nonexistent_file():
    """Test reading a non-existent YAML file raises ValueError."""
    from ..io import read_yaml

    with pytest.raises(ValueError, match="Configuration file not found"):
        read_yaml("/nonexistent/path/to/file.yml")


def test_read_yaml_with_pathlib(tmp_path):
    """Test reading YAML file with pathlib.Path object."""
    from ..io import read_yaml

    test_file = tmp_path / "test_pathlib.yml"
    test_file.write_text("test_key: test_value\n")

    data = read_yaml(test_file)
    assert data["test_key"] == "test_value"


def test_load_ebl(sample_ebl_file):
    """Test loading EBL data file."""
    from ..io import load_ebl

    zz, energies, taus = load_ebl(sample_ebl_file)

    # Check redshift values
    assert len(zz) == 7
    assert zz[0] == 0.01
    assert zz[-1] == 1.0

    # Check energy values (multiplied by 1000 in load_ebl)
    assert len(energies) == 4
    assert energies[0] == 10000.0  # 10.0 * 1000 = 10000.0 GeV

    # Check tau values shape
    assert taus.shape == (4, 7)


def test_load_ebl_nonexistent_file():
    """Test loading non-existent EBL file raises ValueError."""
    from ..io import load_ebl

    with pytest.raises(ValueError):
        load_ebl("/nonexistent/ebl_file.txt")


def test_load_performance_ecsv(tmp_path):
    """Test loading performance data from ECSV file."""
    from ..io import load_performance_ecsv

    # Create a sample ECSV file
    ecsv_file = tmp_path / "test_performance.ecsv"
    table = QTable(
        {
            "min_energy": [0.1, 0.2, 0.5] * u.TeV,
            "max_energy": [0.2, 0.5, 1.0] * u.TeV,
            "gamma_rate": [100, 50, 25] * u.Unit("1/min"),
            "background_rate": [1000, 500, 250] * u.Unit("1/min"),
        }
    )
    table.write(ecsv_file, format="ascii.ecsv")

    # Load the file
    loaded_table = load_performance_ecsv(ecsv_file)

    # Check loaded data
    assert len(loaded_table) == 3
    assert "min_energy" in loaded_table.colnames
    assert "gamma_rate" in loaded_table.colnames
    assert loaded_table["min_energy"].unit == u.TeV
    assert np.allclose(loaded_table["gamma_rate"].value, [100, 50, 25])


def test_save_fits_hdu(tmp_path):
    """Test saving FITS HDU to file."""
    from ..io import save_fits_hdu

    # Create a simple FITS HDU
    data = np.arange(100).reshape(10, 10)
    hdu = PrimaryHDU(data)

    # Save to file
    output_file = tmp_path / "test_output.fits"
    save_fits_hdu(hdu, output_file)

    # Check file exists
    assert output_file.exists()

    # Verify content
    from astropy.io import fits

    with fits.open(output_file) as hdul:
        assert len(hdul) == 1
        assert np.array_equal(hdul[0].data, data)


def test_save_fits_hdu_overwrite(tmp_path):
    """Test saving FITS HDU with overwrite option."""
    from ..io import save_fits_hdu

    data = np.arange(25).reshape(5, 5)
    hdu = PrimaryHDU(data)
    output_file = tmp_path / "test_overwrite.fits"

    # Save first time
    save_fits_hdu(hdu, output_file)

    # Save again with overwrite=True
    new_data = np.arange(25, 50).reshape(5, 5)
    new_hdu = PrimaryHDU(new_data)
    save_fits_hdu(new_hdu, output_file, overwrite=True)

    # Verify new content
    from astropy.io import fits

    with fits.open(output_file) as hdul:
        assert np.array_equal(hdul[0].data, new_data)


def test_save_fits_hdu_invalid_directory():
    """Test saving FITS HDU to non-existent directory raises ValueError."""
    from ..io import save_fits_hdu

    hdu = PrimaryHDU(np.zeros((5, 5)))

    with pytest.raises(ValueError, match="Output directory"):
        save_fits_hdu(hdu, "/nonexistent/dir/file.fits")
