def test_read_yaml(tmp_path):
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
