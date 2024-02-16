def test_version():
    from iact_estimator import __version__

    assert __version__ != "unknown"
    assert __version__ != "0.0.0"
