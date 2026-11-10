"""Smoke tests for the package scaffold."""

import importlib


def test_package_imports_and_has_version():
    pkg = importlib.import_module("hbfgwas")
    assert isinstance(pkg.__version__, str)
    assert pkg.__version__.count(".") == 2


def test_minimal_runtime_dependencies_import():
    for name in ("numpy", "scipy", "pandas", "matplotlib", "yaml"):
        importlib.import_module(name)
