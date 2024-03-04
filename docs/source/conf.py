# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

from pathlib import Path
import tomllib

from iact_estimator.version import __version__

PROJECT_ROOT_DIR = Path(__file__).parent.parent.parent.resolve()

with open(PROJECT_ROOT_DIR / "pyproject.toml", "rb") as pyproject_toml:
    pyproject_cfg = tomllib.load(pyproject_toml)

project = pyproject_cfg["project"]["name"]
authors = "".join([author["name"] for author in pyproject_cfg["project"]["authors"]])
copyright = f"2023, {authors}"
author = authors
version = __version__
release = version

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_nb",
    "numpydoc",
    "sphinxarg.ext",
    "autoapi.extension",
    "sphinx_copybutton",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinxcontrib.towncrier",
]

templates_path = ["_templates"]
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "**.ipynb_checkpoints",
    "changes",
]
nitpicky = True

# Report warnings for all validation checks except "ES01", "SA01", "EX01"
# for more details, please see
# https://numpydoc.readthedocs.io/en/latest/validation.html#built-in-validation-checks
# numpydoc_validation_checks = {"all", "ES01", "SA01", "EX01", "RT02"}
numpydoc_show_class_members = False

# Options for myst-nb
# https://myst-nb.readthedocs.io/
nb_execution_mode = "force"
nb_execution_allow_errors = False
nb_execution_raise_on_error = True
nb_remove_code_source = False

# Options for towncrier extension
# https://sphinxcontrib-towncrier.readthedocs.io/en/latest/
towncrier_draft_config_path = "towncrier.toml"
towncrier_draft_autoversion_mode = "draft"  # or: 'sphinx-version', 'sphinx-release'
towncrier_draft_include_empty = True
towncrier_draft_working_directory = PROJECT_ROOT_DIR

# -- Options for sphinx-autoapi ---------------------------------------
# https://sphinx-autoapi.readthedocs.io/en/latest/index.html

autoapi_dirs = ["../../src"]
autoapi_root = "api"
autoapi_add_toctree_entry = False
autoapi_keep_files = True
autoapi_template_dir = templates_path[0]
autoapi_ignore = ["*scripts*", "*version*"]

# -- Options for intersphinx extension ---------------------------------------
# https://www.sphinx-doc.org/en/master/usage/extensions/intersphinx.html#configuration

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "astropy": ("https://docs.astropy.org/en/latest/", None),
    "gammapy": ("https://docs.python.org/3", None),
    "matplotlib": ("https://matplotlib.org/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "furo"
# html_static_path = ["_static"]
