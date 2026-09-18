# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# autodoc and autosummary import the package, so the repository root has to be
# on sys.path: ``htesp`` is a package now, not a flat ``src/`` directory.
import os
import sys

sys.path.insert(0, os.path.abspath(".."))


# -- Project information -----------------------------------------------------

project = 'HTESP'
copyright = '2024, Iowa State University'
author = 'Niraj K. Nepal, Lin-Lin Wang'

# The full version, including alpha/beta/rc tags
release = '2.0.0'
version = '2.0'

# Produced under U.S. Government contract DE-AC02-07CH11358 for Ames National
# Laboratory, operated by Iowa State University for the U.S. Department of
# Energy.  See LICENSE and docs/license.rst.


# -- General configuration ---------------------------------------------------
# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom ones.

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.extlinks',
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
]

autodoc_default_options = {
    'members': True,
    'member-order': 'bysource',
}

# autosummary lists modules that pull in pymatgen, ase, ifermi and friends.
# Missing optional dependencies should not break the build.
autodoc_mock_imports = [
    "pymatgen", "mp_api", "ase", "bsym", "lmfit", "spglib", "qmpy_rester",
    "ifermi", "plotly", "matminer", "sklearn", "phonopy",
]

autosummary_generate = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

root_doc = 'index'


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.
try:
    import sphinx_material  # noqa: F401
    html_theme = 'sphinx_material'
except ImportError:                      # pragma: no cover - build convenience
    html_theme = 'alabaster'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_show_sphinx = False
