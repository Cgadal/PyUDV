# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

import sphinx_gallery
from sphinx_gallery.sorting import ExplicitOrder, FileNameSortKey

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

sys.path.insert(0, os.path.abspath("../"))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "pyudv"
copyright = "2023, Cyril Gadal"
author = "Cyril Gadal"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.napoleon",
    "numpydoc",
    "sphinx.ext.autodoc",  # Core Sphinx library for auto html doc generation from docstrings
    "sphinx.ext.autosummary",  # Create neat summary tables for modules/classes/methods etc
    # Link to other project's documentation (see mapping below)
    "sphinx.ext.intersphinx",
    # Add a link to the Python source code for classes, functions etc.
    "sphinx.ext.viewcode",
    # Automatically document param types (less noise in class signature)
    "sphinx_autodoc_typehints",
    # 'nbsphinx',  # Integrate Jupyter Notebooks and Sphinx
    "sphinx.ext.doctest",
    "sphinx.ext.mathjax",
    "sphinx_gallery.gen_gallery",
    "sphinx.ext.coverage",
]

autosummary_generate = True  # Turn on sphinx.ext.autosummary
autosummary_imported_members = True  # Also documents imports in __init__.py
# Remove 'view source code' from top of page (for html, not python)
html_show_sourcelink = False
autodoc_inherit_docstrings = True  # If no docstring, inherit from base class
# Enable 'expensive' imports for sphinx_autodoc_typehints
set_type_checking_flag = True
nbsphinx_allow_errors = True  # Continue through Jupyter errors
# autodoc_typehints = "description" # Sphinx-native method. Not as good as sphinx_autodoc_typehints
add_module_names = False  # Remove namespaces from class/method signatures

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"
highlight_language = "python3"

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Options for Sphinx gallery ---------------------------------------------

intersphinx_mapping = {
    "python":
    ("https://docs.python.org/{.major}".format(sys.version_info), None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "matplotlib": ("https://matplotlib.org/stable", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
}

examples_dirs = ["../examples"]
gallery_dirs = ["_examples"]

sphinx_gallery_conf = {
    "examples_dirs": examples_dirs,  # path to your example scripts
    "gallery_dirs":
    gallery_dirs,  # path to where to save gallery generated output
    # directory where function/class granular galleries are stored
    "backreferences_dir": "_gen_modules/backreferences",
    # Modules for which function/class level galleries are created.
    "doc_module": ("pyudv"),
    "reference_url": {
        "pyudv": None,  # The module you locally document uses None
        "numpy": "https://docs.scipy.org/doc/numpy/",
        "scipy": "https://docs.scipy.org/doc/scipy/",
        "matplotlib": "https://matplotlib.org/stable",
    },
    "matplotlib_animations": True,
    "plot_gallery": True,
    "ignore_pattern": "/_",
    "filename_pattern": "plot_",
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# Pydata theme
html_theme = "pydata_sphinx_theme"
html_theme_options = {
    "show_prev_next": False,
    "logo": {
        "text": "pyudv",
    },
    "github_url": "https://github.com/cgadal-pythonpackages/pyudv",
}

html_static_path = ["_static"]
html_css_files = ["pydata-custom.css"]

# If true, links to the reST sources are added to the pages.
html_show_sourcelink = False
html_copy_source = False
