import os
import sys

# Point Sphinx to your local modules if your notebooks import them
sys.path.insert(0, os.path.abspath("../../gpsp-atlas"))

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'gpsp-atlas'
copyright = '2026, Ekrem Oğuzhan Angüner'
author = 'Ekrem Oğuzhan Angüner'
release = '1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.mathjax",
    "myst_nb",
]

myst_enable_extensions = [
    "dollarmath",
    "amsmath",
]

templates_path = ['_templates']
exclude_patterns = []

# Tell myst-nb not to execute the notebooks during the build
# Since you want a static site with pre-computed outputs/images
nb_execution_mode = "off"

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ['_static']

html_theme_options = {
    "navigation_depth": 4
}

# Add this line to remove the default Sphinx/RTD footer text
html_show_sphinx = False

nb_execution_show_tb = False
nb_output_stderr = "remove"
