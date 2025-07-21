# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'simspace'
copyright = '2025, Tianxiao Zhao'
author = 'Tianxiao Zhao'
release = '0.2.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = []

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']

import os
import sys
sys.path.insert(0, os.path.abspath('../../'))

extensions = [
    'nbsphinx',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',  # For Google/NumPy docstring support
    'sphinx_autodoc_typehints',
    'sphinx.ext.viewcode',
]

html_theme = 'sphinx_rtd_theme'

extensions.append("myst_parser")
nbsphinx_execute = 'never' 
exclude_patterns = ['**.ipynb_checkpoints']