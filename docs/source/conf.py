# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
import os
import sys
sys.path.insert(0, os.path.abspath('../../cavsiopy/'))
# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'cavsiopy'
copyright = '2023, RRI Working Group'
author = 'E. Ceren Kalafatoglu Eyiguler, Warren Holley, Andrew D. Howarth, Donald W. Danskin, Kuldeep Pandey, Carley J. Martin, Glenn C. Hussey, Robert G. Gillies, Andrew W. Yau'
release = '1.2.2'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc']
autodoc_mock_imports = ["cartopy", "cavsiopy"]

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
