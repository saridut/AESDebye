# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
# sys.path.append(os.path.abspath('../../PyAESDebye'))
# sys.path.append(os.path.abspath('../../src'))
# sys.path.append(os.path.abspath('../../include'))

# sys.path.insert(0, os.path.abspath('../'))
sys.path.insert(0, os.path.abspath('../modules/'))

project = 'PyAESDebye'
copyright = '2024, Navid Panchi'
author = 'Navid Panchi'
release = '0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

templates_path = ['_templates']
# Add extensions
extensions = [
    # 'breathe',
    # 'exhale',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'nbsphinx'
]

# Paths to templates and static files
templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'Rose-X-ASF.txt']
add_module_names = False

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# HTML output
html_theme = 'furo'
html_static_path = ['_static']

autodoc_member_order = 'bysource'
nbsphinx_execute = 'never'


# # Breathe configuration
# breathe_projects = {"AESDebye": "./doxygen_docs/xml"}
# breathe_default_project = "AESDebye"

# # Exhale configuration
# exhale_args = {
#     "containmentFolder": "./api",
#     "rootFileName": "library_root.rst",
#     "rootFileTitle": "Library API",
#     "doxygenStripFromPath": "../../..",
#     "createTreeView": True,
# }

# # Setup the exhale extension
# primary_domain = 'cpp'
# highlight_language = 'cpp'

