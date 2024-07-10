# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os
import sphinx_rtd_theme
import sys

from datetime import datetime

sys.path.insert(0, os.path.abspath('..'))

pyproject_path = '../qlayers'
readme_path = '../README.md'

project = '3DQLayers'
copyright = f'{datetime.now().year}, Alexander J Daniel'
author = 'Alexander J Daniel'
release = 'v0.0.1-rc.3'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.napoleon',
              'sphinx.ext.viewcode',
              'sphinx.ext.autodoc',
              'sphinx.ext.autosummary',
              'sphinx.ext.todo',
              'myst_parser',
              'nbsphinx',
              'nbsphinx_link']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

napoleon_include_init_with_doc = True

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_context = {
    "sidebar_external_links_caption": "Links",
    "sidebar_external_links"        : [
        (
            '<i class="fa fa-github fa-fw"></i> Source Code',
            "https://github.com/alexdaniel654/qlayers",
        ),
        (
            '<i class="fa fa-bug fa-fw"></i> Issue Tracker',
            "https://github.com/alexdaniel654/qlayers/issues",
        ),
        (
            '<i class="fa fa-file-text fa-fw"></i> Citation',
            "https://doi.org/10.5281/zenodo.12707172",
        ),
    ],
}
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
html_style = os.path.join("css", "custom.css")
html_favicon = os.path.join("_static", "favicon.ico")
html_logo = os.path.join("_static", "logo.png")
