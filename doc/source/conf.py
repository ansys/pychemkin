# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "PyChemkin"
copyright = "2024, Ansys Inc"
author = "ANSYS, Inc. <ansys.support@ansys.com>"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "autoapi.extension",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "nbsphinx",
]

templates_path = ["_templates"]
# exclude_patterns = []
nbsphinx_execute = "never"
autoapi_dirs = ["../../../ansys/chemkin"]
# autoapi_options = [ 'members', 'imported-members', 'inherited-members', 'undoc-memebers', #
#                   'special-members', 'private-members', 'show-inheritance', 'show-module-summary', ]
autoapi_ignore = ["*wrapper*", "*reactormodel*", "*color*"]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'alabaster'
# html_theme = "furo"
html_theme = "ansys_sphinx_theme"
html_static_path = ["_static"]
