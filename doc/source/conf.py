"""Configuration file for the Sphinx documentation builder.

For the full list of built-in configuration values, see the documentation:
https://www.sphinx-doc.org/en/master/usage/configuration.html

"""

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

from datetime import datetime
import os
import pathlib

from ansys_sphinx_theme import get_version_match
from sphinx.builders.latex import LaTeXBuilder
from sphinx_gallery import sorting as sg_sorting

LaTeXBuilder.supported_image_types = ["image/png", "image/pdf", "image/svg+xml"]

project = "PyChemkin"
copyright = f"(c) {datetime.now().year} ANSYS, Inc. All rights reserved"
author = "ANSYS, Inc. <ansys.support@ansys.com>"
cname = os.getenv("DOCUMENTATION_CNAME", default="chemkin.docs.pyansys.com")
switcher_version = get_version_match("0.1.0dev")

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration


html_context = {
    "github_user": "ansys",
    "github_repo": "pychemkin",
    "github_version": "main",
    "doc_path": "doc/source",
}
html_theme = "ansys_sphinx_theme"
html_short_title = html_title = "PyChemkin"
html_static_path = ["_static"]
templates_path = ["_templates"]
html_theme_options = {
    "logo": "pyansys",
    "switcher": {
        "json_url": f"https://{cname}/versions.json",
        "version_match": switcher_version,
    },
    "check_switcher": False,
    "github_url": "https://github.com/ansys/pychemkin",
    "show_prev_next": False,
    "show_breadcrumbs": True,
    "collapse_navigation": True,
    "use_edit_page_button": True,
    "additional_breadcrumbs": [
        ("PyAnsys", "https://docs.pyansys.com/"),
    ],
    "icon_links": [
        {
            "name": "Support",
            "url": "https://github.com/ansys/pychemkin/discussions",
            "icon": "fa fa-comment fa-fw",
        },
        {
            "name": "Download documentation in PDF",
            "url": f"https://{cname}/version/{switcher_version}/_static/assets/download/pychemkin.pdf",  # noqa: E501
            "icon": "fa fa-file-pdf fa-fw",
        },
    ],
    "ansys_sphinx_theme_autoapi": {
        "project": project,
    },
}

extensions = [
    "sphinx.ext.intersphinx",
    "sphinx_copybutton",
    "sphinx_design",
    "sphinx_jinja",
    "ansys_sphinx_theme.extension.autoapi",
]

# The suffix(es) of source filenames.
source_suffix = {
    ".rst": "restructuredtext",
    ".mystnb": "jupyter_notebook",
    ".md": "markdown",
}

# The master toctree document.
master_doc = "index"

# Configuration for Sphinx autoapi
suppress_warnings = [
    "autoapi.python_import_resolution",
]

# -- Sphinx gallery configuration --------------------------------------------
examples_source = pathlib.Path(__file__).parent.parent / "examples"
examples_output = pathlib.Path(__file__) / "examples"
example_subdir_names = ["modeling_features", "workflows", "use_cases"]
explicit_order = [
    "../../examples/chemistry",
    "../../examples/mixture",
    "../../examples/batch",
    "../../examples/engine",
    "../../examples/PSR",
    "../../examples/PFR",
    "../../examples/reactor_network",
    "../../examples/premixed_flame",
]
example_order = sg_sorting.ExplicitOrder(explicit_order)


# sphinx gallery options
sphinx_gallery_conf = {
    # convert rst to md for ipynb
    "pypandoc": True,
    # path to your examples scripts
    "examples_dirs": [str(examples_source / subdir) for subdir in example_subdir_names],
    # path where to save gallery generated examples
    "gallery_dirs": [str(examples_output / subdir) for subdir in example_subdir_names],
    "thumbnail_size": (320, 240),
    "remove_config_comments": True,
}

# -- Lincheck configuration --------------------------------------------------
linkcheck_ignore = [
    r"https://www.ansys.com/*",
    r"https://ansys.com/*",
    r"https://ansyshelp.ansys.com/*",
]
