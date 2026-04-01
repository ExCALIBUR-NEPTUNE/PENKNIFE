# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'PENKNIFE'
copyright = '2026, UKAEA'
author = 'UKAEA'
release = ''

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["breathe"]

#templates_path = ['_templates']
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "pydata_sphinx_theme"
html_static_path = ['_static']
html_css_files = ["custom.css"]

html_sidebars = {
    "**": ["globaltoc.html"]
}

html_theme_options = {
    "navbar_start": ["navbar-logo"],
    "navbar_align": "left",
    "navigation_depth": 0,
}

import os
docs_version = "./docs_version"
if os.path.exists(docs_version):

    with open(docs_version) as fh:
        version = fh.read().strip()

    html_theme_options.update({
        "check_switcher": False,
        "switcher": {
            "json_url": "https://excalibur-neptune.github.io/PENKNIFE/switcher.json",
            "version_match": version,
        },
        "navbar_start": ["navbar-logo", "version-switcher"]
    })

breathe_projects = {"PENKNIFE": "../../build/doxygen/xml"}
breathe_default_project = "PENKNIFE"

# Enable referencing figures etc by number.
numfig = True
