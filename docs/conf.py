import os
import sys
import tomllib

with open(os.path.abspath("../pyproject.toml"), "rb") as f:
    pyproject = tomllib.load(f)

project = pyproject["project"]["name"]
copyright = "2026, Navid Panchi"
author = pyproject["project"]["authors"][0]["name"]
release = pyproject["project"]["version"]
version = ".".join(release.split(".")[:2])

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.githubpages",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Options for HTML output -------------------------------------------------

html_theme = "sphinx_book_theme"
html_static_path = ["_static"]

# Theme options (HOOMD-blue style configuration)
html_theme_options = {
    "repository_url": "https://github.com/navidpanchi/AESDebye",
    "use_repository_button": True,
    "use_issues_button": True,
    "use_edit_page_button": True,
    "path_to_docs": "docs",
    "home_page_in_toc": True,
}

# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = True

# Autodoc settings to show __init__ signature and docstring
autoclass_content = "both"
autodoc_default_options = {
    "members": True,
    "special-members": "__init__",
    "undoc-members": True,
    "inherited-members": True,
}

