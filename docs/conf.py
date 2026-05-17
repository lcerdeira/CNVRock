"""Sphinx configuration for CNVRock documentation."""

from datetime import date

# ── Project ──────────────────────────────────────────────────────────────────
project = "CNVRock"
author = "Louise Cerdeira et al."
copyright = f"{date.today().year}, {author}"
release = "0.1"

# ── Extensions ───────────────────────────────────────────────────────────────
extensions = [
    "myst_parser",         # Markdown support
    "sphinx_copybutton",   # Copy-to-clipboard on code blocks
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
]

myst_enable_extensions = [
    "deflist",
    "colon_fence",
    "tasklist",
    "linkify",
    "substitution",
]

# Allow either .md or .rst sources
source_suffix = {
    ".md": "markdown",
    ".rst": "restructuredtext",
}

master_doc = "index"
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# ── HTML output ──────────────────────────────────────────────────────────────
html_theme = "sphinx_rtd_theme"
html_static_path = []
html_title = "CNVRock — KpSC AMR copy-number detection"

html_theme_options = {
    "navigation_depth": 4,
    "collapse_navigation": False,
    "sticky_navigation": True,
    "titles_only": False,
}

# ── intersphinx ──────────────────────────────────────────────────────────────
intersphinx_mapping = {
    "python":  ("https://docs.python.org/3", None),
    "numpy":   ("https://numpy.org/doc/stable", None),
    "pandas":  ("https://pandas.pydata.org/pandas-docs/stable", None),
    "pytorch": ("https://pytorch.org/docs/stable", None),
}

# Code copy button: drop prompt characters
copybutton_prompt_text = r">>> |\$ |# "
copybutton_prompt_is_regexp = True
