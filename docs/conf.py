"""Sphinx configuration for NovInvenio docs.

Scaffolded for a future Read the Docs build. Not yet wired into CI —
see .readthedocs.yaml at the repo root for the RTD build config this
conf.py is meant to pair with.
"""
import os
import sys

sys.path.insert(0, os.path.abspath('..'))

project = 'NovInvenio'
copyright = '2026, Jason Stajich'
author = 'Jason Stajich'
release = '0.5.0'

extensions = [
    'myst_parser',
]

# MyST lets Sphinx build directly from the existing Markdown docs
# (METHOD_DESCRIPTION.md, docs/*.md, docs/adr/*.md, docs/agents/*.md)
# without converting them to .rst first.
source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}
myst_enable_extensions = [
    'colon_fence',
    'deflist',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

html_theme = 'sphinx_rtd_theme'
html_static_path = []
