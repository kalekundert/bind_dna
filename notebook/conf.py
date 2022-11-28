#!/usr/bin/env python3.6

project = 'DNA-binding proteins'
copyright = '2019, Kale Kundert'
author = 'Kale Kundert'

templates_path = ['.templates']
source_suffix = '.rst'
master_doc = 'index'
exclude_patterns = [ #
        'build',
        'Thumbs.db',
        '.DS_Store',
        'README.*',
]

# The `\pu` MathJax command (provided by the mhchem extension) is only 
# supported in MathJax>=3.0, which is only supported in sphinx>=4.0.
needs_sphinx = '4.0'

extensions = [ #
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.todo',
    'sphinxcontrib.programoutput',
    'exmemo.sphinx.notebook',
    'exmemo.sphinx.biology',
    'exmemo.sphinx.general',
    'myst_parser',
]

suppress_warnings = ['ref.citation']
rst_epilog = """\
.. |br| raw:: html\n\n   <br />
.. |Cq| replace:: :math:`C_q`
"""
pygments_style = 'sphinx'
todo_include_todos = True
todo_link_only = True
myst_enable_extensions = [
        'strikethrough',
        'dollarmath',
        'deflist',
        'tasklist',
]

from sphinx_rtd_theme import get_html_theme_path
from exmemo.sphinx import favicon_path

html_theme = "sphinx_rtd_theme"
html_theme_path = [get_html_theme_path()]
html_favicon = str(favicon_path)
html_theme_options = {}
html_static_path = ['.static']

