# -*- coding: utf-8 -*-
#
# amrex documentation build configuration file, created by
# sphinx-quickstart on Thu Oct 19 14:30:08 2017.
#
# This file is execfile()d with the current directory set to its
# containing dir.
#
# Note that not all possible configuration values are present in this
# autogenerated file.
#
# All configuration values have a default; values that are commented out
# serve to show the default.

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))

import re
import sphinx_rtd_theme
import breathe
from datetime import datetime

def get_IAMR_version():
    today = datetime.today()
    return u'%s.%.2d-dev' % (str(today.year)[-2:], (today.month + 1) % 12)

# -- General configuration ------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.mathjax',
              'sphinxcontrib.bibtex',
              'sphinx.ext.githubpages',
              'sphinx.ext.viewcode',
              'sphinx.ext.intersphinx',
              'sphinx.ext.autosectionlabel',
              'breathe']

# bibtex
bibtex_bibfiles = ["refs.bib"]

# TODO: make IAMR tutorials
# intersphinx_mapping = {
#     'IAMR_tutorials': ('https://IAMR-codes.github.io/IAMR/tutorials_html/', None)
#     # 'IAMR_tutorials': ('../../../sphinx_tutorials/build/html/',
#     #                    '../../sphinx_tutorials/build/html/objects.inv')
# }
intersphinx_mapping = {
    'amrex': ('https://amrex-codes.github.io/amrex/docs_html/', None),
    'amrex_hydro':('https://amrex-codes.github.io/amrex/hydro_html/', None)
}

# Add any paths that contain templates here, relative to this directory.
templates_path = ['ytemplates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = '.rst'

# The main toctree document.
main_doc = 'index'

# General information about the project.
project = u'IAMR'
copyright = u'2017-2018, IAMR Team'
author = u'IAMR Team'

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
version = get_IAMR_version()
# The full version, including alpha/beta/rc tags.
release = get_IAMR_version()

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = None

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This patterns also effect to html_static_path and html_extra_path
exclude_patterns = []

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = False

numfig = True

# -- breathe options ------------------------------------------------------

breathe_projects = {
    "IAMR": "../../../out/docs_xml/doxygen/",
    }

breathe_default_project = "IAMR"

breathe_default_members = ('members', 'undoc-members', 'protected-members',
                           'private-members', 'content-only')

breathe_doxygen_config_options = {'EXTRACT_ALL': 'YES',
                                  'SHOW_USED_FILES': 'YES',
                                  'RECURSIVE': 'YES'}


# -- Options for HTML output ----------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
# html_theme_options = {}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_css_files = [
    'theme_overrides.css',  # overrides for wide tables in RTD theme
]

# Custom sidebar templates, must be a dictionary that maps document names
# to template names.
#
# This is required for the alabaster theme
# refs: http://alabaster.readthedocs.io/en/latest/installation.html#sidebars
html_sidebars = {
    '**': [
        'about.html',
        'navigation.html',
        'relations.html',  # needs 'show_related': True theme option to display
        'searchbox.html',
    ]
}


# -- Options for HTMLHelp output ------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'IAMRdoc'


# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',

    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (main_doc, 'IAMR.tex', u'IAMR Documentation',
     u'IAMR Team', 'manual'),
]

# -- Options for MathJax ---------------------------------------------------
mathjax3_config = {'tex': {'macros': {}}}

with open('mathsymbols.tex', 'r') as f:
    for line in f:
        macros = re.findall(r'\\newcommand{\\(.*?)}(\[(\d)\])?{(.+)}', line)
        for macro in macros:
            if len(macro[1]) == 0:
                mathjax3_config['tex']['macros'][macro[0]
                                                ] = "{" + macro[3] + "}"
            else:
                mathjax3_config['tex']['macros'][macro[0]] = [
                    "{" + macro[3] + "}", int(macro[2])]


numfig = True


# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (main_doc, 'IAMR', u'IAMR Documentation',
     [author], 1)
]


# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (main_doc, 'IAMR', u'IAMR Documentation',
     author, 'IAMR', 'One line description of project.',
     'Miscellaneous'),
]
