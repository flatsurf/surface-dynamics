# -- Project information -----------------------------------------------------

project = 'surface_dynamics'
copyright = '2021, the surface_dynamics authors'
author = 'the surface_dynamics authors'

# The full version, including alpha/beta/rc tags
release = '0.4.7'

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.todo',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.extlinks',
    'sphinx_rtd_theme',
]

autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'show-inheritance': True
}

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = []

intersphinx_mapping = {
    'sage': 'https://doc.sagemath.org/html/en/reference/objects.inv',
}

# Shortcuts for external links
from sage.misc.sagedoc import extlinks

nitpick_ignore = [
        # Something is complaining when rendering a class that's inheriting in
        # certain ways, probably when SageObject is involved.
        # Let's ignore these warnings:
        # WARNING: py:class reference target not found: sage.structure.sage_object.SageObject
        ('py:class', 'sage.structure.sage_object.SageObject'),
        # WARNING: py:class reference target not found: sage.structure.parent.Parent
        ('py:class', 'sage.structure.parent.Parent'),
        # WARNING: py:class reference target not found: sage.structure.element.Element
        ('py:class', 'sage.structure.element.Element'),
        # WARNING: py:class reference target not found: sage.structure.unique_representation.UniqueRepresentation
        ('py:class', 'sage.structure.unique_representation.UniqueRepresentation'),
        # WARNING: py:class reference target not found: surface_dynamics.misc.sql_db.SQLQuery
        ('py:class', 'surface_dynamics.misc.sql_db.SQLQuery'),
        # WARNING: py:class reference target not found: surface_dynamics.misc.sql_db.SQLDatabase
        ('py:class', 'surface_dynamics.misc.sql_db.SQLDatabase'),
]

