import sage_docbuild.conf

# -- Project information -----------------------------------------------------

project = 'surface-dynamics'
copyright = '2021-2023, the surface-dynamics authors'
author = 'the surface-dynamics authors'

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
    'myst_nb'
]

# Extensions when rendering .ipynb/.md notebooks
myst_enable_extensions = [
    "dollarmath",
    "amsmath",
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

# Allow linking to external projects, e.g., SageMath
intersphinx_mapping = {"sage": ("https://doc.sagemath.org/html/en/reference", None)}

# Show tracebacks on mystnb execution errors.
nb_execution_show_tb = True

# Raise an exception on any failed executing in notebook.
nb_execution_raise_on_error = True

# -- Options for HTML output ----------------------------------------------

# Imitate the look of the SageMath documentation.
html_theme = sage_docbuild.conf.html_theme
html_theme_options = sage_docbuild.conf.html_theme_options
pygments_style = sage_docbuild.conf.pygments_style
pygments_dark_style = sage_docbuild.conf.pygments_dark_style
html_css_files = sage_docbuild.conf.html_css_files

if html_css_files != ["custom-furo.css"]:
    raise NotImplementedError(
        "CSS customization has changed in SageMath. The configuration of surface-dynamics documentation build needs to be updated."
    )

html_css_files = ["https://doc.sagemath.org/html/en/reference/_static/custom-furo.css"]

# There is no surface-dynamics logo yet.
html_theme_options["light_logo"] = html_theme_options["dark_logo"] = "logo.svg"

# Output file base name for HTML help builder.
htmlhelp_basename = "sage-surfacedynamicsdoc"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["static"]

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

