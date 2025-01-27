import sage_docbuild.conf

# Use the same link shortening that SageMath uses
from sage.misc.sagedoc import extlinks


project = 'surface-dynamics'
copyright = "2021-2024, the surface-dynamics authors"
author = 'the surface-dynamics authors'
release = '0.6.0'


extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.todo',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.extlinks',
    'myst_nb'
]


myst_enable_extensions = [
    "dollarmath",
    "amsmath",
]


autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'show-inheritance': True
}


intersphinx_mapping = {"sage": ("https://doc.sagemath.org/html/en/reference", None)}


# Show tracebacks on mystnb execution errors.
nb_execution_show_tb = True


# Raise an exception on any failed executing in notebook.
nb_execution_raise_on_error = True


# Imitate the look of the SageMath documentation.
html_theme = sage_docbuild.conf.html_theme
html_theme_options = sage_docbuild.conf.html_theme_options
pygments_style = sage_docbuild.conf.pygments_style
pygments_dark_style = sage_docbuild.conf.pygments_dark_style
html_css_files = sage_docbuild.conf.html_css_files


if html_css_files != ["custom-furo.css", "custom-jupyter-sphinx.css", "custom-codemirror-monokai.css"]:
    raise NotImplementedError(
        "CSS customization has changed in SageMath. The configuration of surface-dynamics documentation build needs to be updated."
    )


html_css_files = [
    "https://doc.sagemath.org/html/en/reference/_static/custom-furo.css",
    "https://doc.sagemath.org/html/en/reference/_static/custom-jupyter-sphinx.css",
    "https://doc.sagemath.org/html/en/reference/_static/custom-codemirror-monokai.css",
]


html_theme_options["light_logo"] = html_theme_options["dark_logo"] = "logo.svg"


# Output file base name for HTML help builder.
htmlhelp_basename = "sage-surfacedynamicsdoc"


html_static_path = ["static"]


# Shortcuts for external links
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


def setup(app):
    app.connect('autodoc-process-docstring', sage_docbuild.conf.skip_TESTS_block)
