# Configuration file for the Sphinx documentation builder.

import os
from pathlib import Path

BUILDDIR = Path(os.environ.get('ABS_BUILDDIR', '.')).absolute()

project = 'surface-dynamics'
copyright = '2021-2025, the surface-dynamics authors'
author = 'the surface-dynamics authors'

release = '0.7.0'

extensions = [
    # Sphinx's own autodoc does not pick up cached functions so we use the one from SageMath.
    'sage_docbuild.ext.sage_autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    'sphinx.ext.todo',
    'sphinx_book_theme',
    'myst_nb',
]

myst_enable_extensions = [
    "dollarmath",
    "amsmath",
]

autodoc_default_options = {
    'members': True,
    'undoc-members': True,
}

templates_path = ['_templates']

exclude_patterns = []

html_theme = 'sphinx_book_theme'

html_logo = 'https://github.com/flatsurf/surface-dynamics/raw/master/doc/source/static/logo.svg?sanitize=true'

html_theme_options = {
    "repository_url": "https://github.com/flatsurf/surface-dynamics",
    "icon_links": [{
        "name": "flatsurf",
        "url": "https://flatsurf.github.io",
        "icon": "https://flatsurf.github.io/assets/logo.svg",
        "type": "url",
    }, {
        "name": "GitHub",
        "url": "https://github.com/flatsurf/surface-dynamics",
        "icon": "fa-brands fa-square-github",
        "type": "fontawesome",
    }, {
        "name": "Zulip",
        "url": "https://sagemath.zulipchat.com/#narrow/channel/271193-flatsurf",
        "icon": "fa-regular fa-comments",
        "type": "fontawesome",
    },
    ],
    "use_edit_page_button": True,
    "repository_branch": "master",
    "path_to_docs": "doc/source",
}

intersphinx_mapping = {
    "sage": ("https://doc.sagemath.org/html/en/reference", None),
    "sage-flatsurf": ("https://flatsurf.github.io/sage-flatsurf/", None),
}

html_static_path = ["static"]
html_css_files = ['extra.css']

nb_execution_show_tb = True
nb_execution_raise_on_error = True
