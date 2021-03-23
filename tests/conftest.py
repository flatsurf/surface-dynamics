import contextlib
import sys
import pytest

# Copied from the MIT licensed setuptools tests:
# https://github.com/pypa/setuptools/blob/main/setuptools/tests/fixtures.py
@pytest.fixture(autouse=True, scope="session")
def workaround_xdist_376(request):
    """
    Workaround pytest-dev/pytest-xdist#376
    ``pytest-xdist`` tends to inject '' into ``sys.path``,
    which may break certain isolation expectations.
    Remove the entry so the import
    machinery behaves the same irrespective of xdist.
    """
    if not request.config.pluginmanager.has_plugin('xdist'):
        return

    with contextlib.suppress(ValueError):
        sys.path.remove('')
