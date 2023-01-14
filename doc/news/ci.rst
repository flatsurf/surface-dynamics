**Removed:**

* Dropped testing for most optional dependencies with SageMath 9.1, 9.2, and 9.3. The compiled dependencies (libeantic, libintervalxt) were not on conda-forge at the time and are therefore not available for the correct version of libflint. (Before we were using our nightly builds from the flatsurf channel but these are fairly unreliable and not maintained.)

**Fixed:**

* Adapted CI setup on GitHub to changes in setup-miniconda.
