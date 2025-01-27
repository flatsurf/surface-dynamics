<p align="center">
    <img alt="logo" src="https://user-images.githubusercontent.com/373765/255104540-451305f4-42e4-4c16-aee1-b38a0e6d41ad.svg">
</p>

<h1><p align="center">surface-dynamics</p></h1>

<p align="center">
  <img src="https://img.shields.io/badge/License-GPL_2.0_or_later-blue.svg" alt="License: GPL 2.0 or later">
  <a href="https://github.com/flatsurf/surface-dynamics/actions/workflows/test.yml"><img src="https://github.com/flatsurf/surface-dynamics/actions/workflows/test.yml/badge.svg" alt="Test"></a>
  <a href="https://doi.org/10.5281/zenodo.13356803"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.13356803.svg" alt="DOI 10.5281/zenodo.13356803"></a>
</p>

<p align="center">TODO</p>
<hr>

The ``surface-dynamics`` package for [SageMath](https://www.sagemath.org)
provides functionality related to interval exchange transformations,
translation surfaces, mapping classes and more. It is based on SageMath and
relies heavily on:

* [GMP](https://gmplib.org/) for arbitrary precision arithmetic
* [PARI/GP](https://pari.math.u-bordeaux.fr/) for number field computations
* [GAP](https://www.gap-system.org/) for finite groups representation and permutation groups
* [PPL](https://www.bugseng.com/ppl) (Parma Polyhedra Library) and
  [LattE](https://www.math.ucdavis.edu/~latte/) (Lattice point Enumeration) for
  polytope computations

## Usage

Here is an example session showcasing some of the computations that are
possible with surface-dynamics. For further examples, please consult [our
documentation](https://flatsurf.github.io/surface-dynamics/).

```python
sage: from surface_dynamics.all import *
sage: o = Origami('(1,2)', '(1,3)')
sage: o
(1,2)(3)
(1,3)(2)
sage: o.sum_of_lyapunov_exponents()
4/3
sage: o.lyapunov_exponents_approx()    # abs tol 0.05
[0.33441823619678734]
sage: o.veech_group()
Arithmetic subgroup with permutations of right cosets
 S2=(2,3)
 S3=(1,2,3)
 L=(1,2)
 R=(1,3)
sage: q = QuadraticStratum(1, 1, 1, 1)
sage: q.orientation_cover()
H_5(2^4)
sage: q.components()
[Q_2(1^4)^hyp]
sage: c = q.components()[0]
sage: c
Q_2(1^4)^hyp
sage: c.orientation_cover_component()
H_5(2^4)^odd

sage: AbelianStrata(genus=3).list()
[H_3(4), H_3(3, 1), H_3(2^2), H_3(2, 1^2), H_3(1^4)]

sage: O = OrigamiDatabase()
sage: q = O.query(("stratum", "=", AbelianStratum(2)), ("nb_squares", "=", 5))
sage: q.number_of()
2
sage: for o in q:
....:     print("%s\n- - - - - - - -" % o)
(1)(2)(3)(4,5)
(1,2,3,4)(5)
- - - - - - - -
(1)(2)(3,4,5)
(1,2,3)(4)(5)
- - - - - - - -
sage: Q12_reg = QuadraticStratum(12).regular_component()
sage: Q12_reg.lyapunov_exponents_H_plus(nb_iterations=2**20)   # abs tol 0.05
[0.6634, 0.4496, 0.2305, 0.0871]
sage: Q12_reg.lyapunov_exponents_H_minus(nb_iterations=2**20)  # abs tol 0.05
[1.0000, 0.3087, 0.1192]
```

## Installation

The easiest and recommended way to install surface-dynamics is to install
sage-flatsurf which includes surface-dynamics. Please follow the [instructions
for Linux and
macOS](https://flatsurf.github.io/sage-flatsurf/install.html#install-with-the-pixi-tarball)
or the [instructions for
Windows](https://flatsurf.github.io/sage-flatsurf/install.html#install-with-the-windows-installer).

If you have a working copy of SageMath already you can also try to install
surface-dynamics from PyPI using `pip`.

    $ sage -pip install surface-dynamics

## Build and Develop surface-dynamics with pixi

While surface-dynamics can be built and developed with `pip` like any other
Python package, we strongly recommend that you install [pixi](https://pixi.sh)
to get all the dependencies right.

Once you have cloned this source repository, you can use the following commands:

* `pixi run test` to build surface-dynamics and run its test suites
* `pixi run sage` to build surface-dynamics and spawn SageMath with the local surface-dynamics available
* `pixi run doc` to build and preview the documentation

<details>
<summary>What is pixi?</summary>

pixi is a tool for developers based on
[conda](https://en.wikipedia.org/wiki/Conda_(package_manager)) &
[conda-forge](https://conda-forge.org) so that we can all use the same
workflows in the same defined environments.

pixi allows us to ship a very opinionated setup to developers, namely a number
of opinionated scripts with corresponding tested (and opinionated)
dependencies.

This makes the whole development experience much more reliable and
reproducible, e.g., the CI on GitHub Pull Requests runs with the exact same
setup, so if something fails there, you can just run the CI command to
hopefully get exactly the same behavior locally.
</details>

<details>
<summary>How do I use pixi?</summary>

If you have not used pixi before, the most relevant pixi command is:

```sh
pixi run TASK
```

Run `pixi task list` to see the available tasks.

All tasks are defined in the `pyproject.toml` file and most are used somewhere in
our GitHub Continuous Integration setup, see .github/workflows/.
</details>

<details>
<summary>Why don't we add all these dependencies normally to pyproject.toml?</summary>

The dependency handling that Python provides when it comes to binary
dependencies is not very robust. At the moment, pixi/conda solve this problem
in a much better way.
</details>

<details>
<summary>Can I use pip and other tools with pixi?</summary>

More experienced developers may not want to use the provided tasks. You can
also just use the curated list of dependencies that pixi provides and drop into
a shell with these dependencies installed. For example, to run the doctests
directly, you could:

```sh
pixi shell -e dev
pip install -e .
sage -tp surface_dynamics
```
</details>

## Feedback and Contributions

If you have tried out surface-dynamics, we are thrilled to learn about your
experiences. If you ran into any problems or have suggestions for improvement,
please [create an issue](https://github.com/flatsurf/surface-dynamics/issues).

If you want to contribute to surface-dynamics, [pull
requests](https://github.com/flatsurf/surface-dynamics/pulls) are always
welcome :heart:

We are also happy to walk you through the process personally if you are unsure
how to get started. Feel free to reach out in the [#flatsurf stream on
Zulip](https://sagemath.zulipchat.com/#narrow/channel/271193-flatsurf) in any
case.

## Authors

See [AUTHORS](./AUTHORS) for a list of authors or visit [our zenodo
page](https://zenodo.org/badge/latestdoi/347440823).

## License

surface-dynamics is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License (GPL) as published by the
Free Software Foundation; either version 2.0 of the License, or (at your
option) any later version. See https://www.gnu.org/licenses/.

## How to cite this project

If you have used this project to prepare a publication please cite us as
described [on our zenodo page](https://zenodo.org/badge/latestdoi/347440823).

## Versions

The first release of ``surface-dynamics`` as a SageMath SPKG happened on the
30th of July 2015. Refer to our [Releases
Page](https://github.com/flatsurf/surface-dynamics/releases) for the latest
releases.

## Acknowledgements

* Julian Rüth's contributions to this project have been supported by the Simons
  Foundation Investigator grant of Alex Eskin.
