---
jupytext:
  formats: ipynb,md:myst,sage:light
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.14.6
kernelspec:
  display_name: SageMath 9.7
  language: sage
  name: sagemath
---

# Square-tiled surfaces

Square-tiled surfaces or orgamis are translation surfaces obtained by gluing unit squares
by translation.

## Building square-tiled surfaces

Square-tiled surfaces are encoded by a pair of permutations and this is how you
build them in surface-dynamics
```{code-cell}
from surface_dynamics import Origami
o = Origami("(1,2)", "(1,3)")
```
The shortcut to get back the two permutations are r (for right) and u (for up)
```{code-cell}
print(o.r())
print(o.u())
```

Some well-know square-tiled surfaces are ready-made
```{code-cell}
ew = origamis.EierlegendeWollmilchsau()
print(ew)
print(ew.r())
print(ew.u())
```
and it is also possible to get representative or exhaustive list in a given stratum
or stratum component
```{code-cell}
H2_hyp = Stratum([2]).hyperelliptic_component()
print(H2_hyp.origamis(4))
```

## Topological invariants

Many topological invariants of square-tiled surfaces and their Teichmüller curves
are available. The stratum and stratum component
```{code-cell}
print(o.stratum())
print(o.stratum_component())
print(ew.stratum())
print(ew.stratum_component()
```

## Veech group and Teichmüller curve

Veech group
```{code-cell}
print(o.veech_group())
print(o.veech_group().is_congruence())
print(ew.veech_group())
```
Approximation of Lyapunov exponents and (exact rational) sum
```{code-cell}
print(o.lyapunov_exponents_approx())
print(o.sum_of_lyapunov_exponents())
print(ew.lyapunov_exponents_approx())
print(ew.sum_of_lyapunov_exponents())
```

If you are interested in some statistics of a Teichmüller curve you can iterate
through the origamis it contains. For example we study the distribution of the
number of cylinders in all Teichmüller curves of the component (genus 3) with
11 squares
```{code-cell}
cc = AbelianStratum(4).odd_component()
for T in cc.arithmetic_teichmueller_curves(11):
....:     cyls = [0]*3
....:     for o in T:
....:         n = len(o.cylinder_decomposition())
....:         cyls[n-1] += 1
....:     print cyls
```

## The origami database

The origami database is a database that contains the list of all arithmetic
Teichmüller curves (up to some number of squares). It is a standard sqlite
database and can also be read from other programs.
```{code-cell}
from surface_dynamics import OrigamiDatabase
D = OrigamiDatabase()
q = D.query(stratum=AbelianStratum(2), nb_squares=9)
print(q.number_of())
o1,o2 = q.list()
print(o1)
print(o2)
```

To get the list of columns available in the database you can do
```{code-cell}
D.cols()
```
Each of these columns is available for display
```{code-cell}
q = D.query(stratum=AbelianStratum(2))
D = OrigamiDatabase()
q = D.query(('stratum', '=', AbelianStratum(2)), ('nb_squares', '<', 15))
q.cols('nb_squares', 'veech_group_level', 'teich_curve_nu2', 'teich_curve_nu3', 'teich_curve_genus', 'monodromy_name')
q.show()
```

You can also get some information about the filling of the database with
```{code-cell}
D.info(genus=3)
```
