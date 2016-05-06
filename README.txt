= Surface dynamics =

The flatsurf package for Sage adds functionality related to
interval exchange transformations and translation surfaces.
It is distributed as an external package for Sage. It
installs on top of an existing Sage installation.

This is mainly developed by Vincent Delecroix (vincent.delecroix@labri.fr)
with contributions from Samuel Lelievre (samuel.lelievre@math.u-psud.fr).

Prerequisites
-------------

Installing flatsurf requires a working Sage installation
and a working gcc.
Under most Linux distribution, just install the package gcc
(it might even be already installed by default).
Under Mac OS X, get a working gcc by installing XCode from
the Mac App Store.

Installation
------------

Install flatsurf by typing the following in a terminal window

    $ sage -p http://www.labri.fr/perso/vdelecro/flatsurf-{VERSION}.spkg

(provided the command "sage" calls your version of Sage). You can also
download the last version of the source code with git

    $ git clone https://daemon@git.math.cnrs.fr/anon/plm/delecroix/flatsurf

or

    $ git clone https://github.com/videlec/flatsurf-package.git

Versions
--------

flatsurf 0.3 is in version beta

flatsurf 0.2 was released on 2015-11-15.

flatsurf 0.1 was released on 2015-07-30.
It should work on Sage version 6.0 or later.

Check
-----

After installing flatsurf, check that it works by launching Sage
and typing the following commands. You should get the same
output as below.

sage: from surface_dynamics.all import *
sage: o = Origami('(1,2)','(1,3)')
sage: print o
(1,2)(3)
(1,3)(2)
sage: o.sum_of_lyapunov_exponents()
4/3
sage: o.lyapunov_exponents_approx()
[0.33441823619678734]
sage: o.veech_group()
Arithmetic subgroup with permutations of right cosets
 S2=(2,3)
 S3=(1,2,3)
 L=(1,2)
 R=(1,3)
sage: QuadraticStratum(1,1,1,1).orientation_cover()
H_5(2^4)^odd

sage: AbelianStrata(genus=3).list()
[H_3(4), H_3(3, 1), H_3(2^2), H_3(2, 1^2), H_3(1^4)]

sage: O = OrigamiDatabase()
sage: q = O.query(("stratum","=",AbelianStratum(2)), ("nb_squares","=",5))
sage: q.number_of()
2
sage: for o in q: print o, "\n"
(1)(2)(3)(4,5)
(1,2,3,4)(5)

(1)(2)(3,4,5)
(1,2,3)(4)(5)

Contact
-------

Your comments and help are welcome: vincent.delecroix@labri.fr
For problems with Mac OS X: samuel.lelievre@gmail.com
