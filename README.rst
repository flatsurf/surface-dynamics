================
Surface dynamics
================

The flatsurf package for SageMath adds functionality related to interval exchange
transformations, translation surfaces, mapping classes and more. It is
distributed as an external Python package. It installs on top of an existing
Sage installation.

This package is based on `SageMath <http://www.sagemath.org>`_ and relies heavily on:

* gmp or mpir for arbitrary precision arithmetic
* PARI/GP for number field computations
* GAP for finite groups representation and permutation groups

Prerequisites
-------------

Installing flatsurf requires a working Sage installation and a working gcc.
Under most Linux distribution, just install the package gcc
(it might even be already installed by default).
Under Mac OS X, get a working gcc by installing XCode from
the Mac App Store.

Installation
------------

TO BE UPDATED WITH MAKING FLATSURF A PURE PYTHON PACKAGE

Check
-----

After installing flatsurf, check that it works by launching Sage
and typing the following commands. You should get the same
output as below.::

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

Authors
-------

* Vincent Delecroix: maintainer
* Samuel Lelièvre: contribution for origamis and permutation representative 
  of quadratic strata
* Charles Fougeron: Lyapunov exponents for strata coverings

Versions
--------

* flatsurf 0.3 is in version beta
* flatsurf 0.2 was released on 2015-11-15.
* flatsurf 0.1 was released on 2015-07-30.