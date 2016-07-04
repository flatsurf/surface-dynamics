r"""
Databases for translation surfaces.

This file contain different databases (with different implementations) for
algorithms related to translation surfaces:

- structure of Strebel differentials for quadratic differentials in order to
  differentiate
- database of separatrix and cylinder diagrams up to isomorphism
- database of volume of connected components of Abelian strata
"""

from __future__ import absolute_import

import os
import sage.misc.misc

from sage.env import SAGE_SHARE

from . import __path__ as db_path
if len(db_path) != 1:
    raise RuntimeError("problem with setting paths")
db_path = db_path[0]

def line_count(filename):
    r"""
    Returns the number of lines in the file whose name is filename.
    """
    from sage.rings.integer import Integer
    f = open(filename)
    lines = 0
    read_f = f.read # loop optimization

    buf = read_f(0X100000)
    while buf:
        lines += buf.count('\n')
        buf = read_f(0X100000)  # 1024 x 1024

    f.close()

    return Integer(lines)

class GenericRepertoryDatabase:
    r"""
    Database that consists of a list of files in a repertory.
    """
    default_name = "generic_db"
    default_path = sage.misc.misc.SAGE_TMP

    def __init__(self, path=None, name=None, read_only=True):
        r"""
        INPUT:

        - ``path`` - string (default is SAGE_TMP) - path to the database. If the
          repertory does not exists, it is created.

        -  ``name`` - string (default is 'generic_db') - name of the repertory
           that will contain files of the database.

        - ``read_only`` - boolean
        """
        if path is None:
            path = self.default_path

        if name is None:
            name = self.default_name

        full_path = os.path.join(path,name)

        if not os.path.isdir(full_path):
            if read_only:
                raise ValueError("you must set read_only to `False` if the database does not already exist")
            try:
                os.makedirs(full_path)
            except OSError:
                raise ValueError("not able to create the database at {}".format(full_path))

        self.path = os.path.abspath(full_path)
        self.read_only = read_only

    def __eq__(self, other):
        r"""
        TESTS::

            sage: from surface_dynamics.databases.flat_surfaces import CylinderDiagrams
            sage: import os

            sage: rep1 = os.path.join(tmp_dir())
            sage: rep2 = os.path.join(tmp_dir())

            sage: C1 = CylinderDiagrams(rep1,read_only=False)
            sage: C2 = CylinderDiagrams(rep1,read_only=False)
            sage: C3 = CylinderDiagrams(rep2,read_only=False)
            sage: C1 == C1
            True
            sage: C1 == C2
            True
            sage: C1 == C3
            False
        """
        return self.__class__ is other.__class__ and other.path == self.path

    def __ne__(self, other):
        r"""
        TESTS::

            sage: from surface_dynamics.databases.flat_surfaces import CylinderDiagrams
            sage: import os

            sage: rep1 = os.path.join(tmp_dir())
            sage: rep2 = os.path.join(tmp_dir())


            sage: C1 = CylinderDiagrams(rep1,read_only=False)
            sage: C2 = CylinderDiagrams(rep1,read_only=False)
            sage: C3 = CylinderDiagrams(rep2,read_only=False)
            sage: C1 != C1
            False
            sage: C1 != C2
            False
            sage: C1 != C3
            True
        """
        return not self.__eq__(other)

    def clean(self):
        r"""
        Clean the database.

        EXAMPLES::

            sage: from surface_dynamics.all import *
            sage: from surface_dynamics.databases.flat_surfaces import CylinderDiagrams

            sage: C = CylinderDiagrams(tmp_dir(), "cylinder_diagrams", read_only=False)
            sage: C.update(AbelianStratum(4))
            sage: import os
            sage: os.listdir(C.path)
            ['cyl_diags-4-odd-1', 'cyl_diags-4-hyp-2', 'cyl_diags-4-odd-2',
             'cyl_diags-4-odd-3', 'cyl_diags-4-hyp-3', 'cyl_diags-4-hyp-1']
            sage: C.clean()
            sage: os.listdir(C.path)
            []
        """
        assert not self.read_only
        for filename in os.listdir(self.path):
            os.remove(os.path.join(self.path,filename))

class IrregularComponentTwins(GenericRepertoryDatabase):
    r"""
    Twin data of generalized permutation of irregular components of strata of
    Abelian differentials.
    """
    default_name = "generalized_permutation_twins"
    default_path = db_path

    def __repr__(self):
        r"""
        String representation

        TESTS::

            sage: from surface_dynamics.databases.flat_surfaces import IrregularComponentTwins
            sage: D = IrregularComponentTwins()
            sage: D.__repr__()
            'Database of twins of irregular components at ...'
        """
        return "Database of twins of irregular components at %s"%self.path

    def filename(self, stratum):
        r"""
        Returns the name of the file for the given component.

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: from surface_dynamics.databases.flat_surfaces import IrregularComponentTwins
            sage: D = IrregularComponentTwins()
            sage: D.filename(QuadraticStratum(12))
            'twins-12-irr'
            sage: D.filename(QuadraticStratum(2,2))
            Traceback (most recent call last):
            ...
            AssertionError: the stratum has no irregular component
        """
        from surface_dynamics.flat_surfaces.quadratic_strata import QuadraticStratum
        assert isinstance(stratum, QuadraticStratum)
        assert stratum.has_regular_and_irregular_components(), "the stratum has no irregular component"
        return ('twins-' +
                '_'.join(str(z) for z in stratum.zeros(poles=False)) +
                '_p' * stratum.nb_poles() +
                '-irr')

    def has_stratum(self, stratum):
        r"""
        Test whether the component is in the database.

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: from surface_dynamics.databases.flat_surfaces import IrregularComponentTwins
            sage: D = IrregularComponentTwins()
            sage: D.has_stratum(QuadraticStratum(12))
            True
        """
        return os.path.isfile(os.path.join(self.path,self.filename(stratum)))

    def update(self, stratum=None):
        r"""
        Update the database with the component comp.
        If comp is None update the whole database.

        The database should be not in read only mode.
        """
        assert not self.read_only

        f = open(os.path.join(selfr.path, self.filename(stratum)))

        p = stratum.irregular_component().permutation_representative()
        res = []
        for q in p.rauzy_diagram(symmetric=True):
            if q.is_cylindric():
                res.append(q._twin)
        res.sort()

        for twin in res:
            f.write(str(twin) + '\n')

    def list_strata(self):
        r"""
        Returns the list of components for which the list of twins is stored.

        EXAMPLES::

            sage: from surface_dynamics.databases.flat_surfaces import IrregularComponentTwins
            sage: G = IrregularComponentTwins()
            sage: G.list_strata()
            [Q_3(9, -1), Q_3(6, 3, -1), Q_4(12), Q_3(3^3, -1), Q_4(6^2), Q_4(9, 3), Q_4(6, 3^2), Q_4(3^4)]
        """
        from surface_dynamics.flat_surfaces.quadratic_strata import QuadraticStratum
        from sage.rings.integer import Integer
        s = set()
        for f in os.listdir(self.path):
            if f.startswith('twins-'):
                g = f[6:]
                i = g.index('-')
                comp = g[:i].replace('p','-1')
                s.add(comp)

        return sorted(QuadraticStratum(map(Integer,g.split('_'))) for g in s)

    def get(self, stratum):
        r"""
        Get the list of twins for the stratum of quadratic differentials ``stratum``.

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: from surface_dynamics.databases.flat_surfaces import IrregularComponentTwins
            sage: D = IrregularComponentTwins()
            sage: l = D.get(QuadraticStratum(6,3,-1))
            sage: l[0][0]
            [(1, 8), (1, 5), (1, 4), (0, 4), (0, 3), (1, 7), (1, 6)]
            sage: l[0][1]
            [(1, 2), (1, 3), (1, 0), (1, 1), (0, 2), (0, 1), (0, 6), (0, 5), (0, 0)]
            sage: len(l)
            1634
        """
        assert self.has_stratum(stratum)
        f = open(os.path.join(self.path, self.filename(stratum)))

        res = []
        s = f.readline()
        while s:
            res.append(eval(s))
            s = f.readline()
        return res

    def count(self, stratum):
        r"""
        Returns the number of twins for that stratum.

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: from surface_dynamics.databases.flat_surfaces import IrregularComponentTwins
            sage: D = IrregularComponentTwins()
            sage: Q = QuadraticStratum(12)
            sage: len(D.get(Q))
            6993
            sage: D.count(Q)
            6993
        """
        return line_count(os.path.join(self.path, self.filename(stratum)))

class CylinderDiagrams(GenericRepertoryDatabase):
    r"""
    Database of cylinder diagrams.

    The database consists of several files with the following name convention:
    "stratum-component-ncyls". As an example, the list of 3-cylinder diagrams in
    the odd component of H(2,2) is named "2_2-odd-3".

    EXAMPLES::

        sage: from surface_dynamics.all import *
        sage: from surface_dynamics.databases.flat_surfaces import CylinderDiagrams
        sage: import os
        sage: C = CylinderDiagrams()
        sage: a = AbelianStratum(3,1,1,1).unique_component()
        sage: C.filename(a, 2)
        'cyl_diags-3_1_1_1-c-2'
        sage: os.path.isfile(os.path.join(C.path, C.filename(a, 2)))
        True
        sage: l = list(C.get_iterator(a, 2))
        sage: l[0]
        (0,9)-(0,5,7,8,6) (1,6,2,5,3,8,4,7)-(1,2,3,9,4)
        sage: l[0].ncyls()
        2
        sage: l[0].stratum()
        H_4(3, 1^3)
    """
    default_name = "cylinder_diagrams"
    default_path = db_path

    def __repr__(self):
        r"""
        String representation.

        TEST::

            sage: from surface_dynamics.databases.flat_surfaces import CylinderDiagrams
            sage: C = CylinderDiagrams(tmp_dir(), read_only=False)
            sage: print C    # indirect doctest
            Database of cylinder diagrams at ...
        """
        return "Database of cylinder diagrams at %s"%self.path

    def filename(self, comp, ncyls):
        r"""
        Returns the name of the file for the given component ``comp`` and the
        given of number of cyliders ``ncyls``.

        EXAMPLES::

            sage: from surface_dynamics.all import *
            sage: from surface_dynamics.databases.flat_surfaces import CylinderDiagrams
            sage: C = CylinderDiagrams()
            sage: C.filename(AbelianStratum(4).odd_component(), 3)
            'cyl_diags-4-odd-3'
            sage: C.filename(AbelianStratum(3,3).hyperelliptic_component(), 2)
            'cyl_diags-3_3-hyp-2'
        """
        return ('cyl_diags-' +
                '_'.join(str(z) for z in comp.stratum().zeros()) +
                '-' + comp._name +
                '-' + str(ncyls))

    def has_component(self, comp):
        r"""
        Test wheter the database has the component ``comp``.

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: from surface_dynamics.databases.flat_surfaces import CylinderDiagrams
            sage: import os

            sage: C = CylinderDiagrams(tmp_dir(), read_only=False)
            sage: C.clean()

            sage: a1 = AbelianStratum(4).odd_component()
            sage: a2 = AbelianStratum(1,1,1,1).unique_component()

            sage: C.has_component(a1)
            False
            sage: C.has_component(a2)
            False
            sage: C.update(AbelianStratum(4))
            sage: C.has_component(a1)
            True

            sage: C.has_component(a2)
            False

            sage: C.has_component(-19)
            Traceback (most recent call last):
            ...
            AssertionError: the argument must be a component of stratum of Abelian differentials
        """
        from surface_dynamics.flat_surfaces.abelian_strata import AbelianStratumComponent

        assert isinstance(comp, AbelianStratumComponent), "the argument must be a component of stratum of Abelian differentials"
        return os.path.isfile(os.path.join(self.path,self.filename(comp, 1)))

    def has_stratum(self, stratum):
        r"""
        Test whether the database contains the data for a given stratum.

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: from surface_dynamics.databases.flat_surfaces import CylinderDiagrams
            sage: import os

            sage: C = CylinderDiagrams(tmp_dir(), read_only=False)
            sage: C.clean()

            sage: a1 = AbelianStratum(4)
            sage: a2 = AbelianStratum(1,1,1,1)

            sage: C.has_stratum(a1)
            False
            sage: C.has_stratum(a2)
            False
            sage: C.update(AbelianStratum(4))
            sage: C.has_stratum(a1)
            True

            sage: C.has_stratum(a2)
            False

            sage: C.has_stratum(1)
            Traceback (most recent call last):
            ...
            AssertionError: the argument must be a stratum of Abelian differential
        """
        from surface_dynamics.flat_surfaces.abelian_strata import AbelianStratum

        assert isinstance(stratum, AbelianStratum), "the argument must be a stratum of Abelian differential"
        return self.has_component(stratum.one_component())

    def _files(self):
        r"""
        Iterator over the list of files.

        TESTS::

            sage: from surface_dynamics.databases.flat_surfaces import CylinderDiagrams
            sage: C = CylinderDiagrams()
            sage: C._files()
            <generator object _files at ...>
            sage: it = C._files()
            sage: it.next()  # random
            'cyl_diags-...'
            sage: it.next()  # random
            'cyl_diags-...'
        """
        for f in os.listdir(self.path):
            if f.startswith('cyl_diags-'):
                yield f

    def list_strata(self):
        r"""
        List available strata in that database.

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: from surface_dynamics.databases.flat_surfaces import CylinderDiagrams
            sage: import os

            sage: C = CylinderDiagrams(tmp_dir(), read_only=False)
            sage: C.clean()

            sage: C.list_strata()
            []
            sage: C.update(AbelianStratum(1,1))
            sage: C.list_strata()
            [H_2(1^2)]
            sage: C.update(AbelianStratum(2))
            sage: C.list_strata()
            [H_2(2), H_2(1^2)]
        """
        from surface_dynamics.flat_surfaces.abelian_strata import AbelianStratum
        from sage.rings.integer import Integer
        s = set()
        for f in self._files():
            g = f[10:]
            s.add(g[:g.index('-')])

        return [AbelianStratum(map(Integer, g.split('_'))) for g in s]

    def get_iterator(self, comp, ncyls=None):
        r"""
        Returns an iterator over the list of cylinder diagrams for the component
        ``comp`` read from the database.

        INPUT:

        - ``comp`` - a component of stratum

        - ``ncyls`` - number of cylinders

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: from surface_dynamics.databases.flat_surfaces import CylinderDiagrams
            sage: import os

            sage: C = CylinderDiagrams(tmp_dir(), read_only=False)

            sage: A = AbelianStratum(2)
            sage: a = A.unique_component()
            sage: C.update(A)
            sage: list(C.get_iterator(a))
            [(0,2,1)-(0,2,1), (0,1)-(0,2) (2)-(1)]

            sage: C.clean()
            sage: C.get_iterator(a)
            Traceback (most recent call last):
            ...
            ValueError: not in the database
        """
        from surface_dynamics.flat_surfaces.strata import Stratum, StratumComponent

        if not isinstance(comp,StratumComponent):
            raise TypeError("comp should be a component of stratum")

        if ncyls is None:
            from itertools import chain
            g = comp.stratum().genus()
            s = comp.stratum().nb_zeros()
            return chain(*(self.get_iterator(comp,i) for i in xrange(1,g+s)))

        filename = os.path.join(self.path, self.filename(comp,ncyls))
        if not os.path.isfile(filename):
            raise ValueError("not in the database")
        return self._one_file_iterator(filename)

    def _one_file_iterator(self, filename):
        f = open(filename)
        from surface_dynamics.flat_surfaces.separatrix_diagram import CylinderDiagram
        line = f.readline()
        while line:
            yield CylinderDiagram(line[:-1])
            line = f.readline()
        f.close()

    def _check_symmetries(self, filename):
        r"""
        TESTS::

            sage: from surface_dynamics.databases.flat_surfaces import CylinderDiagrams
            sage: C = CylinderDiagrams()
            sage: C._check_symmetries('cyl_diags-4_2-odd-5')   # not tested -- to be fixed
            sage: for filename in C._files():          # not tested -- very long time, ~1h
            ....:     C._check_symmetries(filename)    # not tested -- very long time, ~1h
        """
        filename = os.path.join(self.path, filename)
        for c in self._one_file_iterator(filename):
            if c != c.canonical_label(inplace=False):
                raise ValueError("in file {}, {} does not have canonical labels".format(
                    filename, c))
            c1 = c.horizontal_symmetry()
            c1.canonical_label(inplace=True)
            if c1 < c:
                raise ValueError("in file {}, {} is not minimal... {} is the good one!".format(
                    filename, c,c1))
            c1 = c.vertical_symmetry()
            c1.canonical_label(inplace=True)
            if c1 < c:
                raise ValueError("in file {}, {} is not minimal... {} is the good one!".format(
                    filename, c,c1))
            c1 = c.inverse()
            c1.canonical_label(inplace=True)
            if c1 < c:
                raise ValueError("in file {}, {} is not minimal... {} is the good one!".format(
                    filename,c,c1))

    def _remove_symmetries(self, filename):
        r"""
        Use this at your own risk! It modifies the file of the database
        inplace!!!
        """
        filename = os.path.join(self.path, filename)
        n = 0
        s = set()
        garbage = set()
        for c in self._one_file_iterator(filename):
            n += 1
            if c in garbage:
                continue

            cc = [c]

            c1 = c.horizontal_symmetry()
            c1.canonical_label(inplace=True)
            cc.append(c1)

            c1 = c.vertical_symmetry()
            c1.canonical_label(inplace=True)
            cc.append(c1)

            c1 = c.inverse()
            c1.canonical_label(inplace=True)
            cc.append(c1)

            s.add(min(cc))
            garbage.update(cc)

        print "old cardinality:", n
        print "new cardinality:", len(s)
        print "gain           :", n-len(s)

        f = open(filename, "w")
        for c in sorted(s):
            f.write(str(c) + "\n")
        f.close()

    def update(self, stratum, verbose=False):
        r"""
        Compute once for all the given cylinder diagrams of the given
        ``stratum``.

        Warning::

            Depending on the dimension of the stratum, it may be very long!

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: from surface_dynamics.databases.flat_surfaces import CylinderDiagrams

            sage: C = CylinderDiagrams(tmp_dir(), read_only=False)
            sage: C.update(AbelianStratum(4), verbose=True)
            computation for H_3(4)
             ncyls = 1
             2 cyl. diags for H_3(4)^odd
             1 cyl. diags for H_3(4)^hyp
             ncyls = 2
             4 cyl. diags for H_3(4)^odd
             2 cyl. diags for H_3(4)^hyp
             ncyls = 3
             4 cyl. diags for H_3(4)^odd
             2 cyl. diags for H_3(4)^hyp
            sage: sorted(os.listdir(C.path))
            ['cyl_diags-4-hyp-1',
             'cyl_diags-4-hyp-2',
             'cyl_diags-4-hyp-3',
             'cyl_diags-4-odd-1',
             'cyl_diags-4-odd-2',
             'cyl_diags-4-odd-3']
        """
        import sys

        if verbose:
            print "computation for %s"%stratum
            sys.stdout.flush()

        for ncyls in xrange(1, stratum.genus() + stratum.nb_zeros()):
            if verbose:
                print " ncyls = %d"%ncyls
                sys.stdout.flush()
            d = stratum.cylinder_diagrams_by_component(ncyls, force_computation = True)
            for c in d:
                if verbose:
                    print " %d cyl. diags for %s"%(len(d[c]),c)
                f = open(os.path.join(self.path,self.filename(c,ncyls)), "w")
                for cyl in sorted(d[c]):
                    f.write(str(cyl) + "\n")
                f.close()

    def count(self, comp, ncyls=None):
        r"""
        Returns the number of cylinder diagrams for a stratum or a component of
        stratum with given number of cylinders.
        """
        from surface_dynamics.flat_surfaces.abelian_strata import AbelianStratum

        if isinstance(comp, AbelianStratum):
            return sum(self.count(cc, ncyls) for cc in comp.components())

        if ncyls is None:
            g = comp.stratum().genus()
            s = comp.stratum().nb_zeros()
            return sum((self.count(comp,i) for i in xrange(1,g+s)))

        return line_count(os.path.join(self.path, self.filename(comp, ncyls)))
