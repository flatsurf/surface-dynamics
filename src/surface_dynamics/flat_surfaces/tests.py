r"""
Testing abelian strata and related combinatorics
"""
import sage.misc.prandom as prandom
from sage.misc.misc import cputime

import random

from sage.combinat.permutation import Permutations

from copy import copy

from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.groups.perm_gps.permgroup import PermutationGroup
from sage.combinat.partition import Partitions

from separatrix_diagram import canonical_perm,SeparatrixDiagram,CylinderDiagram

class Test:
    r"""
    Testing class for ribbon graphs.
    """
    def __init__(self, min_num_seps=5, max_num_seps=10):
        r"""
        Create a ribbon graph testing object.

        INPUT:

        - ``num_darts`` - the number of darts of the randomly generated ribbon
          graphs.
        """
        self.min_num_seps=min_num_seps
        self.max_num_seps=max_num_seps

    def __repr__(self):
        r"""
        Return the string representation of self
        """
        return "Ribbon graph testing class"

    def _do(self, name):
        """
        Perform the test 'test_name', where name is specified as an
        argument. This function exists to avoid a call to eval.
        """

        print "test_%s"%name
        Test.__dict__["test_%s"%name](self)

    def _get_random_ribbon_graph(self):
        r"""
        Return a random ribbon graph with right parameters.
        """
        n = random.randint(self.min_num_seps,self.max_num_seps)
        S = SymmetricGroup(2*n)

        e = S([(2*i+1,2*i+2) for i in xrange(n)])
        f = S.random_element()
        P = PermutationGroup([e,f])

        while not P.is_transitive():
            f = S.random_element()
            P = PermutationGroup([e,f])

        return RibbonGraph(
                 edges=[e(i+1)-1 for i in xrange(2*n)],
                 faces=[f(i+1)-1 for i in xrange(2*n)])

    def _get_random_cylinder_diagram(self):
        r"""
        Return a random cylinder diagram with right parameters
        """
        test = False
        while test:
            n = random.randint(self.min_num_seps,self.max_num_seps)
            S = SymmetricGroup(2*n)

            bot = S.random_element()
            b = [[i-1 for i in c] for c in bot.cycle_tuples(singletons=True)]

            p = Partitions(2*n,length=len(b)).random_element()
            top = S([i+1 for i in canonical_perm(p)])
            t = [[i-1 for i in c] for c in top.cycle_tuples(singletons=True)]
            prandom.shuffle(t)

            c = CylinderDiagram(zip(b,t))
            test = c.is_connected()

        return c

    def random(self, seconds=0):
        """
        Perform random tests for a given number of seconds, or
        indefinitely if seconds is not specified.
        """
        self.test("random", seconds)

    def test(self, name, seconds=0):
        """
        Repeatedly run 'test_name', where name is passed as an
        argument. If seconds is nonzero, run for that many seconds. If
        seconds is 0, run indefinitely.
        """
        seconds = float(seconds)
        total = cputime()
        n = 1
        while seconds == 0 or cputime(total) < seconds:
            s = "** test_dimension: number %s"%n
            if seconds > 0:
                s += " (will stop after about %s seconds)"%seconds
            t = cputime()
            self._do(name)
            print "\ttime=%s\telapsed=%s"%(cputime(t),cputime(total))
            n += 1

    def test_random(self):
        """
        Do a random test from all the possible tests.

        EXAMPLES::

            sage: from sage.modular.arithgroup.tests import Test
            sage: Test().test_random() #random
            Doing random test
        """
        tests = [a for a in Test.__dict__.keys() if a[:5] == "test_" and a != "test_random"]
        name = prandom.choice(tests)
        print "Doing random test %s"%name
        Test.__dict__[name](self)

    def test_chain_complex(self):
        r"""
        Test dimensions of chain complex as well as euler characteristic and
        rank of the intersection form.
        """
        g = self._get_random_ribbon_graph()
        C = g.chain_complex()

        assert C.chain_space(0).rank() == g.num_vertices()
        assert C.chain_space(1).rank() == g.num_edges()
        assert C.chain_space(2).rank() == g.num_faces()

        assert C.cycle_space(0).rank() - C.boundary_space(0).rank() == 0  # this is reduced homology
        assert C.cycle_space(1).rank() - C.boundary_space(1).rank() == 2*g.genus()
        assert C.cycle_space(2).rank() - C.boundary_space(2).rank() == 1

        c,m = g.cycle_basis(intersection_matrix=True)
        assert m.is_skew_symmetric()
        assert m.rank() == 2*g.genus()

    def test_homological_dimension(self):
        from sage.rings.finite_rings.finite_field_constructor import GF

        g = self._get_random_cylinder_diagram()
        C = g.to_ribbon_graph().chain_complex(GF(2))
        Z1 = C.cycle_space(1)
        B1 = C.boundary_space(1)
        V = Z1.submodule(g.circumferences_of_cylinders(GF(2)))

        assert g.homological_dimension_of_cylinders() ==  V.rank() - V.intersection(B1).rank()
