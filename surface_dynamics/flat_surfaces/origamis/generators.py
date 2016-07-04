#encoding=utf8
r"""
List of common origamis

Usage
=====

To see a list of all origami constructors, type "origamis." or "pillowcases."
and then press the tab key. The documentation for each constructor includes
information about each origami, which provides a useful reference.

Organization
============

The constructors available in this database are organized as follows

- :meth:`Escalators <OrigamiGenerators.Escalator>`
- :meth:`EierlegendeWollmilchsau <OrigamiGenerators.EierlegendeWollmilchsau>`
- :meth:`Podium <OrigamiGenerators.Podium>`
- :meth:`Stair <OrigamiGenerators.Stair>`
- :meth:`ProjectiveLine <OrigamiGenerators.ProjectiveLine>`
"""

from origami import Origami
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ

class OrigamiGenerators():
    r"""
    Examples of origamis.
    """
    def __repr__(self):
        r"""
        String representation
        """
        return "Origami generators"

    @staticmethod
    def Escalator(n):
        r"""
        Escalator origamis

        The escalator origami has 2n squares and is defined by

        r = (1 2)(3 4) ... (2n-1 2n)
        u = (2 3)(4 5) ... (2n 1)

        EXAMPLES::

            sage: from surface_dynamics.all import *
            sage: o = origamis.Escalator(3)
            sage: o
            Escalator with 3 steps
            sage: print str(o)
            (1,2)(3,4)(5,6)
            (1,6)(2,3)(4,5)
            sage: o.veech_group().index()
            3
        """
        lr = [(i,i+1) for i in xrange(1,2*n+1,2)]
        lu = [(1,2*n)] + [(i,i+1) for i in xrange(2,2*n,2)]
        positions = [((i+1)//2,i//2) for i in xrange(2*n)]
        name = "Escalator with %d steps" %n
        o = Origami(lr,lu,positions=positions,name=name)
        return o

    @staticmethod
    def EierlegendeWollmilchsau():
        r"""
        Eierlegende Wollmilchsau origami.

        This origami is associated to the quaternion group and was discovered by
        Gaby Schmithüsen.

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: o = origamis.EierlegendeWollmilchsau()
            sage: o
            Eierlegende Wollmilchsau
            sage: print str(o)
            (1,2,3,4)(5,6,7,8)
            (1,5,3,7)(2,8,4,6)
            sage: o.veech_group().index()
            1
            sage: o.sum_of_lyapunov_exponents()
            1

            sage: oo = origamis.CyclicCover([1,1,1,1])
            sage: o.is_isomorphic(oo)
            True
        """
        r = [1,2,3,0,5,6,7,4]
        u = [4,7,6,5,2,1,0,3]
        positions = [(0,0),(1,0),(2,0),(3,0),(0,1),(1,1),(2,1),(3,1)]
        name = "Eierlegende Wollmilchsau"
        o = Origami(r,u,check=False,as_tuple=True,positions=positions,name=name)
        return o

    @staticmethod
    def CyclicCover(a, M=None):
        r"""
        Return the Abelian differential associated to the quadruple ``a``
        
        INPUT:
            
            - ``a`` - quadruple of integers
            
            - ``M`` - positive integer (default: None)

        EXAMPLES:

            sage: from surface_dynamics.all import *

        There are two well known cyclic covers with completely degenerated
        Lyapunov spectrum::

            sage: o = origamis.CyclicCover([1,1,1,1])
            sage: o.sum_of_lyapunov_exponents()
            1

            sage: o = origamis.CyclicCover([1,1,1,3])
            sage: o.sum_of_lyapunov_exponents()
            1
        """
        from sage.arith.all import gcd

        if M is None:
            M = sum(a)

        a = map(Integer,a)
        if len(a) != 4:
            raise ValueError("a should be of length 4")
            
        if any(ai%2==0 for ai in a):
            raise ValueError("ai should be odd")
            
        if gcd([M]+a) != 1:
            raise ValueError("gcd(M,a) should be 0")
            
        if M%2:
            raise ValueError("M should be even")

        if sum(a) % M:
            raise ValueError("the sum of ai should be 0 mod M")
            
        r = [None] * (2*M)
        u = [None] * (2*M)
        
        # initialize horizontal cylinders
        i = 0
        seen = set([])
        pos = set([0])
        neg = set([])
        while pos or neg:
            # initialize pos 2i and 2i+1
            if pos:
                i = pos.pop()

                seen.add(i)
                r[2*i] = 2*i+1
                r[2*i+1] = (2*(i+a[0]+a[3])) % (2*M)  # sigma_A + sigma_D
                j = (i + a[0] + a[3]) % M
                if j not in seen and j not in pos:
                    pos.add(j)
        
                u[2*i] = (2*(i-a[3])+1) % (2*M) # sigma_D^-1
                j = (i-a[3]) % M
                if j not in seen and j not in neg:
                    neg.add(j)
                u[2*i+1] = (2*(i+a[3])) % (2*M) # sigma_D
                j = (i+a[3]) % M
                if j not in seen and j not in neg:
                    neg.add(j)
            
            if neg:
                i = neg.pop()
                seen.add(i)
                r[2*i+1] = 2*i
                r[2*i] = (2*(i+a[1]+a[2]) + 1) % (2*M) # sigma_B + sigma_C
                j = (i+a[1]+a[2]) % M
                if j not in seen and j not in neg:
                    neg.add(j)
                    
                u[2*i] = (2*(i+a[2]) + 1) % (2*M) # sigma_C
                j = (i+a[2]) % M
                if j not in seen and j not in neg:
                    pos.add(j)
                u[2*i+1] = (2*(i-a[2])) % (2*M) # sigma_C^-1
                j = (i-a[2]) %M
                if j not in seen and j not in neg:
                    pos.add(j)
        
        return Origami(r,u,as_tuple=True,name="M_%d(%d,%d,%d,%d)" %(M,a[0],a[1],a[2],a[3]))

    @staticmethod
    def Stair(n):
        r"""
        Stair origamis

        The stair origami with n=2k squares is defined by the permutations 

        r = (1 2)(3 4) ... (n-1 n)
        u = (1)(2 3) ... (n-2 n-1)(n)

        The stair origamis with n=2k+1 squares is defined by the permutations

        r = (1 2)(3 4) ... (n-2 n-1)(n)
        u = (1)(2 3) ... (n-1 n-2)

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: o = origamis.Stair(4)
            sage: o
            Stair origami with 4 squares
            sage: print str(o)
            (1,2)(3,4)
            (1)(2,3)(4)
            sage: o = origamis.Stair(5)
            sage: o
            Stair origami with 5 squares
            sage: print str(o)
            (1,2)(3,4)(5)
            (1)(2,3)(4,5)
        """
        r = []
        u = [0]

        for i in xrange(1,n-1,2):
            r.append(i)
            r.append(i-1)
            u.append(i+1)
            u.append(i)

        if n%2:
            r.append(n-1)
        else:
            r.append(n-1)
            r.append(n-2)
            u.append(n-1)

        positions = [((i+1)//2,i//2) for i in xrange(n)]
        o = Origami(r,u,as_tuple=True,positions=positions,name="Stair origami with %d squares" %n)
        return o

    @staticmethod
    def Podium(data):
        r"""
        If ``data`` is an integer than the standard podium with ``data`` steps is
        returned. Otherwise, ``data`` should be a weakly decreasing list of integers
        (i.e. a integer partition).

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: o = origamis.Podium([3,3,2,1])
            sage: o
            Podium origami with partition [3, 3, 2, 1]
            sage: print o
            (1,2,3)(4,5,6)(7,8)(9)
            (1,4,7,9)(2,5,8)(3,6)
        """
        from sage.combinat.partition import Partition

        if isinstance(data, (int,Integer)):
            p = Partition([i for i in xrange(data,0,-1)])
        else:
            p = Partition(data)

        p = Partition(data)
        q = p.conjugate()

        r=[]
        positions = []
        i = 0
        for j,jj in enumerate(p):
            r.extend(xrange(i+1,i+jj))
            r.append(i)
            i += jj
            positions.extend((k,j) for k in xrange(jj))

        u = [None]*sum(p)
        for j in xrange(len(q)):
            k = j
            for jj in xrange(q[j]-1):
                u[k] = k+p[jj]
                k += p[jj]
            u[k] = j

        return Origami(r,u,positions=positions,name="Podium origami with partition %s" %str(p),as_tuple=True)

    @staticmethod
    def ProjectiveLine(p, r=None, u=None):
        r"""
        Return the projective line with action given by the matrices ``r``
        and ``u`` on the projective line over the field with ``p`` elements.

        If ``r`` and ``u`` are None, then r and u are choosen to be

        .. MATH::

            r: z \mapsto z+1   \qquad    u: z \mapsto 1/z

        The monodromy of this origami is included in `PGL(2,\ZZ/p\ZZ)`.

        EXAMPLES::

            sage: from surface_dynamics.all import *

        By default, the matrix used to generate the projective line origami are
        given by `z -> z+1` and `z -> 1/z` which gives a family of origamis
        which belongs to the strata ``H(2^k)``::

            sage: o = origamis.ProjectiveLine(3); o
            Projective line origami on GF(3)
             r = [1 1]    u = [0 1]
                 [0 1]        [1 0]
            sage: o.stratum()
            H_2(2)
            sage: o.veech_group().index()
            9
            sage: o.sum_of_lyapunov_exponents()
            4/3

            sage: o = origamis.ProjectiveLine(5); o
            Projective line origami on GF(5)
             r = [1 1]    u = [0 1]
                 [0 1]        [1 0]
            sage: o.stratum()
            H_3(2^2)
            sage: o.veech_group().index()
            9
            sage: o.sum_of_lyapunov_exponents()
            5/3

            sage: o = origamis.ProjectiveLine(7); o
            Projective line origami on GF(7)
             r = [1 1]    u = [0 1]
                 [0 1]        [1 0]
            sage: o.stratum()
            H_3(2^2)
            sage: o.veech_group().index()
            45
            sage: o.sum_of_lyapunov_exponents()
            5/3

        Any invertible matrix mod p can be choosed to generate the origami::

            sage: r = matrix([[1,3],[0,1]])
            sage: u = matrix([[1,0],[3,1]])
            sage: o = origamis.ProjectiveLine(5,r,u); o
            Projective line origami on GF(5)
             r = [1 3]    u = [1 0]
                 [0 1]        [3 1]
            sage: o.veech_group().index()
            10
            sage: o.sum_of_lyapunov_exponents()
            9/5
            sage: o.stratum_component()
            H_3(4)^hyp
        """
        from sage.arith.all import is_prime
        from sage.rings.finite_rings.finite_field_constructor import GF
        from sage.groups.matrix_gps.linear import GL
        from sage.modules.free_module import VectorSpace

        p = ZZ(p)
        if not p.is_prime():
            raise ValueError("p (={}) must be a prime number".format(p))

        G = GL(2,GF(p))
        V = VectorSpace(GF(p),2)

        if r is None:
            mr = G([[1,1],[0,1]])
        else:
            mr = G(r)
        if u is None:
            mu = G([[0,1],[1,0]])
        else:
            mu = G(u)

        sr = str(mr).split('\n')
        su = str(mu).split('\n')

        mr = mr
        mu = mu

        if r is None and u is None:
            positions = [(i,0) for i in xrange(p+1)]
        else:
            positions = None

        r = []
        u = []
        for i in xrange(p):
            v = V((i,1)) * mr
            if v[1]:
                r.append(int(v[0]/v[1]))
            else:
                r.append(p)

            v = V((i,1)) * mu
            if v[1]:
                u.append(int(v[0]/v[1]))
            else:
                u.append(p)

        v = V((1,0)) * mr
        if v[1]:
            r.append(int(v[0]/v[1]))
        else:
            r.append(p)

        v = V((1,0)) * mu
        if v[1]:
            u.append(int(v[0]/v[1]))
        else:
            u.append(p)

        o = Origami(r,u,
            as_tuple=True,
            positions=positions,
            name="Projective line origami on GF(%d)\n r = %s    u = %s\n     %s        %s"%(p,sr[0],su[0],sr[1],su[1]))
        return o

    @staticmethod
    def Heisenberg(p):
        r"""
        Return the Heisenberg origami.

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: h2 = origamis.Heisenberg(2)
            sage: h2.stratum_component()
            H_3(1^4)^c
            sage: h2.veech_group().index()
            3
            sage: h2.sum_of_lyapunov_exponents()
            2

            sage: h3 = origamis.Heisenberg(3)
            sage: h3.stratum_component()
            H_10(2^9)^even
            sage: h3.veech_group().index()
            1
            sage: h3.sum_of_lyapunov_exponents()
            5

            sage: h5 = origamis.Heisenberg(5)
            sage: h5.stratum_component()
            H_51(4^25)^even
            sage: h5.veech_group().index()
            1
            sage: h5.sum_of_lyapunov_exponents()
            15
        """
        p = ZZ(p)
        if not p.is_prime():
            raise ValueError("p (={}) must be prime".format(p))
        N = p**3  # map (a,b,d) -> a*p**2 + b*p + d
        pp = p*p
        r = [None]*N
        u = [None]*N
        for a in xrange(p):
            for b in xrange(p):
                for d in xrange(p):
                    n = a*pp + b*p + d
                    r[n] = ((a+1)%p)*pp + b*p + d
                    u[n] = a*pp + ((b+1)%p)*p + ((a+d)%p)
        return Origami(r,u,
                as_tuple=True,
                name="Heisenberg origami on GF(%d)"%p)

origamis = OrigamiGenerators()
