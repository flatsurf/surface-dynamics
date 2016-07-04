r"""
Marked partition

The suspensions of interval exchanges have a natural geometric structure which
is a translation surface. The singularities of the interval exchanges yield
singularities on the surface which gives an integer partition. As the left and
the right of the interval plays a special role, we consider partitions with
marking. Either the points at the left and the right of the interval coincide in
which case the marking is of type one, otherwise they are two different parts
and the marking is of type two.
"""
#*****************************************************************************
#       Copyright (C) 2012 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sage_object import SageObject
from sage.combinat.partition import Partition
from sage.rings.integer import Integer

class Marking(SageObject):
    r"""
    EXAMPLES::

        sage: from surface_dynamics.interval_exchanges.marked_partition import Marking
        sage: Marking('2o3')
        2o3
        sage: Marking('3|0')
        3|0
        sage: Marking(1,(2,0))
        2|0
        sage: Marking(2,(2,1))
        2o1
    """
    def __init__(self, *args):
        r"""
        TESTS::

            sage: from surface_dynamics.interval_exchanges.marked_partition import Marking
            sage: m = Marking('2o3')
            sage: loads(dumps(m)) == m
            True
        """
        if len(args) == 1 and isinstance(args[0], Marking):
            self.t = args[0].t
            self.data = args[0].data
        elif len(args) == 1 and isinstance(args[0], str):
            tt = args[0]
            if "|" in tt:
                self.t = 1
                self.data = tuple(map(Integer,tt.split('|')))
            elif "o" in tt:
                self.t = 2
                self.data = tuple(map(Integer,tt.split('o')))
        elif len(args) == 1 and isinstance(args[0], (list,tuple)) and len(args[0]) == 2:
            self.t = Integer(args[0][0])
            self.data = tuple(map(Integer,args[0][1]))
        elif len(args) == 2:
            self.t = Integer(args[0])
            self.data = tuple(map(Integer,args[1]))
        else:
            raise ValueError("can not build marking from given data")

        if self.t != 1 and self.t != 2:
            raise ValueError("the first argument should be 1 or 2")
        if len(self.data) != 2:
            raise ValueError("the second argument should have length 2")
        if self.data[0] < 0 or self.data[1] < 0:
            raise ValueError("the second argument should be a 2-tuple of non negative integers")
        if self.t == 1 and self.data[1] > self.data[0]:
            raise ValueError("for type 1, the 2-tuple should have a first argument greater or equal to the second")

    def __hash__(self):
        r"""
        TESTS::

            sage: from surface_dynamics.interval_exchanges.marked_partition import Marking
            sage: hash(Marking('9o4'))  #indirect doctest
            3713072971714925208
        """
        return hash(self.t) + hash(self.data)

    def _repr_(self):
       r"""
       String representation.

       TESTS::

        sage: from surface_dynamics.interval_exchanges.marked_partition import Marking
        sage: Marking(1, (2,0))._repr_()
        '2|0'
        sage: Marking(2, (3,2))._repr_()
        '3o2'
       """
       if self.t == 1:
           return "%d|%d"%self.data
       else:
           return "%do%d"%self.data

    def __eq__(self,other):
        r"""
        TESTS::

            sage: from surface_dynamics.interval_exchanges.marked_partition import Marking
            sage: Marking('2|1') == Marking('2|1')
            True
            sage: Marking('1o2') == Marking('2o1')
            False
        """
        if not isinstance(other,Marking):
            return False
        return self.t == other.t and self.data == other.data

    def __ne__(self,other):
        r"""
        TESTS::

            sage: from surface_dynamics.interval_exchanges.marked_partition import Marking
            sage: Marking('3|2') != Marking('2o3')
            True
        """
        return not self.__eq__(other)

    def left(self):
        r"""
        Return the part that is marked on the left.

        EXAMPLES::

            sage: from surface_dynamics.interval_exchanges.marked_partition import Marking
            sage: Marking(1,(2,0)).left()
            2
            sage: Marking(2,(1,3)).left()
            1
        """
        return self.data[0]

    def right(self):
        r"""
        Return the part that is marked on the right.

        EXAMPLES::

            sage: from surface_dynamics.interval_exchanges.marked_partition import Marking
            sage: Marking(1,(2,0)).right()
            2
            sage: Marking(2,(1,3)).right()
            3
        """
        if self.t == 1: return self.data[0]
        else: return self.data[1]

class MarkedPartition(SageObject):
    r"""
    Marked partition

    A marked partition is an integer partition ``p`` with the additional data of
    either a couple ``(m,a)`` where ``m`` is an element of ``p`` and ``a`` is an
    integer lesser than ``p`` or a couple ``(m_l,m_r)`` where ``m_l`` and
    ``m_r`` are distinct element of ``p``. In the first case the type is ``1``
    while in the second the type is ``2``.

    EXAMPLES::

        sage: from surface_dynamics.interval_exchanges.marked_partition import MarkedPartition
        sage: m1 = MarkedPartition([3,2,1],1,(2,0)); m1
        2|0 [3, 2, 1]
        sage: m2 = MarkedPartition([3,2,1],2,(2,3)); m2
        2o3 [3, 2, 1]

    It also posisble to initialize it from a string::

        sage: MarkedPartition('2|1 [3,2,2]')
        2|1 [3, 2, 2]

    TESTS::

        sage: MarkedPartition(m1)
        2|0 [3, 2, 1]
        sage: MarkedPartition(str(m1))
        2|0 [3, 2, 1]
        sage: MarkedPartition(m2)
        2o3 [3, 2, 1]
        sage: MarkedPartition(str(m2))
        2o3 [3, 2, 1]
    """
    def __init__(self, *args, **kwds):
        r"""
        TESTS::

            sage: from surface_dynamics.interval_exchanges.marked_partition import MarkedPartition
            sage: p = MarkedPartition([3,2,1], 1,(2,0))
            sage: loads(dumps(p)) == p
            True
        """
        from sage.combinat.partition import Partition

        if len(args) == 1 and isinstance(args[0], MarkedPartition):
            from copy import copy
            self.p = copy(args[0].p)
            self.m = copy(args[0].m)
        elif len(args) == 1 and isinstance(args[0], str):
            import re
            x = re.compile("(?P<mark>[^[]*)[^[]*\[(?P<parts>[^]]*)]")
            m = x.match(args[0])
            self.m = Marking(m.groupdict()["mark"])
            self.p = Partition(sorted(map(Integer,m.groupdict()["parts"].split(',')),reverse=True))
        elif len(args) == 1 and isinstance(args[0], (list,tuple)) and len(args[0]) == 3:
            self.p = Partition(sorted(args[0][0],reverse=True))
            self.m = Marking(args[0][1],args[0][2])
        elif len(args) == 3:
            self.p = Partition(sorted(args[0],reverse=True))
            self.m = Marking(args[1],args[2])
        else:
            raise ValueError("can not build marked partition from given data")

        if kwds.get("check",True):
            if self.p == Partition([]):
                if self.m.t != 2 or self.m.data[0] != 0 or self.m.data[1] != 0:
                    raise ValueError("empty partition has only type 2 with data (0,0)")

            elif self.m.t == 1:
                if self.m.data[0] not in self.p:
                    raise ValueError("%d is not an element of parts"%data[0])
            else:
                if ((self.m.data[0] not in self.p or self.m.data[1] not in self.p) or
                    (self.m.data[0] == self.m.data[1] and list(self.p).count(self.m.data[0]) < 2)):
                    raise ValueError("parts do not contains (m_l,m_r) = (%d,%d)"%(self.m.data[0],self.m.data[1]))

    def _repr_(self):
        r"""
        String representation.

        EXAMPLES::

            sage: from surface_dynamics.interval_exchanges.marked_partition import MarkedPartition
            sage: MarkedPartition([3,1],1,(3,2))._repr_()
            '3|2 [3, 1]'
            sage: MarkedPartition([3,1],2,(3,1))._repr_()
            '3o1 [3, 1]'
        """
        return "%s %s" %(self.m,self.p)

    def __eq__(self, other):
        r"""
        Equality test.

        TESTS::

            sage: from surface_dynamics.interval_exchanges.marked_partition import MarkedPartition
            sage: MarkedPartition([3,2,1],1,(2,0)) == MarkedPartition([3,2,1],1,(2,0))
            True
            sage: MarkedPartition([3,2,1],1,(2,0)) == MarkedPartition([3,2,1],1,(2,1))
            False
        """
        if not isinstance(other,MarkedPartition):
            return False

        return self.m == other.m and self.p == other.p

    def __ne__(self,other):
        r"""
        Difference test.

        TESTS::

            sage: from surface_dynamics.interval_exchanges.marked_partition import MarkedPartition
            sage: MarkedPartition([3,2,1],1,(2,0)) != MarkedPartition([3,2,1],1,(2,0))
            False
            sage: MarkedPartition([3,2,1],1,(2,0)) != MarkedPartition([3,2,1],1,(2,1))
            True
        """
        return not self.__eq__(other)

    def __hash__(self):
        r"""
        TESTS::

            sage: from surface_dynamics.interval_exchanges.marked_partition import MarkedPartition
            sage: MarkedPartition([3,1],1,(3,2)).__hash__()
            7426167593987224238
        """
        return hash(self.m) + hash(self.p)

    def left(self):
        r"""
        Return the part that is marked on the left.

        EXAMPLES::

            sage: from surface_dynamics.interval_exchanges.marked_partition import MarkedPartition
            sage: MarkedPartition([3,2,1],1,(2,0)).left()
            2
            sage: MarkedPartition([3,2,1],2,(1,3)).left()
            1
        """
        return self.m.left()

    def right(self):
        r"""
        Return the part that is marked on the right.

        EXAMPLES::


            sage: from surface_dynamics.interval_exchanges.marked_partition import MarkedPartition
            sage: MarkedPartition([3,2,1],1,(2,0)).right()
            2
            sage: MarkedPartition([3,2,1],2,(1,3)).right()
            3
        """
        return self.m.right()

    def partition(self):
        r"""
        Returns the underlying partition.

        EXAMPLES::

            sage: from surface_dynamics.interval_exchanges.marked_partition import MarkedPartition
            sage: MarkedPartition([1,3,2],1,(2,0)).partition()
            [3, 2, 1]
            sage: MarkedPartition([2,3,1],2,(1,3)).partition()
            [3, 2, 1]
        """
        return self.p

    def marking(self):
        r"""
        Return a 3-tuple ``(type,d0,d1)`` where ``type`` is the type of the
        marking and ``(d0,d1)`` corresponds to the data associated to the
        marking.

        EXAMPLES::

            sage: from surface_dynamics.interval_exchanges.marked_partition import MarkedPartition
            sage: MarkedPartition([3,2,1],1,(2,0)).marking()
            2|0
        """
        return self.m

    def is_odd(self):
        r"""
        Returns True if all terms of p are odd

        EXAMPLES::

            sage: from surface_dynamics.interval_exchanges.marked_partition import MarkedPartition
            sage: MarkedPartition([3,1,1],1,(3,1)).is_odd()
            True
            sage: MarkedPartition([3,2,2],1,(3,1)).is_odd()
            False
        """
        return all(k%2 for k in self.p)
