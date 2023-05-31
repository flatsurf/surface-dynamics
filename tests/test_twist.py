import pytest

from surface_dynamics import CylinderDiagram, Twist

def test_find_homologous_cylinders():
    cd = CylinderDiagram('(0,2,1)-(4,5) (3,5)-(0,2,1) (4)-(3)')
    tw = Twist(cd)
    assert set(tw.find_homologous_cylinders()[0]) == {0, 1}
    cd = CylinderDiagram('(0,3)-(0,5) (1,2)-(1,4) (4,6)-(3,7) (5,7)-(2,6)')
    tw = Twist(cd)
    assert set(tw.find_homologous_cylinders()[0]) == {2, 3}
    cd = CylinderDiagram('(0)-(2) (1,2,3)-(4,5) (4)-(3) (5)-(0,1)')
    tw = Twist(cd)
    assert tw.find_homologous_cylinders() == []

def test_twist_dimension():
    cd = CylinderDiagram('(0)-(1) (1,3,4,2)-(5,6) (5)-(0,4) (6)-(2,3)')
    tw = Twist(cd)
    assert tw.twist_dimension() == 3
    cd = CylinderDiagram('(0,2,1)-(5,6) (3,4)-(0,2,1) (5)-(4) (6)-(3)')
    tw = Twist(cd)
    assert tw.twist_dimension() == 2
