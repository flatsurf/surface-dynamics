import pytest

from surface_dynamics import CylinderDiagram
from surface_dynamics.flat_surfaces.twist_space import TwistSpace

def test_homologous_cylinders():
    cd = CylinderDiagram('(0,2,1)-(4,5) (3,5)-(0,2,1) (4)-(3)')
    assert cd.homologous_cylinders() == [[0, 1]]
    cd = CylinderDiagram('(0,3)-(0,5) (1,2)-(1,4) (4,6)-(3,7) (5,7)-(2,6)')
    assert cd.homologous_cylinders() == [[2, 3]]
    cd = CylinderDiagram('(0)-(2) (1,2,3)-(4,5) (4)-(3) (5)-(0,1)')
    assert cd.homologous_cylinders() == []

def test_cylinder_dimension():
    cd = CylinderDiagram('(0)-(1) (1,3,4,2)-(5,6) (5)-(0,4) (6)-(2,3)')
    tw = TwistSpace(cd)
    assert tw.cylinder_dimension() == 3
    cd = CylinderDiagram('(0,2,1)-(5,6) (3,4)-(0,2,1) (5)-(4) (6)-(3)')
    tw = TwistSpace(cd)
    assert tw.cylinder_dimension() == 2
