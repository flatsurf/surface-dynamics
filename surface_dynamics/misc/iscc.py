#*****************************************************************************
#       Copyright (C) 2019 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from sage.cpython.string import bytes_to_str, str_to_bytes

def iscc_rg_string(r):
    r"""
    EXAMPLES::

        sage: from surface_dynamics import *
        sage: from surface_dynamics.misc.iscc import iscc_rg_string
        sage: r = RibbonGraph(faces='(0,2,4,6,7)(8,9,5,3,1)', edges='(0,1)(2,3)(4,5)(6,7)(8,9)')
        sage: iscc_rg_string(r)
        '[b0, b1] -> { [l0, l1, l2, l3, l4, l5, l6, l7, l8, l9] : l0 >= 0 and l1 >= 0 and l2 >= 0 and l3 >= 0 and l4 >= 0 and l5 >= 0 and l6 >= 0 and l7 >= 0 and l8 >= 0 and l9 >= 0 and l0 = l1 and l2 = l3 and l4 = l5 and l6 = l7 and l8 = l9 and b0 = l0 + l2 + l4 + l6 + l7 and b1 = l1 + l8 + l9 + l5 + l3 }'
    """
    boundary_vars = ', '.join('b%d'%i for i in range(r.num_faces()))
    lengths_vars = ', '.join('l%d'%i for i in r.darts())
    
    ieqs = ' and '.join('l%d >= 0' % i for i in r.darts())
    edges = ' and '.join('l%d = l%d' %(i,j) for i,j in r.edges())
    
    eqns = []
    for i,f in enumerate(r.faces()):
        eqns.append('b%d = %s' % (i, ' + '.join('l%d' %j for j in f)))
    eqns = ' and '.join(eqns)
    
    return "[{boundary_vars}] -> {{ [{lengths_vars}] : {ieqs} and {edges} and {eqns} }}".format(
               boundary_vars=boundary_vars,
               lengths_vars=lengths_vars,
               edges=edges,
               ieqs=ieqs,
               eqns=eqns)

def iscc_cd_string(cd):
    r"""
    Given a cylinder diagram produces a parametrized polytopes (by the widths
    of cylinders)

    EXAMPLES::

        sage: from surface_dynamics.misc.iscc import iscc_cd_string
        sage: from surface_dynamics import CylinderDiagram
        sage: cd = CylinderDiagram('(0,1)-(0,2) (2)-(1)')
        sage: iscc_cd_string(cd)
        '[w0, w1] -> { [l0, l1, l2] : l0 >= 0 and l1 >= 0 and l2 >= 0 and w0 = l0 + l1 = l0 + l2 and w1 = l2 = l1 }'
    """
    widths_vars = ', '.join('w%d'%i for i in range(cd.ncyls()))
    lengths_vars = ', '.join('l%d'%i for i in range(cd.num_edges()))
    
    ieqs = []
    for i in range(cd.num_edges()):
        ieqs.append('l%d >= 0' % i)
    ieqs = ' and '.join(ieqs)
    
    eqns = []
    for i,(top,bot) in enumerate(cd.cylinders()):
        eqns.append('w%d = %s = %s' % (i,
                        ' + '.join('l%d' %j for j in top),
                        ' + '.join('l%d' %j for j in bot)))
    eqns = ' and '.join(eqns)
    
    return "[{widths_vars}] -> {{ [{lengths_vars}] : {ieqs} and {eqns} }}".format(
               widths_vars = widths_vars,
               lengths_vars=lengths_vars,
               ieqs=ieqs,
               eqns=eqns)

def iscc_card(arg):
    r"""
    EXAMPLES::

        sage: from surface_dynamics import *
        sage: from surface_dynamics.misc.iscc import iscc_card
        sage: cd = CylinderDiagram('(0,1)-(0,2) (2)-(1)')
        sage: iscc_card(cd)  # optional - barvinok
        '[w0, w1] -> { 1 : 0 <= w1 <= w0 }'

        sage: r = RibbonGraph(faces='(0,2,4,6,7)(8,9,5,3,1)', edges='(0,1)(2,3)(4,5)(6,7)(8,9)')
        sage: iscc_card(r)   # optional - barvinok
        '[b0, b1] -> { (((1 - 1/8 * b0) + 17/12 * b1 + 5/8 * b1^2 + 1/12 * b1^3) + 1/4 * floor((2b0)/4)) : (b0 + b1) mod 2 = 0 and 0 <= b1 <= -2 + b0; ((1 + 31/24 * b0 + 5/8 * b0^2 + 1/12 * b0^3) + 1/4 * floor((2b0)/4)) : (b0 + b1) mod 2 = 0 and b0 >= 0 and b1 >= b0 }'
    """
    from surface_dynamics.flat_surfaces.separatrix_diagram import CylinderDiagram
    from surface_dynamics.flat_surfaces.homology import RibbonGraph

    if isinstance(arg, CylinderDiagram):
        poly = iscc_cd_string(arg)
    elif isinstance(arg, RibbonGraph):
        poly = iscc_rg_string(arg)

    from subprocess import Popen, PIPE
    
    cmd = 'P := {};\n'.format(poly)
    cmd += 'q := card P;\n'
    cmd += 'q;'
    
    proc = Popen(['iscc'], stdin=PIPE, stdout=PIPE, stderr=PIPE)
    ans, err = proc.communicate(str_to_bytes(cmd))
    
    return bytes_to_str(ans).strip()

def parse_iscc_result(ans, d):
    R = PolynomialRing(QQ, ['w%d'%i for i in range(d)], d)
    gens = {'w%d'%i: R.gen(i) for i in range(d)}
    prefix = '[' + ', '.join('w%d'%i for i in range(d)) + ']'
    prefix += ' -> { '
    suffix = ' }\n'
    if not ans.startswith(prefix):
        raise ValueError('bad prefix')
    if not ans.endswith(suffix):
        raise ValueError('bad suffix')
    pieces = ans[len(prefix):len(ans)-len(suffix)].split(';')
    res = []
    for piece in pieces:
        if piece.count(':') != 1:
            raise ValueError('parsing error!')
        i = piece.find(':')
        poly = R(sage_eval(piece[:i], locals=gens))
        cond = piece[i+1:]
        res.append((cond,poly))
    return res

def higher_order_terms(p):
    R = p.parent()
    p_trunc = R.zero()
    deg = p.degree()
    for coeff, exp in zip(p.coefficients(), p.exponents()):
        if sum(exp) == deg:
            p_trunc += R({exp:coeff})
    return p_trunc
