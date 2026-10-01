#!/usr/bin/env python
r"""
Check the dynamic cycle cover
"""
# *****************************************************************************
#       Copyright (C) 2026 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************


def check_cycle_cover(graph, cycle_cover):
    # check that `cycle_cover` is made of cycles in `graph`
    # and that it is indeed a cover
    from collections import defaultdict

    edges = defaultdict(int)
    for u_v_label in graph.edges():
        edges[u_v_label] += 1
    for cycle in cycle_cover:
        for edge in cycle:
            assert isinstance(edge, tuple)
            assert len(edge) == 3
            assert graph.has_edge(*edge)
            edges[edge] -= 1
        for i in range(len(cycle)):
            e0 = cycle[i]
            e1 = cycle[(i + 1) % len(cycle)]
            assert e0[1] == e1[0]

    # check that it is indeed a cover
    assert all(value <= 0 for value in edges.values()), (graph, cycle_cover, edges)


def test_cycle_cover():
    from sage.all import DiGraph
    from surface_dynamics.misc.cycle_cover import dynamic_cycle_cover

    Glist = []

    G = DiGraph(multiedges=True, loops=True)
    G.add_edge(0, 1)
    G.add_edge(1, 0)
    G.add_edge(1, 2)
    G.add_edge(1, 2)
    G.add_edge(2, 3)
    G.add_edge(3, 1)
    Glist.append(G)

    G = DiGraph(multiedges=False, loops=False)
    G.add_edge(0, 1)
    G.add_edge(0, 2)
    G.add_edge(0, 3)
    G.add_edge(0, 4)
    G.add_edge(1, 2)
    G.add_edge(1, 3)
    G.add_edge(1, 4)
    G.add_edge(2, 3)
    G.add_edge(2, 4)
    G.add_edge(3, 4)
    G.add_edge(4, 0)
    Glist.append(G)

    for G in Glist:
        cycles = dynamic_cycle_cover(G)
        check_cycle_cover(G, dynamic_cycle_cover(G))
