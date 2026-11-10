# ****************************************************************************
#       Copyright (C) 2026 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from collections import defaultdict
from sage.all import DiGraph


def dynamic_cycle_cover(G):
    r"""
    Return a cycle cover of the multigraph ``G``.

    The algorithm is greedy and not guaranteed to be minimal.

    EXAMPLES::

        sage: from surface_dynamics.misc.cycle_cover import dynamic_cycle_cover
        sage: G = DiGraph(multiedges=True, loops=True)
        sage: G.add_edge(0, 0, 'a')
        sage: G.add_edge(0, 1, 'b')
        sage: G.add_edge(0, 1, 'c')
        sage: G.add_edge(1, 0, 'd')
        sage: dynamic_cycle_cover(G)  # random
        [[(0, 0, 'a')], [(0, 1, 'c'), (1, 0, 'd')], [(0, 1, 'b'), (1, 0, 'd')]]
    """
    if not G.is_strongly_connected():
        raise ValueError("the input graph must be strongly connected")

    cycle_list = []
    endpoints_to_labels = defaultdict(list)
    H = DiGraph()
    H.add_vertices(G.vertices())
    todo = defaultdict(int)
    for u_v_label in G.edges():
        u, v, label = u_v_label
        if u == v:
            cycle_list.append([u_v_label])
            continue
        endpoints_to_labels[u,v].append(label)
        H.add_edge(u, v, 1)
        todo[u_v_label] += 1

    m = H.num_edges()

    while todo:
        u, v, label = next(iter(todo))
        assert H.edge_label(u, v) == 1

        vertex_cycle = [u] + H.shortest_path(v, u, by_weight=True)
        cycle = []
        for i in range(len(vertex_cycle) - 1):
            u0 = vertex_cycle[i]
            u1 = vertex_cycle[i + 1]
            labels = endpoints_to_labels[u0, u1]
            if len(labels) == 1:
                label, = labels
                H.set_edge_label(u0, u1, m + 1)
            else:
                label = labels.pop()
            edge = (u0, u1, label)
            if edge in todo:
                todo[edge] -= 1
                if not todo[edge]:
                    del todo[edge]
            cycle.append(edge)

        cycle_list.append(cycle)

    return cycle_list
