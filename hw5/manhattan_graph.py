#!/usr/bin/env python3

import bisect


def main():
    h, w = list(map(int, input().split()))

    matrix_down = []
    for _ in range(h):
        matrix_down.append(list(map(int, input().split())))
    input()
    matrix_right = []
    for _ in range(h + 1):
        matrix_right.append(list(map(int, input().split())))

    graph = form_graph(h, w, matrix_down, matrix_right)
    paths = dijkstra_search(graph, 0)

    print(-paths[graph.vi(w, h)])


class AdjVertex:
    """
    An adjacent vertex (index + distance)
    """
    def __init__(self, index, distance):
        self.i = index
        self.d = distance

    def __str__(self):
        return '{}({})'.format(self.i, self.d)

    def __repr__(self):
        return self.__str__()


class Vertex:
    """
    A Vertex (index + list of adjacent vertices)
    """
    def __init__(self, index=None):
        self.i = index
        self.adj = []

    def __str__(self):
        result = ''
        if self.i is not None:
            result = result + '{}'.format(self.i)
        result = result + '[{}]'.format(', '.join(list(map(str, self.adj))))
        return result

    def __repr__(self):
        return self.__str__()


class Graph:
    """
    A graph representation
    """
    def __init__(self, width: int, height: int):
        """
        Construct a new "square" graph with the given 'width' and 'height'
        """
        self.w = width
        self.h = height
        self.g = [Vertex() for _ in range(self.h * self.w)]

    def vi(self, width: int, height: int) -> int:
        """
        :return: a vertex index for the given 'width' and 'height'
        """
        return height * self.w + width

    def __str__(self):
        result_h = 'Graph({}, {})'.format(self.w, self.h)
        result_g = '\n'.join([' '.join([str(self.g[self.vi(w, h)]) for w in range(self.w)]) for h in range(self.h)])
        return result_h + '\n' + result_g + '\n'

    def __repr__(self):
        return self.__str__()


def form_graph(h: int, w: int, matrix_down: list, matrix_right: list) -> Graph:
    """
    Form a graph from a "Right-Down matrix" representation
    """
    graph = Graph(w + 1, h + 1)

    for w_index in range(graph.w):
        for h_index in range(graph.h - 1):
            graph.g[graph.vi(w_index, h_index)].adj.append(
                AdjVertex(graph.vi(w_index, h_index + 1), -matrix_down[h_index][w_index])
            )

    for h_index in range(graph.h):
        for w_index in range(graph.w - 1):
            graph.g[graph.vi(w_index, h_index)].adj.append(
                AdjVertex(graph.vi(w_index + 1, h_index), -matrix_right[h_index][w_index])
            )

    return graph


def dijkstra_search(gr: Graph, v_start: int) -> list:
    """
    Find the SHORTEST paths from 'v_start' to all vertices in the given 'graph'
    :param gr: a graph to examine
    :param v_start: a vertex index to start search from
    :return: a list of path lengths. The vertices' layout is the same as in the original 'graph'

    The implicit assumption of this algorithm is that the graph is acyclic
    """
    paths = [None for _ in range(gr.h * gr.w)]
    paths[v_start] = 0

    examine_list = [(0, v_start), ]

    while len(examine_list) > 0:
        _, v_curr = examine_list.pop(0)
        dist_curr = paths[v_curr]
        for v_next in gr.g[v_curr].adj:
            if paths[v_next.i] is None:
                paths[v_next.i] = dist_curr + v_next.d
            else:
                paths[v_next.i] = min(paths[v_next.i], dist_curr + v_next.d)
            bisect.insort_right(examine_list, (dist_curr + v_next.d, v_next.i))

    return paths


if __name__ == '__main__':
    main()
