#!/usr/bin/env python3


def main():
    peptide1 = list(input())
    peptide2 = list(input())

    graph = form_graph(peptide1, peptide2)
    ford_bellman_search(graph)
    path = form_peptide_path(graph)
    peptide1_r, peptide2_r = path_peptides(peptide1, peptide2, path, graph)

    print(-graph.v[graph.vi(graph.w - 1, graph.h - 1)].d, ''.join(peptide1_r), ''.join(peptide2_r), sep='\n')


# A scoring matrix. Insertion-deletion is penalized using a dedicated constant
SCORE_MATRIX = (
    (2, -2, 0, 0, -3, 1, -1, -1, -1, -2, -1, 0, 1, 0, -2, 1, 1, 0, -6, -3),
    (-2, 12, -5, -5, -4, -3, -3, -2, -5, -6, -5, -4, -3, -5, -4, 0, -2, -2, -8, 0),
    (0, -5, 4, 3, -6, 1, 1, -2, 0, -4, -3, 2, -1, 2, -1, 0, 0, -2, -7, -4),
    (0, -5, 3, 4, -5, 0, 1, -2, 0, -3, -2, 1, -1, 2, -1, 0, 0, -2, -7, -4),
    (-3, -4, -6, -5, 9, -5, -2, 1, -5, 2, 0, -3, -5, -5, -4, -3, -3, -1, 0, 7),
    (1, -3, 1, 0, -5, 5, -2, -3, -2, -4, -3, 0, 0, -1, -3, 1, 0, -1, -7, -5),
    (-1, -3, 1, 1, -2, -2, 6, -2, 0, -2, -2, 2, 0, 3, 2, -1, -1, -2, -3, 0),
    (-1, -2, -2, -2, 1, -3, -2, 5, -2, 2, 2, -2, -2, -2, -2, -1, 0, 4, -5, -1),
    (-1, -5, 0, 0, -5, -2, 0, -2, 5, -3, 0, 1, -1, 1, 3, 0, 0, -2, -3, -4),
    (-2, -6, -4, -3, 2, -4, -2, 2, -3, 6, 4, -3, -3, -2, -3, -3, -2, 2, -2, -1),
    (-1, -5, -3, -2, 0, -3, -2, 2, 0, 4, 6, -2, -2, -1, 0, -2, -1, 2, -4, -2),
    (0, -4, 2, 1, -3, 0, 2, -2, 1, -3, -2, 2, 0, 1, 0, 1, 0, -2, -4, -2),
    (1, -3, -1, -1, -5, 0, 0, -2, -1, -3, -2, 0, 6, 0, 0, 1, 0, -1, -6, -5),
    (0, -5, 2, 2, -5, -1, 3, -2, 1, -2, -1, 1, 0, 4, 1, -1, -1, -2, -5, -4),
    (-2, -4, -1, -1, -4, -3, 2, -2, 3, -3, 0, 0, 0, 1, 6, 0, -1, -2, 2, -4),
    (1, 0, 0, 0, -3, 1, -1, -1, 0, -3, -2, 1, 1, -1, 0, 2, 1, -1, -2, -3),
    (1, -2, 0, 0, -3, 0, -1, 0, 0, -2, -1, 0, 0, -1, -1, 1, 3, 0, -5, -3),
    (0, -2, -2, -2, -1, -1, -2, 4, -2, 2, 2, -2, -1, -2, -2, -1, 0, 4, -6, -2),
    (-6, -8, -7, -7, 0, -7, -3, -5, -3, -2, -4, -4, -6, -5, 2, -2, -5, -6, 17, 0),
    (-3, 0, -4, -4, 7, -5, 0, -1, -4, -1, -2, -2, -5, -4, -4, -3, -3, -2, 0, 10),
)


# A mapping to encode amino acid as an index in the SCORE_MATRIX
AMINO_TO_SCORE_INDEX = {
    'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4,
    'G': 5, 'H': 6, 'I': 7, 'K': 8, 'L': 9,
    'M': 10, 'N': 11, 'P': 12, 'Q': 13, 'R': 14,
    'S': 15, 'T': 16, 'V': 17, 'W': 18, 'Y': 19
}


class Edge:
    """
    An edge of a graph
    """
    def __init__(self, f: int, t: int, d: int, mark=None):
        self.f = f
        self.t = t
        self.d = d
        self.mark = mark

    def __str__(self):
        return 'E[{}>({})>{}{}]'.format(
            self.f, self.d, self.t,
            " '{}'".format(str(self.mark)) if self.mark is not None else ''
        )

    def __repr__(self):
        return self.__str__()


class Vertex:
    """
    A Vertex of a graph
    """
    def __init__(self):
        self.d = None
        self.p_edge = None

    def __str__(self):
        return 'V[d{}, {}]'.format(self.d, self.p_edge)

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
        self.v = [Vertex() for _ in range(self.w * self.h)]
        self.e = []

    def vi(self, width: int, height: int) -> int:
        """
        :return: a vertex index for the given 'width' and 'height'
        """
        return height * self.w + width

    def reverse_vi(self, index: int) -> tuple:
        """
        :return: a pair (width, height) of the vertex at the given index
        """
        return index % self.w, index // self.w

    def __str__(self):
        result_h = 'Graph({}, {})'.format(self.w, self.h)
        return result_h + '\n' + str(self.v) + '\n' + str(self.e)

    def __repr__(self):
        return self.__str__()


# The insertion-deletion penalty
INDEL_PENALTY = 5


def form_graph(peptide1: list, peptide2: list) -> Graph:
    """
    Form an alignment-lookup graph from two peptides
    """
    graph = Graph(len(peptide1) + 1, len(peptide2) + 1)

    for h_i in range(graph.h - 1):
        p2_index = AMINO_TO_SCORE_INDEX[peptide2[h_i]]
        for w_i in range(graph.w - 1):
            p1_index = AMINO_TO_SCORE_INDEX[peptide1[w_i]]
            diff_letter_score = SCORE_MATRIX[p1_index][p2_index]

            # Skips
            graph.e.append(Edge(graph.vi(w_i, h_i), graph.vi(w_i + 1, h_i), INDEL_PENALTY, 'w'))
            graph.e.append(Edge(graph.vi(w_i, h_i), graph.vi(w_i, h_i + 1), INDEL_PENALTY, 'h'))
            # Both
            graph.e.append(Edge(graph.vi(w_i, h_i), graph.vi(w_i + 1, h_i + 1), -diff_letter_score, 'b'))
            # Jumps
            graph.e.append(Edge(graph.vi(0, 0), graph.vi(w_i, h_i), 0, 's'))
            graph.e.append(Edge(graph.vi(w_i, h_i), graph.vi(graph.w - 1, graph.h - 1), 0, 'e'))

    graph.v[0].d = 0

    return graph


def ford_bellman_search(graph: Graph):
    """
    Find the SHORTEST paths from a single vertex in 'graph' with the preset distance
    :param graph: a graph to examine. Must meet the following requirements:
        * No cycles
        * 'graph.v[i].d is None' for any i except for exactly one (for that vertex, d must be set to some value)
    As a result, the graph vertices are set with proper distances and parent vertices
    """
    paths_undergone_modification = True
    while paths_undergone_modification:
        paths_undergone_modification = False
        for edge in graph.e:
            if graph.v[edge.f].d is None:
                continue
            new_distance = graph.v[edge.f].d + edge.d
            if graph.v[edge.t].d is None or new_distance < graph.v[edge.t].d:
                paths_undergone_modification = True
                graph.v[edge.t].d = new_distance
                graph.v[edge.t].p_edge = edge


def form_peptide_path(graph: Graph) -> list:
    """
    Form a path of peptides from the given 'graph'
    """
    result = []
    curr_e = graph.v[graph.vi(graph.w - 1, graph.h - 1)].p_edge
    while curr_e is not None:
        result.insert(0, curr_e)
        curr_e = graph.v[curr_e.f].p_edge
    return result


# A spacing symbol
SPACING = '-'


def path_peptides(peptide_w: list, peptide_h: list, path: list, graph: Graph) -> tuple:
    """
    Patch peptides using the given 'path'
    :param peptide_w: the peptide that forms the "width" of the path matrix
    :param peptide_h: the peptide that forms the "height" of the path matrix
    :param path: a path generated by 'form_peptide_path'
    :param graph: a graph to examine (used to calculate width and height of a vertex)
    :return: patched 'peptide_w' and 'peptide_h'
    """
    peptide_w_r = []
    peptide_w_i = 0
    peptide_h_r = []
    peptide_h_i = 0
    for edge in path:
        if edge.mark == 'h':
            peptide_w_r.append(SPACING)
            peptide_h_r.append(peptide_h[peptide_h_i])
            peptide_h_i += 1
        elif edge.mark == 'w':
            peptide_w_r.append(peptide_w[peptide_w_i])
            peptide_w_i += 1
            peptide_h_r.append(SPACING)
        elif edge.mark == 'b':
            peptide_w_r.append(peptide_w[peptide_w_i])
            peptide_w_i += 1
            peptide_h_r.append(peptide_h[peptide_h_i])
            peptide_h_i += 1
        elif edge.mark == 's':
            w, h = graph.reverse_vi(edge.t)
            peptide_w_i += w
            peptide_h_i += h
        elif edge.mark == 'e':
            break

    return peptide_w_r, peptide_h_r


if __name__ == '__main__':
    main()
