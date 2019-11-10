#!/usr/bin/env python3


def main():
    peptide1 = list(input())
    peptide2 = list(input())

    graph = form_graph(peptide1, peptide2)
    paths = ford_bellman_search(graph, 0)
    path = form_peptide_path(paths)
    peptide1_r, peptide2_r = path_peptides(peptide1, peptide2, graph, path)

    print(-paths[-1][0], ''.join(peptide1_r), ''.join(peptide2_r), sep='\n')


# A scoring matrix. Insertion-deletion is penalized using a dedicated constant
SCORE_MATRIX = (
    (4, 0, -2, -1, -2, 0, -2, -1, -1, -1, -1, -2, -1, -1, -1, 1, 0, 0, -3, -2),
    (0, 9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2),
    (-2, -3, 6, 2, -3, -1, -1, -3, -1, -4, -3, 1, -1, 0, -2, 0, -1, -3, -4, -3),
    (-1, -4, 2, 5, -3, -2, 0, -3, 1, -3, -2, 0, -1, 2, 0, 0, -1, -2, -3, -2),
    (-2, -2, -3, -3, 6, -3, -1, 0, -3, 0, 0, -3, -4, -3, -3, -2, -2, -1, 1, 3),
    (0, -3, -1, -2, -3, 6, -2, -4, -2, -4, -3, 0, -2, -2, -2, 0, -2, -3, -2, -3),
    (-2, -3, -1, 0, -1, -2, 8, -3, -1, -3, -2, 1, -2, 0, 0, -1, -2, -3, -2, 2),
    (-1, -1, -3, -3, 0, -4, -3, 4, -3, 2, 1, -3, -3, -3, -3, -2, -1, 3, -3, -1),
    (-1, -3, -1, 1, -3, -2, -1, -3, 5, -2, -1, 0, -1, 1, 2, 0, -1, -2, -3, -2),
    (-1, -1, -4, -3, 0, -4, -3, 2, -2, 4, 2, -3, -3, -2, -2, -2, -1, 1, -2, -1),
    (-1, -1, -3, -2, 0, -3, -2, 1, -1, 2, 5, -2, -2, 0, -1, -1, -1, 1, -1, -1),
    (-2, -3, 1, 0, -3, 0, 1, -3, 0, -3, -2, 6, -2, 0, 0, 1, 0, -3, -4, -2),
    (-1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2, 7, -1, -2, -1, -1, -2, -4, -3),
    (-1, -3, 0, 2, -3, -2, 0, -3, 1, -2, 0, 0, -1, 5, 1, 0, -1, -2, -2, -1),
    (-1, -3, -2, 0, -3, -2, 0, -3, 2, -2, -1, 0, -2, 1, 5, -1, -1, -3, -3, -2),
    (1, -1, 0, 0, -2, 0, -1, -2, 0, -2, -1, 1, -1, 0, -1, 4, 1, -2, -3, -2),
    (0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1, 0, -1, -1, -1, 1, 5, 0, -2, -2),
    (0, -1, -3, -2, -1, -3, -3, 3, -2, 1, 1, -3, -2, -2, -3, -2, 0, 4, -3, -1),
    (-3, -2, -4, -3, 1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11, 2),
    (-2, -2, -3, -2, 3, -3, 2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1, 2, 7)
)


# A mapping to encode amino acid as an index in the SCORE_MATRIX
AMINO_TO_SCORE_INDEX = {
    'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4,
    'G': 5, 'H': 6, 'I': 7, 'K': 8, 'L': 9,
    'M': 10, 'N': 11, 'P': 12, 'Q': 13, 'R': 14,
    'S': 15, 'T': 16, 'V': 17, 'W': 18, 'Y': 19
}


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

    def reverse_vi(self, index: int) -> tuple:
        """
        :return: a pair (width, height) of the vertex at the given index
        """
        return index % self.w, index // self.w

    def __str__(self):
        result_h = 'Graph({}, {})'.format(self.w, self.h)
        result_g = '\n'.join([' '.join([str(self.g[self.vi(w, h)]) for w in range(self.w)]) for h in range(self.h)])
        return result_h + '\n' + result_g + '\n'

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
        # same_letter_p2_score = SCORE_MATRIX[p2_index][p2_index]
        for w_i in range(graph.w - 1):
            p1_index = AMINO_TO_SCORE_INDEX[peptide1[w_i]]
            # same_letter_p1_score = SCORE_MATRIX[p1_index][p1_index]
            diff_letter_score = SCORE_MATRIX[p1_index][p2_index]

            # Skips
            graph.g[graph.vi(w_i, h_i)].adj.append(AdjVertex(graph.vi(w_i + 1, h_i), INDEL_PENALTY))  # same_letter_p1_score
            graph.g[graph.vi(w_i, h_i)].adj.append(AdjVertex(graph.vi(w_i, h_i + 1), INDEL_PENALTY))  # same_letter_p2_score
            # Both letters
            graph.g[graph.vi(w_i, h_i)].adj.append(AdjVertex(graph.vi(w_i + 1, h_i + 1), -diff_letter_score))

    return graph


def ford_bellman_search(gr: Graph, v_start: int) -> list:
    """
    Find the SHORTEST paths from 'v_start' to all vertices in the given 'graph'.

    :param gr: a graph to examine. Must have no cycles
    :param v_start: a vertex index to start search from
    :return: a list of path lengths. The vertices' layout is the same as in the original 'graph'
    """
    paths = [None for _ in range(gr.h * gr.w)]
    paths[v_start] = (0, None)

    paths_undergone_modification = True
    while paths_undergone_modification:
        paths_undergone_modification = False
        for v_start in range(len(gr.g)):
            if paths[v_start] is None:
                continue
            curr_dist = paths[v_start][0]
            for edge in gr.g[v_start].adj:
                edge_end = edge.i
                edge_dist = curr_dist + edge.d
                if paths[edge_end] is None or edge_dist < paths[edge_end][0]:
                    paths[edge_end] = (edge_dist, v_start)
                    paths_undergone_modification = True

    return paths


def form_peptide_path(paths: list) -> list:
    """
    Form a path of peptides using 'paths' obtained from 'dijkstra_search'
    """
    result = []
    curr_v = paths[-1][1]
    while curr_v is not None:
        result.append(curr_v)
        curr_v = paths[curr_v][1]
    return result


def _resolve_path(a: tuple, b: tuple) -> str:
    """
    Resolve the path from 'a' to 'b' and give a direction on it
    """
    if a[0] == b[0] and a[1] + 1 == b[1]:
        return 'h'
    if a[0] + 1 == b[0] and b[1] == a[1]:
        return 'w'
    if a[0] + 1 == b[0] and a[1] + 1 == b[1]:
        return 'b'
    raise Exception("Unexpected path")


# A spacing symbol
SPACING = '-'


def path_peptides(peptide_w: list, peptide_h: list, gr: Graph, path: list) -> tuple:
    """
    Patch peptides using the given 'path'
    :param peptide_w: the peptide that forms the "width" of the path matrix
    :param peptide_h: the peptide that forms the "height" of the path matrix
    :param gr: a graph (used only to calculate vertices' indexes properly)
    :param path: a prepared path (not a whole path matrix returned by 'dijkstra_search') starting at the last amino
    :return: patched 'peptide_w' and 'peptide_h'
    """
    # Form a list of directions
    curr_v = gr.vi(gr.w - 1, gr.h - 1)
    prev_v = None
    directions = []
    for path_i in range(len(path)):
        prev_v = curr_v
        curr_v = path[path_i]
        if curr_v is None:
            break
        prev_v_real = gr.reverse_vi(prev_v)
        curr_v_real = gr.reverse_vi(curr_v)
        directions.insert(0, _resolve_path(curr_v_real, prev_v_real))

    # Follow directions
    peptide_w_r = []
    peptide_w_i = 0
    peptide_h_r = []
    peptide_h_i = 0
    for direction in directions:
        if direction == 'h':
            peptide_w_r.append(SPACING)
            peptide_h_r.append(peptide_h[peptide_h_i])
            peptide_h_i += 1
        elif direction == 'w':
            peptide_w_r.append(peptide_w[peptide_w_i])
            peptide_w_i += 1
            peptide_h_r.append(SPACING)
        elif direction == 'b':
            peptide_w_r.append(peptide_w[peptide_w_i])
            peptide_w_i += 1
            peptide_h_r.append(peptide_h[peptide_h_i])
            peptide_h_i += 1

    return peptide_w_r, peptide_h_r


if __name__ == '__main__':
    main()
