#!/usr/bin/env python3

import copy


def main():
    k = int(input())

    patterns = form_binary_patterns(k)
    result = circular_string(patterns)
    print(result)


def _find_euler_cycle(graph_: list, start_vertex_i: int) -> list:
    """
    Find an eulerian path in the given 'graph'
    :param graph_: the graph object to examine. Graph is a list of tuples, whose [0] element is a list of indexes of
    vertices connected to the given one by an edge
    :param start_vertex_i: the index of a vertex to start eulerian path from
    :return: a list of indexes of vertices forming an eulerian path

    Note this function does NOT check balancity (eulericity) of the given graph
    """
    graph = []
    for vertex in graph_:
        graph.append(copy.copy(vertex[0]))

    cycle = [start_vertex_i]

    vertices_to_visit = [0]  # Indexes of vertices in a 'cycle' list
    while len(vertices_to_visit) > 0:
        current_cycle_start_i = vertices_to_visit.pop(0)
        current_start_vertex_i = cycle[current_cycle_start_i]
        if len(graph[current_start_vertex_i]) == 0:
            continue
        current_cycle = []  # Start vertex is already in the cycle

        current_vertex_i = current_start_vertex_i
        current_edges = graph[current_vertex_i]
        while len(current_edges) > 0:
            current_vertex_i = current_edges.pop(0)
            if len(current_edges) > 0:
                vertices_to_visit.append(current_cycle_start_i + len(current_cycle))
            current_cycle.append(current_vertex_i)
            current_edges = graph[current_vertex_i]

        for v in current_cycle[::-1]:
            cycle.insert(current_cycle_start_i + 1, v)

    return cycle


def circular_string(patterns: list) -> str:
    """
    Reconstruct the DNA from its k-mers (in the given list 'patterns')
    """
    # A graph representation. Each vertex is a tuple: ([], ...); the list at 0 position is a list of connected vertices
    graph = []
    # Map of graph vertices (keys) to graph index
    graph_map = {}

    # Create a graph
    for pattern in patterns:
        v_from = pattern[:-1]
        if v_from not in graph_map.keys():
            graph.append(([], v_from))
            graph_map[v_from] = len(graph) - 1
        v_from_index = graph_map[v_from]

        v_to = pattern[1:]
        if v_to not in graph_map.keys():
            graph.append(([], v_to))
            graph_map[v_to] = len(graph) - 1
        v_to_index = graph_map[v_to]

        graph[v_from_index][0].append(v_to_index)

    # Find an eulerian cycle and an eulerian path
    euler_cycle = _find_euler_cycle(graph, 0)

    # Convert the path found to mers
    # euler_mers = [graph[euler_cycle[0]][1]]
    # for path_i in euler_cycle[1:-1]:
    euler_mers = []
    for path_i in euler_cycle[:-1]:
        current_vertex_mer = graph[path_i][1]
        euler_mers.append(current_vertex_mer[-1])

    # Form the result
    result = ''.join(euler_mers)

    return result


def form_binary_patterns(k: int) -> list:
    """
    Return a list of strings containing all binary numbers of length not exceeding 'k' (with leading zeroes)
    """
    result = []
    format_string = '{{:0{}b}}'.format(k)

    for n in range(2 ** k):
        result.append(format_string.format(n))

    return result


if __name__ == '__main__':
    main()
