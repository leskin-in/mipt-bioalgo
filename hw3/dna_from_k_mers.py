#!/usr/bin/env python3

import sys
import copy


def main():
    k = None
    patterns = []
    for l in sys.stdin:
        if k is None:
            k = int(l[:-1])
            continue
        patterns.append(l[:-1])

    result = reconstruct_dna(patterns)
    print(result)


def _check_balancity_and_turn_to_euler_cycle(graph: list) -> tuple:
    """
    Check the graph is balanced (except for exactly two vertices). If so, add an edge connecting the two vertices that
    make the graph unbalanced. If the graph is unbalanced, throws an exception.
    :returns: start_vertex_index, end_vertex_index: indexes of path start and path end vertices
    """
    # A map: vertex ID -> balance of incoming and outgoing edges (outgoing = -1, incoming = +1)
    vertex_balances = {}
    for i in range(len(graph)):
        vertex_balances[i] = 0

    for vertex_i in range(len(graph)):
        vertex = graph[vertex_i]
        vertex_balances[vertex_i] -= len(vertex[0])
        for vertex_to in vertex[0]:
            vertex_balances[vertex_to] += 1

    start_vertex_i = None
    end_vertex_i = None
    for vertex_i in range(len(graph)):
        if vertex_balances[vertex_i] < 0:
            if start_vertex_i is not None:
                raise Exception(
                    "The given graph is not balanced. Two vertices with negative balances: {} and {}".format(
                        graph[start_vertex_i], graph[vertex_i]
                    )
                )
            start_vertex_i = vertex_i
        elif vertex_balances[vertex_i] > 0:
            if end_vertex_i is not None:
                raise Exception(
                    "The given graph is not balanced. Two vertices with positive balances: {} and {}".format(
                        graph[end_vertex_i], graph[vertex_i]
                    )
                )
            end_vertex_i = vertex_i
    if start_vertex_i is None:
        raise Exception("No start vertex found")
    if end_vertex_i is None:
        raise Exception("No end vertex found")

    graph[end_vertex_i][0].append(start_vertex_i)

    return start_vertex_i, end_vertex_i


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


def _euler_cycle_to_path(cycle: list, cycling_edge: tuple) -> list:
    """
    Convert an eulerian cycle into a path which starts at 'cycling_edge[1]' and ends at 'cycling_edge[0]'
    """
    for i in range(len(cycle) - 1):
        start_i = cycle[i]
        end_i = cycle[i + 1]
        if start_i == cycling_edge[1] and end_i == cycling_edge[0]:
            result = cycle[i + 1:]
            result.pop()  # Remove the "connecting" last vertex
            result.extend(cycle[:i + 1])
            return result

    raise Exception("The given cycling edge {} was not found in the given cycle {}".format(cycling_edge, cycle))


def reconstruct_dna(patterns: list) -> str:
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
    cycling_edge = _check_balancity_and_turn_to_euler_cycle(graph)
    euler_cycle = _find_euler_cycle(graph, cycling_edge[0])
    euler_path = _euler_cycle_to_path(euler_cycle, cycling_edge)

    # Convert the path found to mers
    euler_mers = [graph[euler_path[0]][1]]
    for path_i in euler_path[1:]:
        current_vertex_mer = graph[path_i][1]
        euler_mers.append(current_vertex_mer[-1])

    # Form the result
    result = ''.join(euler_mers)

    return result


if __name__ == '__main__':
    main()
