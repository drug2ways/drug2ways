# -*- coding: utf-8 -*-

"""Alternative methods for path calculations in graphs (depricated and not used in the package)."""

from typing import List

from networkx import DiGraph


def enumerate_paths(graph, source, targets, lmax, cycle_free=False):
    """Enumerate paths.

    :param graph: graph
    :param source: source node
    :param targets: target nodes
    :param lmax: lmax
    :param cycle_free:
    :return:
    """
    history = {}
    _enumerate_paths(
        graph=graph,
        source=source,
        targets=targets,
        history=history,
        lmax=lmax,
        cycle_free=cycle_free,
    )
    paths = []
    for length in range(lmax):
        for path in history[source][length]:
            paths.append(path)
    return paths


def _enumerate_paths(
    graph: DiGraph,
    source: int,
    targets: List[int],
    lmax: int,
    history,
    cycle_free: bool,
):
    """Enumerate all paths.

    :param graph: graph
    :param source: source node
    :param targets: target nodes
    :param lmax: maximum number of edges to find a path
    :param history:
    :param cycle_free:
    :return:
    """
    # dfs to traverse graph and store history
    # for each visited node
    # containing number of paths already found

    if lmax == 0:
        return

    if graph.nodes[source]['isTarget']:
        history[source] = history.get(source, [[] for _ in range(lmax)])
        if source in targets or not targets:
            history[source][0] = [[source]]
            return
        else:
            return

    paths = [[] for _ in range(lmax)]
    for neighbor in graph.neighbors(source):
        # if not graph.nodes[neighbor]['isTarget'] and lmax == 1:
        #     continue
        if neighbor not in history or neighbor in history and len(history[neighbor]) < lmax - 1:
            _enumerate_paths(
                graph=graph,
                source=neighbor,
                targets=targets,
                lmax=lmax - 1,
                history=history,
                cycle_free=cycle_free,
            )

        # After a node has been visited,
        # history[neighbor] is a dictionary with
        # all paths found for lenghts < lmax.
        # Different keys in dict groups paths
        # by their size (in number of nodes).
        # If there are no paths of length k, k < lmax
        # entry will be None.
        if neighbor in history:
            for length in range(len(history[neighbor][:lmax - 1])):
                for hist_path in history[neighbor][length]:
                    if cycle_free and source in hist_path:
                        continue
                    paths[length + 1].append([source] + hist_path)

    history[source] = paths
