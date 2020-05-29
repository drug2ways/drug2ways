# -*- coding: utf-8 -*-

"""Graph traversal methods."""

from typing import Any, Dict, Iterable, Tuple, List, Union

from networkx import DiGraph

__all__ = [
    'compute_all_paths_multitarget',
    'compute_all_paths_multitarget_dict',
]


def _compute_all_paths_multitarget(graph, source, targets, history, lmax):
    """Compute all paths and store history for each visited node containing number of paths already found.

    :param graph: graph
    :param source: source node
    :param targets: target nodes
    :param history:
    :param lmax:
    :return:
    """
    if lmax == 0:
        return

    paths = [[[0, 0] for _ in targets] for _ in range(lmax)]
    for neighbor in graph.neighbors(source):
        if neighbor in targets or (not targets and graph.nodes[neighbor]['isTarget']):
            target_id = targets.index(neighbor)
            if graph[source][neighbor]['polarity'] == 1:
                paths[0][target_id][0] += 1
            else:
                paths[0][target_id][1] += 1

            continue

        if neighbor not in history or neighbor in history and len(history[neighbor]) < lmax - 1:
            _compute_all_paths_multitarget(graph, neighbor, targets, history, lmax - 1)

        # After a node has been visited,
        # history[neighbor] is a dictionary with
        # all paths found for lenghts < lmax.
        # Different keys in dict groups paths
        # by their size (in number of nodes).

        if neighbor in history:
            for length in range(len(history[neighbor][:lmax - 1])):
                if graph[source][neighbor]['polarity'] == 1:
                    for target_id in range(len(targets)):
                        paths[length + 1][target_id][0] += history[neighbor][length][target_id][0]
                        paths[length + 1][target_id][1] += history[neighbor][length][target_id][1]
                else:
                    for target_id in range(len(targets)):
                        paths[length + 1][target_id][0] += history[neighbor][length][target_id][1]
                        paths[length + 1][target_id][1] += history[neighbor][length][target_id][0]

    history[source] = paths


def _compute_all_paths_multitarget_dict(graph, source, targets, history, lmax):
    """Compute all paths and store history for each visited node containing number of paths already found.

    :param graph: graph
    :param source: source node
    :param targets: target nodes
    :param lmax: lmax
    :param history:
    """
    if lmax == 0:
        return

    paths = [{} for _ in range(lmax)]
    for neighbor in graph.neighbors(source):
        if graph.nodes[neighbor]['isTarget']:
            if neighbor in targets or not targets:
                paths[0][neighbor] = [0, 0]
                if graph[source][neighbor]['polarity'] == 1:
                    paths[0][neighbor][0] += 1
                else:
                    paths[0][neighbor][1] += 1
            continue

        if neighbor not in history or neighbor in history and len(history[neighbor]) < lmax - 1:
            _compute_all_paths_multitarget_dict(graph, neighbor, targets, history, lmax - 1)

        # After a node has been visited,
        # history[neighbor] is a dictionary with
        # all paths found for lenghts < lmax.
        # Different keys in dict groups paths
        # by their size (in number of nodes).

        if neighbor in history:
            for length in range(len(history[neighbor][:lmax - 1])):
                relation_polarity = graph[source][neighbor]['polarity']
                for target in history[neighbor][length].keys():
                    if target not in paths[length + 1]:
                        paths[length + 1][target] = [0, 0]
                    if relation_polarity == 1:
                        paths[length + 1][target][0] += history[neighbor][length][target][0]
                        paths[length + 1][target][1] += history[neighbor][length][target][1]
                    else:
                        paths[length + 1][target][0] += history[neighbor][length][target][1]
                        paths[length + 1][target][1] += history[neighbor][length][target][0]

    history[source] = paths


def get_paths_to_target(graph, source, target, visited_nodes, history, lmax):
    """Return all paths from source to target avoiding nodes in visitied_nodes.

    Paths are returned as a tuple with
        (1) a pair of integers with the total number of paths source to target
        (2) a dictionary of all nodes and how many nodes cross through them
    """
    if lmax == 0:
        if graph[source][target]['polarity'] == 1:
            return [1, 0], {}
        else:
            return [0, 1], {}

    # pre conditions:
    # 1. len(history[source]) is at least lmax
    # 2. target in history[source][lmax]

    valid_neighbors = []
    for node in history[source][lmax][target][1]:
        if graph.has_edge(source, node):
            valid_neighbors.append(node)

    paths = [0, 0]
    node_count = {}
    for neighbor in valid_neighbors:
        # 1. len(history[neighbor]) is at least lmax-1
        if len(history[neighbor]) < lmax - 1:
            continue
        # 2. target in history[neighbor][lmax - 1]
        if target not in history[neighbor][lmax - 1]:
            continue

        # if neighbor is in list of visited_nodes, we're not interested in its paths.
        if neighbor in visited_nodes:
            continue

        # if neighbor has paths involving any of the visited_nodes,
        # recursive call will return w/o them.
        nodes_in_subpath = [node for node in history[neighbor][lmax - 1][target][1].keys()]
        if any(node in visited_nodes + [source] for node in nodes_in_subpath):
            neighbor_paths, neighbor_node_count = get_paths_to_target(
                graph, neighbor, target, visited_nodes + [source], history, lmax - 1)

        else:
            neighbor_paths = history[neighbor][lmax - 1][target][0]
            neighbor_node_count = history[neighbor][lmax - 1][target][1]

        relation_polarity = graph[source][neighbor]['polarity']
        activation_index = 0 if relation_polarity == 1 else 1
        inhibition_index = 1 if relation_polarity == 1 else 0

        for node, count in neighbor_node_count.items():
            if node not in node_count:
                node_count[node] = [0, 0]

            node_count[node][0] += count[activation_index]
            node_count[node][1] += count[inhibition_index]

        if sum(neighbor_paths):
            if neighbor not in node_count:
                node_count[neighbor] = [0, 0]
            node_count[neighbor][0] += neighbor_paths[activation_index]
            node_count[neighbor][1] += neighbor_paths[inhibition_index]
            paths[0] += neighbor_paths[activation_index]
            paths[1] += neighbor_paths[inhibition_index]

    return paths, node_count


def get_paths_through(graph, source, target, crossnode, visited_nodes, history, lmax):
    """Return all paths crossing crossnode.

    Paths are returned as a tuple with
        (1) a pair of integers with the total number of paths source to target.
        (2) a dictionary of all nodes and how many nodes cross through them

    If source == crossnode, call get_paths_to_target to get all paths from here.
    Otherwise, select neighbors that include crossnode in their history of length lmax - 1
    """
    # if source == crossnode, we need to get all paths but source won't be involved anymore,
    # therefore, we have to switch to get_paths_to_target function
    if source == crossnode:
        paths, node_count = get_paths_to_target(
            graph, source, target, visited_nodes + [source], history, lmax - 1
        )

        return paths, node_count

    valid_neighbors = []
    for node in history[source][lmax][target][1]:
        if graph.has_edge(source, node):
            valid_neighbors.append(node)

    paths = [0, 0]
    node_count = {}
    for neighbor in valid_neighbors:
        # 1. len(history[neighbor]) is at least lmax-1
        if len(history[neighbor]) < lmax - 1:
            continue
        # 2. target in history[neighbor][lmax - 1]
        if target not in history[neighbor][lmax - 1]:
            continue
        # 3. neighbor not in visited_nodes
        if neighbor in visited_nodes:
            continue

        elif neighbor == crossnode:
            neighbor_paths, neighbor_node_count = get_paths_to_target(
                graph, neighbor, target, visited_nodes + [source], history, lmax - 1)

        # crossnode must be in intermediate paths for neighbor to be considered
        elif crossnode not in history[neighbor][lmax - 1][target][1]:
            continue
        else:
            neighbor_paths, neighbor_node_count = get_paths_through(
                graph, neighbor, target, crossnode, visited_nodes + [source], history, lmax - 1)

        relation_polarity = graph[source][neighbor]['polarity']
        activation_index = 0 if relation_polarity == 1 else 1
        inhibition_index = 1 if relation_polarity == 1 else 0

        for node, count in neighbor_node_count.items():
            if node not in node_count:
                node_count[node] = [0, 0]

            node_count[node][0] += count[activation_index]
            node_count[node][1] += count[inhibition_index]

        if sum(neighbor_paths):
            if neighbor not in node_count:
                node_count[neighbor] = [0, 0]
            node_count[neighbor][0] += neighbor_paths[activation_index]
            node_count[neighbor][1] += neighbor_paths[inhibition_index]
            paths[0] += neighbor_paths[activation_index]
            paths[1] += neighbor_paths[inhibition_index]

    return paths, node_count


def _compute_all_simple_paths_multitarget_dict(graph, source, targets, lmax, history, cycle_history={}):
    """Compute all simple paths and store history for each visited node containing number of paths already found.

    :param graph: graph
    :param source: source node
    :param targets: target nodes
    :param lmax: maximum number of edges to find a path
    :param history: maximum number of edges to find a path
    :param cycle_history:
    :return:
    """
    if lmax == 0:
        return

    _history_source = [{} for _ in range(lmax)]
    neighbor_nodes = [neighbor for neighbor in graph.neighbors(source)]

    for neighbor in neighbor_nodes:
        if graph.nodes[neighbor]['isTarget']:
            if neighbor in targets or not targets:
                _history_source[0][neighbor] = [[0, 0], {}]
                if graph[source][neighbor]['polarity'] == 1:
                    _history_source[0][neighbor][0][0] = 1
                else:
                    _history_source[0][neighbor][0][1] = 1

            continue

        if neighbor not in history or (neighbor in history and len(history[neighbor]) < lmax - 1):
            _compute_all_simple_paths_multitarget_dict(
                graph=graph,
                source=neighbor,
                targets=targets,
                lmax=lmax - 1,
                history=history,
                cycle_history=cycle_history,
            )

        if neighbor in history:

            for length in range(len(history[neighbor][:lmax - 1])):

                relation_polarity = graph[source][neighbor]['polarity']
                for target in history[neighbor][length].keys():

                    if target not in _history_source[length + 1]:
                        _history_source[length + 1][target] = [[0, 0], {}]

                    activation_index = 0 if relation_polarity == 1 else 1
                    inhibition_index = 1 if relation_polarity == 1 else 0

                    # Get paths to target
                    paths_to_target = [history[neighbor][length][target][0][0], history[neighbor][length][target][0][1]]

                    # If there're no cycles, intermediate_nodes will not be modified
                    # and there's no need for a costly deepcopy operation
                    intermediate_nodes = history[neighbor][length][target][1]

                    # Get paths starting from neighbor that have source as intermediate node
                    paths_in_cycle = intermediate_nodes.get(source, [0, 0])

                    if paths_in_cycle[0] or paths_in_cycle[1]:
                        # Cycle detected.
                        # Find all paths / nodes involved in cycle so they can be removed
                        # from neighbor's list before being added to source's.
                        nodes_in_cycles = {}
                        number_of_paths_in_cycles = []
                        if (
                            source in cycle_history[0]
                            and neighbor in cycle_history[0][source]
                            and target in cycle_history[0][source][neighbor]
                            and length in cycle_history[0][source][neighbor][target]
                        ):
                            number_of_paths_in_cycles, nodes_in_cycles = cycle_history[0][source][neighbor][target][
                                length]

                        else:
                            number_of_paths_in_cycles, nodes_in_cycles = get_paths_through(
                                graph, neighbor, target, source, [], history, length,
                            )

                            if source not in cycle_history[0]:
                                cycle_history[0][source] = {}
                            if neighbor not in cycle_history[0][source]:
                                cycle_history[0][source][neighbor] = {}
                            if target not in cycle_history[0][source][neighbor]:
                                cycle_history[0][source][neighbor][target] = {}
                            if length not in cycle_history[0][source][neighbor][target]:
                                cycle_history[0][source][neighbor][target][length] = [number_of_paths_in_cycles,
                                                                                      nodes_in_cycles]

                        if target not in cycle_history[1]:
                            cycle_history[1][target] = 0
                        cycle_history[1][target] += 1

                        # Manually copy intermediate_nodes to avoid calling deepcopy
                        intermediate_nodes = {
                            node: [paths[0], paths[1]]
                            for node, paths in history[neighbor][length][target][1].items()
                        }

                        for node, paths in nodes_in_cycles.items():
                            intermediate_nodes[node][0] -= paths[0]
                            intermediate_nodes[node][1] -= paths[1]

                            if intermediate_nodes[node] == [0, 0]:
                                intermediate_nodes.pop(node, None)

                        paths_to_target[0] -= number_of_paths_in_cycles[0]
                        paths_to_target[1] -= number_of_paths_in_cycles[1]

                    # Copy map with path count per node
                    for node, paths in intermediate_nodes.items():
                        if node not in _history_source[length + 1][target][1]:
                            _history_source[length + 1][target][1][node] = [0, 0]

                        _history_source[length + 1][target][1][node][0] += paths[activation_index]
                        _history_source[length + 1][target][1][node][1] += paths[inhibition_index]

                    # Add neighbor to the intermediate node's list
                    if paths_to_target[0] or paths_to_target[1]:
                        _history_source[length + 1][target][1][neighbor] = _history_source[length + 1][target][1].get(
                            neighbor,
                            [0, 0],
                        )
                        _history_source[length + 1][target][1][neighbor][0] += paths_to_target[activation_index]
                        _history_source[length + 1][target][1][neighbor][1] += paths_to_target[inhibition_index]

                    # Update count of paths to target from source
                    _history_source[length + 1][target][0][0] += paths_to_target[activation_index]
                    _history_source[length + 1][target][0][1] += paths_to_target[inhibition_index]

    history[source] = _history_source


def compute_all_paths_multitarget(
    graph: DiGraph,
    source: Iterable[Any],
    targets,
    lmax: int,
    previous_history: Dict,
) -> Tuple[int, List[List[Union[float, int]]]]:
    """Compute all paths to all the targets, separately.

    :param graph: graph
    :param source: source nodes
    :param targets: target nodes
    :param lmax: lmax
    :param previous_history: previous history of visited nodes
    :return:
    """
    _compute_all_paths_multitarget(
        graph=graph,
        source=source,
        targets=targets,
        history=previous_history,
        lmax=lmax,
    )

    # Initialize [relative_activations, relative_inhibitions, number_of_paths_to_target]
    count = [[0.0, 0.0, 0] for _ in targets]
    total_number_of_paths = 0
    # Loop up in the history the results in each of the paths
    for t in range(len(targets)):
        for length in range(lmax):
            source_to_target = previous_history[source][length][t]
            count[t][0] += source_to_target[0]
            count[t][1] += source_to_target[1]
            count[t][2] += source_to_target[0] + source_to_target[1]

        total_number_of_paths = _count_paths(count, t, total_number_of_paths)

    return total_number_of_paths, count


def compute_all_paths_multitarget_dict(
    graph: DiGraph,
    source,
    targets,
    lmax: int,
    previous_history: Dict,
    cycle_history: Union[Dict, Dict],
    simple_paths: bool = False,
) -> Tuple[int, List[List[Union[float, int]]]]:
    """Compute all paths to all the targets, separately. Uses dict to store target path count instead of array.

    :param graph: graph
    :param source: source node
    :param targets: target nodes
    :param lmax: lmax
    :param previous_history: history of visited nodes
    :param cycle_history: history of cycles
    :param simple_paths: simple paths mode
    """
    if simple_paths:
        _compute_all_simple_paths_multitarget_dict(
            graph=graph,
            source=source,
            targets=targets,
            history=previous_history,
            lmax=lmax,
            cycle_history=cycle_history,
        )
    else:
        _compute_all_paths_multitarget_dict(
            graph=graph,
            source=source,
            targets=targets,
            history=previous_history,
            lmax=lmax,
        )

    # Initialize [relative_activations, relative_inhibitions, number_of_paths_to_target]
    count = [[0.0, 0.0, 0] for _ in targets]
    total_number_of_paths = 0

    # Loop up in the history the results in each of the paths
    for t in range(len(targets)):
        target_id = targets[t]
        for length in range(lmax):
            if target_id not in previous_history[source][length]:
                continue

            source_to_target = previous_history[source][length][target_id]
            if simple_paths:
                source_to_target = source_to_target[0]
            count[t][0] += source_to_target[0]
            count[t][1] += source_to_target[1]
            count[t][2] += source_to_target[0] + source_to_target[1]

        total_number_of_paths = _count_paths(count, t, total_number_of_paths)

    return total_number_of_paths, count


def _count_paths(count, t, total_number_of_paths):
    """Count ratio and number of paths."""
    total_number_of_paths += count[t][2]
    # Calculate ratio of positive and negative paths
    count[t][0] = 0 if not count[t][2] else round(count[t][0] / count[t][2], 2)
    count[t][1] = 0 if not count[t][2] else round(count[t][1] / count[t][2], 2)
    return total_number_of_paths
