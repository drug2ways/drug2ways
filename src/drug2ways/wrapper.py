# -*- coding: utf-8 -*-

"""Wrapper methods to calculate paths for a set of source(s) and target(s)."""

import itertools as itt
import json
import logging
import os
import time
from collections import defaultdict
from itertools import combinations
from typing import Any, Iterable, List, Mapping, Optional

from networkx import DiGraph
from tqdm import tqdm

from .alternative_graph_traversal import enumerate_paths
from .graph_processing import generate_reduced_graph
from .graph_traversal import compute_all_paths_multitarget_dict
from .pathway import analyze_paths

__all__ = [
    'wrapper_explore',
    'wrapper_optimize',
    'wrapper_combine',
    'wrapper_pathway_enrichment',
]

logger = logging.getLogger(__name__)


def _check_generic_input(graph: DiGraph, source_nodes: List[Any], target_nodes: List[Any]) -> None:
    """Check source and targets are in the graph.

    :param graph: directed graph
    :param source_nodes: iterable with sources nodes (usually drugs)
    :param target_nodes: iterable with target nodes (usually diseases)
    """
    # Ensure both source node and target nodes are in graph
    if not all(node in graph for node in source_nodes):
        nodes_not_in_graph = {node for node in source_nodes if node not in graph}
        raise ValueError(
            f'The following source nodes are not in the graph ({len(nodes_not_in_graph)}/{len(source_nodes)}):',
            f'({nodes_not_in_graph})',
        )

    if not all(node in graph for node in target_nodes):
        nodes_not_in_graph = {node for node in target_nodes if node not in graph}
        raise ValueError(
            f'The following target nodes are not in the graph ({len(nodes_not_in_graph)}/{len(target_nodes)}):',
            f'({nodes_not_in_graph})',
        )


def _check_optimize_input(
    graph: DiGraph,
    source_nodes: List[Any],
    activate_targets: List[Any],
    inhibit_targets: List[Any],
) -> None:
    """Check source and targets are in the graph.

    :param graph: directed graph
    :param source_nodes: iterable with sources nodes (usually drugs)
    :param activate_targets: iterable with nodes to be activated
    :param inhibit_targets: iterable with nodes to be inhibited
    """
    # Ensure both source node and activation/inhibition nodes are in graph
    if not all(node in graph for node in source_nodes):
        nodes_not_in_graph = {node for node in source_nodes if node not in graph}
        raise ValueError(f'The following source nodes are not in the graph: ({nodes_not_in_graph})')

    if not all(node in graph for node in activate_targets):
        raise ValueError(f'At least one of the target nodes to activate ({activate_targets}) is not in the graph')

    if not all(node in graph for node in inhibit_targets):
        raise ValueError(f'At least one of the target nodes to inhibit ({inhibit_targets}) is not in the graph')


def wrapper_explore(graph: DiGraph, source_nodes: List[Any], target_nodes: List[Any], lmax: int, simple_paths: bool):
    """Calculate the effect of sources nodes over multiple target nodes on a directed graph.

    :param graph: directed graph
    :param source_nodes: iterable with sources nodes (usually drugs)
    :param target_nodes: iterable with target nodes (usually diseases)
    :param lmax: maximum length of the path allowed
    :param simple_paths: if true, only simple paths are calculated
    :return: effect of the source node on each target node
    """
    _check_generic_input(graph, source_nodes, target_nodes)

    results_by_source: List = []
    time_cache = defaultdict(dict)

    # Store history of the visited path for the given node
    previous_history = {}
    # Cycle History:
    # [0] Dict to store pre-calculated cycles
    # [1] Dict to store number of cycles from source to target
    cycle_history = [{}, {}]

    # Path count to all targets, by source
    count_by_source = {}
    cycles_by_source = {}

    # Get the reduced version of the graph and the node2id mapping
    reduced_graph, node2id = generate_reduced_graph(graph, target_nodes)

    # Initialize MPI environment and variables, if found
    number_of_processes = 1
    process_id = 0
    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        process_id = comm.Get_rank()
        number_of_processes = comm.Get_size()
    except ImportError:
        pass

    _target_nodes = [node2id[target_node] for target_node in target_nodes]
    _source_nodes = [source_node for source_node in source_nodes]
    if process_id > 0 or number_of_processes == 1:

        source_index = 0
        if process_id > 0:
            source_node = comm.recv(source=0, tag=process_id)
        else:
            source_node = _source_nodes[source_index]

        work_done = 0
        while source_node != -1:
            # Get node identifiers
            _source_node = node2id[source_node]

            exe_t_0 = time.time()

            # Calculate all paths between source and target
            _, count = compute_all_paths_multitarget_dict(
                graph=reduced_graph,
                source=_source_node,
                targets=_target_nodes,
                lmax=lmax,
                previous_history=previous_history,
                cycle_history=cycle_history,
                simple_paths=simple_paths
            )

            # Cache time needed
            exe_t_f = time.time()
            time_cache[source_node] = exe_t_f - exe_t_0

            count_by_source[source_node] = count
            cycles_by_source[source_node] = cycle_history[1]
            cycle_history[1] = {}
            if process_id > 0:
                comm.send(process_id, dest=0, tag=0)
                source_node = comm.recv(source=0, tag=process_id)
            else:
                source_index += 1
                if source_index >= len(_source_nodes):
                    source_node = -1
                else:
                    source_node = _source_nodes[source_index]
            work_done += 1

    else:
        # process_id = 0 and number_of_processes > 1
        free_workers = [worker for worker in range(1, number_of_processes)]
        status = MPI.Status()  # Get MPI status object
        for source_index, source_node in tqdm(enumerate(source_nodes), total=len(source_nodes)):
            if free_workers:
                worker = free_workers.pop(0)
                comm.send(source_node, dest=worker, tag=worker)

            else:
                # Wait until a worker finishes
                worker = comm.recv(
                    source=MPI.ANY_SOURCE, tag=0, status=status)
                req = comm.isend(source_node, dest=worker, tag=worker)

        # Wait until all workers have finished their work
        while len(free_workers) < number_of_processes - 1:
            worker = comm.recv(
                source=MPI.ANY_SOURCE, tag=0, status=status)
            free_workers.append(worker)

        # After all source nodes are processed, notify workers that there's no more work.
        for worker in range(1, number_of_processes):
            code = -1
            comm.send(code, dest=worker, tag=worker)

    # Master (process_id == 0) receives partial results from other processes
    if process_id == 0 and number_of_processes > 1:
        # print(f'Results from master: {len(count_by_source)}')
        print(f'Waiting to receive partial results from {number_of_processes - 1} other processes.')
        status = MPI.Status()  # Get MPI status object

        for i in range(number_of_processes - 1):
            partial_results = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
            source = status.Get_source()
            print(f'Received {len(partial_results)} from process {source}.')
            count_by_source.update(partial_results[0])
            cycles_by_source.update(partial_results[1])

        # Order results by source id. This is only to make validation easier
        # for source_node in source_nodes:
        #     results_by_source.append(recv_results[source_node])

    elif number_of_processes > 1:
        print(f'[{process_id}] Sending results to master.')
        results = [count_by_source, cycles_by_source]
        comm.send(results, dest=0, tag=process_id)

    if process_id == 0:
        for target_index, target_node in enumerate(target_nodes):
            # target_id = node2id[target_node]
            for source_node in source_nodes:
                if source_node not in count_by_source:
                    continue
                relative_activation, relative_inhibition, paths_to_target = count_by_source[source_node][target_index]

                # number_of_cycles = 0
                # if source_node in cycles_by_source:
                #     if target_id in cycles_by_source[source_node]:
                #         number_of_cycles = cycles_by_source[source_node][target_id]

                if not paths_to_target:
                    continue

                results_by_source.append(
                    dict(
                        source=source_node,
                        target=target_node,
                        relative_activation=relative_activation,
                        relative_inhibition=relative_inhibition,
                        number_of_paths=paths_to_target,
                        # number_of_cycles=number_of_cycles,
                    )
                )
    if not results_by_source:
        logger.warning('There are no paths between the sources and any targets')

    return results_by_source, time_cache


def wrapper_optimize(
    graph: DiGraph,
    source_nodes: List[Any],
    activate_targets: List[Any],
    inhibit_targets: List[Any],
    activation_threshold: Optional[float],
    lmax: int,
    simple_paths: bool,
):
    """Optimize the effect of sources nodes over multiple target nodes on a directed graph.

    :param graph: directed graph
    :param source_nodes: iterable with sources nodes (usually drugs)
    :param activate_targets: iterable with target nodes to activate (usually diseases)
    :param inhibit_targets: iterable with target nodes to inhibit (usually diseases)
    :param activation_threshold: threshold to consider a target as activated (>=) or inhibited (<).
    :param lmax: maximum length of the path allowed
    :param simple_paths: if true, only simple paths are calculated
    :return: set of sources optimizing the specified targets and their effects on each target node.
    """
    _check_optimize_input(graph, source_nodes, activate_targets, inhibit_targets)

    results_by_source: List = []

    target_nodes = activate_targets + inhibit_targets

    # Get the reduced version of the graph and the node2id mapping
    reduced_graph, node2id = generate_reduced_graph(graph, target_nodes)

    # Store history of the visited path for the given node
    previous_history = {}
    # Cycle History:
    # [0] Dict to store pre-calculated cycles
    # [1] Dict to store number of cycles from source to target
    cycle_history = [{}, {}]

    # TODO: Use same parallelization mechanism as in explore

    _target_nodes = [node2id[target_node] for target_node in target_nodes]
    for source_node in source_nodes:

        # Get node identifiers
        _source_node = node2id[source_node]

        # Calculate all paths between source and target
        number_of_paths, count = compute_all_paths_multitarget_dict(
            graph=reduced_graph,
            source=_source_node,
            targets=_target_nodes,
            lmax=lmax,
            previous_history=previous_history,
            cycle_history=cycle_history,
            simple_paths=simple_paths,
        )

        # Skip if there are no paths
        if not number_of_paths:
            continue

        # Check if targets are optimized
        optimized = True

        # Filter results by activation threshold
        if activation_threshold:
            # TODO: What if there are no paths to a target?
            for target in range(len(activate_targets)):
                if count[target][0] < activation_threshold:
                    optimized = False
                    break

            if optimized:
                for target in range(len(inhibit_targets)):
                    if count[target + len(activate_targets)][0] >= activation_threshold:
                        optimized = False
                        break

        if activation_threshold and not optimized:
            continue

        # Store the results for this given target if it is interesting
        results_by_source.append(
            dict(
                source=source_node,
                activate_targets=activate_targets,
                inhibit_targets=inhibit_targets,
                relative_activations=[target[0] for target in count],
                relative_inhibitions=[target[1] for target in count],
                number_of_paths=[target[2] for target in count],
                total_number_of_paths=number_of_paths,
            )
        )

    if not results_by_source:
        logger.warning('There are no paths between the sources and any targets')

    return results_by_source


def wrapper_combine(
    graph: DiGraph,
    source_nodes: List[Any],
    activate_targets: List[Any],
    inhibit_targets: List[Any],
    activation_threshold: Optional[float],
    combination_length: int,
    lmax: int,
    simple_paths: bool,
):
    """Calculate the effect of a combination of source nodes over multiple target nodes on a directed graph.

    :param graph: directed graph
    :param source_nodes: iterable with sources nodes (usually drugs)
    :param activate_targets: iterable with target nodes to activate (usually diseases)
    :param inhibit_targets: iterable with target nodes to inhibit (usually diseases)
    :param activation_threshold: threshold to consider a target as activated (>=) or inhibited (<).
    :param combination_length: number of source nodes of each combination.
    :param lmax: maximum length of the path allowed
    :param simple_paths: if true, only simple paths are calculated
    :return: set of sources optimizing the specified targets and their effects on each target node.
    """
    _check_optimize_input(graph, source_nodes, activate_targets, inhibit_targets)

    results_by_combination: List = []

    target_nodes = activate_targets + inhibit_targets

    # Get the reduced version of the graph and the node2id mapping
    reduced_graph, node2id = generate_reduced_graph(graph, target_nodes)

    # Store history of the visited path for the given node
    previous_history = {}
    # Cycle History:
    # [0] Dict to store pre-calculated cycles
    # [1] Dict to store number of cycles from source to target
    cycle_history = [{}, {}]
    # TODO: Check best way to parallelize this function.

    _target_nodes = [node2id[target_node] for target_node in target_nodes]

    for combination in combinations(source_nodes, combination_length):
        combination_count = [
            [0, 0, 0]
            for _ in target_nodes
        ]
        number_of_paths_in_combination = 0

        for source_node in combination:

            # Get node identifiers
            _source_node = node2id[source_node]

            # Calculate all paths between source and target
            number_of_paths, count = compute_all_paths_multitarget_dict(
                graph=reduced_graph,
                source=_source_node,
                targets=_target_nodes,
                lmax=lmax,
                cycle_history=cycle_history,
                previous_history=previous_history,
                simple_paths=simple_paths,
            )

            # TODO: Skip combination if there are no paths from one source?
            if not number_of_paths:
                break

            number_of_paths_in_combination += number_of_paths
            for target_index, count_to_target in enumerate(count):
                combination_count[target_index][0] += count_to_target[0] * count_to_target[2]
                combination_count[target_index][1] += count_to_target[1] * count_to_target[2]
                combination_count[target_index][2] += count_to_target[2]

        if not number_of_paths_in_combination:
            continue

        # Calculate relative activation/inhibition for the combination
        for count_to_target in combination_count:
            # Skip if there are no paths to target
            if not count_to_target[2]:
                continue

            count_to_target[0] = count_to_target[0] / count_to_target[2]
            count_to_target[1] = count_to_target[1] / count_to_target[2]

        # Check if targets are optimized in the combination
        optimized = True

        if activation_threshold:
            # TODO: What if there are no paths to a target?
            for target_index in range(len(activate_targets)):
                if combination_count[target_index][0] < activation_threshold:
                    optimized = False
                    break

            for _target, combination in zip(activate_targets, combination_count):
                if combination[0] < activation_threshold:
                    optimized = False
                    break

            if optimized:
                for target_index in range(len(inhibit_targets)):
                    if combination_count[target_index + len(activate_targets)][0] >= activation_threshold:
                        optimized = False
                        break

        # Skip this combination if it is not optimizing the desired results
        if activation_threshold and not optimized:
            continue

        # Store the results for this given target if it is interesting
        results_by_combination.append(
            dict(
                combination=combination,
                activate_targets=activate_targets,
                inhibit_targets=inhibit_targets,
                relative_activations=[target[0] for target in combination_count],
                relative_inhibitions=[target[1] for target in combination_count],
                number_of_paths_to_targets=[target[2] for target in combination_count],
                total_number_of_paths=number_of_paths_in_combination,
            )
        )

    if not results_by_combination:
        logger.warning('There are no paths between the sources and any targets')

    return results_by_combination


def wrapper_pathway_enrichment(
    graph: DiGraph,
    source_nodes: List[Any],
    target_nodes: List[Any],
    lmax: int,
    simple_paths: bool,
    export: bool = True,
    output: Optional[str] = None,
    genesets: Mapping[str, Iterable[str]] = None,
):
    """Conduct pathway enrichment on the paths.

    :param graph: directed graph
    :param source_nodes: iterable with sources nodes (usually drugs)
    :param target_nodes: iterable with target nodes (usually diseases)
    :param lmax: maximum length of the path allowed
    :param simple_paths: if true, only simple paths are calculated
    :param output: output directory
    :param export: export df and results
    :param genesets: genesets
    :return: results of the pathway enrichment.
    """
    _check_generic_input(graph, source_nodes, target_nodes)

    # Get the reduced version of the graph and the node2id mapping
    reduced_graph, node2id = generate_reduced_graph(graph, target_nodes)

    id2node = {
        v: k
        for k, v in node2id.items()
    }

    _target_nodes = [node2id[target_node] for target_node in target_nodes]

    results = {}

    for source_node, target_node in itt.product(source_nodes, target_nodes):

        target_id = node2id[target_node]

        # Calculate all paths between source and target
        paths = enumerate_paths(
            graph=reduced_graph,
            source=node2id[source_node],
            targets=[target_id],
            lmax=lmax,
            cycle_free=simple_paths,
        )

        # Skip if there are no paths
        if not paths:
            logger.warning(f'No paths between {source_node} and {target_node}')
            continue

        # Get summarized results for export
        summarized_results, paths_summary, enrichment_results = analyze_paths(
            reduced_graph=reduced_graph,
            paths=paths,
            id2node=id2node,
            genesets=genesets,
        )

        results[(source_node, target_node)] = (summarized_results, paths_summary, enrichment_results)

        if export:
            # export_results
            logger.info(f'Exporting results for pair {source_node} - {target_node}')
            # TODO Fix the lmax later since now it is hacked
            summarized_results.to_csv(
                os.path.join(output, f'overview-{lmax - 1}-{source_node}_{target_node}.tsv'),
                sep='\t',
                index=None,
            )

            with open(
                os.path.join(output, f'paths-detailed{lmax - 1}-{source_node}_{target_node}.tsv'), 'w'
            ) as file:
                json.dump(paths_summary, file, indent=2)

    return results
