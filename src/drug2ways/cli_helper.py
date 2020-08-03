# -*- coding: utf-8 -*-

"""Helper for cli."""

import json
import logging
import os
import time
from typing import List, Optional, Union

import click
from networkx import DiGraph

from drug2ways.bel_helper import get_candidate_drugs, get_candidate_targets
from drug2ways.constants import BEL_FORMATS, EMOJI
from drug2ways.graph_reader import load_graph, read_nodes, read_nodes_to_optimize
from drug2ways.pathway import get_genesets
from drug2ways.wrapper import wrapper_explore, wrapper_optimize, wrapper_combine, wrapper_pathway_enrichment

logger = logging.getLogger(__name__)


def _handle_lmax_parameter(lmax) -> Optional[List[int]]:
    """Prepare iterable for lmax."""
    if isinstance(lmax, int):
        return [lmax]
    elif isinstance(lmax, List):
        return lmax
    else:
        raise ValueError('Unknown type for lmax')


def check_graph_input(sources: str, targets: str, fmt: str, drug_search_bel: bool) -> None:
    """Validate input from cli."""
    if fmt not in BEL_FORMATS and drug_search_bel:
        raise ValueError(
            f'{EMOJI} The Drug-Search functionality only works for BEL graphs. '
            f'Please try with a BEL graph or do not use this functionality {EMOJI}'
        )

    elif drug_search_bel and fmt in BEL_FORMATS:
        pass

    elif not sources:
        raise ValueError(
            f'{EMOJI} You have not provided sources nodes or selected the Drug Search option. {EMOJI}'
        )

    elif not targets:
        raise ValueError(
            f'{EMOJI} You have not provided target nodes or selected the Drug Search option {EMOJI}'
        )


def _prepare_graph(graph: str, fmt: str, sources: str, targets: str):
    """Parallelization of the calculations to improve efficiency."""
    # Initialize MPI environment and variables, if found.
    number_of_processes = 1
    process_id = 0
    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        process_id = comm.Get_rank()
        number_of_processes = comm.Get_size()
        click.echo("Parallelization on")

    except ImportError:
        pass
    # Ensure file is valid
    check_graph_input(sources, targets, fmt, False)
    if process_id == 0:
        click.secho(f'{EMOJI} Loading the graph from {graph} {EMOJI}')
    # Load graph
    directed_graph: DiGraph = load_graph(graph, fmt)
    # Get source nodes from file
    source_nodes = read_nodes(sources)
    # Get target nodes from file
    activate_targets, inhibit_targets = read_nodes_to_optimize(targets)
    if process_id == 0:
        number_of_targets = len(activate_targets + inhibit_targets)
        click.echo(
            f"{EMOJI} A total of {len(source_nodes)} sources and {number_of_targets} targets will be used. {EMOJI}"
        )

        if number_of_processes > 1:
            click.echo(f"{EMOJI} Distributing work among {number_of_processes} processes. {EMOJI}")

    return directed_graph, source_nodes, activate_targets, inhibit_targets, process_id


def _explore_helper(
    graph: str,
    fmt: str,
    sources: str,
    targets: str,
    lmax: Union[int, List[int]],
    simple_paths: bool,
    log: bool,
    export_time: bool,
    output: str,
    name: str,
    drug_search_bel: bool,
) -> None:
    """Wrap explore command in cli.

    :param graph: path to the graph
    :param fmt: graph format
    :param sources: path to the source nodes
    :param targets: path to the target nodes
    :param lmax: max length of path
    :param simple_paths: use simple paths or cycles
    :param log: debug mode
    :param export_time: export time for each pair
    :param output: output directory
    :param name: name of the graph for output purposes
    :param drug_search_bel: search drugs automatically in BEL
    """
    """Parallelization of the calculations to improve efficiency"""
    # Initialize MPI environment and variables, if found.
    number_of_processes = 1
    process_id = 0
    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        process_id = comm.Get_rank()
        number_of_processes = comm.Get_size()
        click.echo("Parallelization on")
    except ImportError:
        pass

    # Setup logging.
    # MPI ranks other than master will have it set to WARNING
    _setup_logging(log, process_id)

    # Ensure file is valid
    check_graph_input(sources, targets, fmt, drug_search_bel)
    if process_id == 0:
        click.secho(f'{EMOJI} Loading the graph from {graph} {EMOJI}')
    # Load graph
    directed_graph: DiGraph = load_graph(graph, fmt)
    # Find drugs (sources) and pathologies (targets) if selected and grpah is in BEL fmt
    if drug_search_bel:
        # Get source nodes (drugs)
        source_nodes = get_candidate_drugs(directed_graph)
        # Get target nodes (pathologies)
        targets_nodes = get_candidate_targets(directed_graph)
    else:
        # Get source nodes from file
        source_nodes = read_nodes(sources)
        # Get target nodes from file
        targets_nodes = read_nodes(targets)
    if process_id == 0:
        click.echo(
            f"{EMOJI} A total of {len(source_nodes)} sources and {len(targets_nodes)} targets will be used. {EMOJI}"
        )
        if number_of_processes > 1:
            click.echo(f"{EMOJI} Distributing work among {number_of_processes} processes. {EMOJI}")

    path_mode = "simple-paths" if simple_paths else "cycles"

    """Export results for each lmax"""
    for lmax in _handle_lmax_parameter(lmax):
        click.secho(
            f'{EMOJI} Calculating paths with lmax ({lmax}) on {path_mode} mode. This might take a while... {EMOJI}',
        )

        # Warn user if lmax larger than 12
        if lmax > 12:
            logger.warning(
                f"Note that the selected Lmax '{lmax}' might converge results if your graph is not large enough"
            )

        # Track the time it takes to run
        exe_t_0 = time.time()
        # Call main function
        results, time_cache = wrapper_explore(
            graph=directed_graph,
            source_nodes=source_nodes,
            target_nodes=targets_nodes,
            lmax=lmax,
            simple_paths=simple_paths,
        )
        # Finished time
        exe_t_f = time.time()
        running_time = exe_t_f - exe_t_0
        if process_id == 0:
            click.secho(f'{EMOJI} Finished in {running_time} seconds {EMOJI}')

            # Export results
            with open(os.path.join(output, f'{name}all_against_all_lmax_{lmax}.json'), 'w') as f:
                json.dump(results, f, indent=2)

            # Export time
            if export_time:
                with open(os.path.join(output, f'{name}time_cache_{lmax}.json'), 'w') as f:
                    json.dump(time_cache, f, indent=2)

            click.secho(f'{EMOJI} Results exported to {output} {EMOJI}')


def _optimize_helper(
    graph: str,
    fmt: str,
    sources: str,
    targets: str,
    lmax: int,
    simple_paths: bool,
    log: bool,
    output: str,
    activation_threshold: float,
) -> None:
    """Wrap optimize command in cli.

    :param graph: graph
    :param fmt: format
    :param sources: source nodes
    :param targets: target nodes
    :param lmax: max length of path
    :param simple_paths: use simple paths or cycles
    :param log: debug mode
    :param output: output directory
    :param activation_threshold: filtering
    """
    _setup_logging(log, lmax)

    directed_graph, source_nodes, activate_targets, inhibit_targets, process_id = _prepare_graph(
        graph=graph,
        fmt=fmt,
        sources=sources,
        targets=targets,
    )

    path_mode = "simple-paths" if simple_paths else "cycles"
    click.secho(
        f'{EMOJI} Calculating paths with lmax ({lmax}) on {path_mode} mode. This might take a while... {EMOJI}',
    )

    # Track the time it takes to run
    exe_t_0 = time.time()
    # Call main function
    results = wrapper_optimize(
        graph=directed_graph,
        source_nodes=source_nodes,
        activate_targets=activate_targets,
        inhibit_targets=inhibit_targets,
        activation_threshold=activation_threshold,
        lmax=lmax,
        simple_paths=simple_paths,
    )
    # Finished time
    exe_t_f = time.time()
    running_time = exe_t_f - exe_t_0
    if process_id == 0:
        click.secho(f'{EMOJI} Finished in {running_time} seconds {EMOJI}')

        # Export results
        with open(os.path.join(output, f'optimize_lmax_{lmax}.json'), 'w') as f:
            json.dump(results, f, indent=2)

        click.secho(f'{EMOJI} Results exported to {output} {EMOJI}')


def _combine_helper(
    graph: str,
    fmt: str,
    sources: str,
    targets: str,
    lmax: int,
    simple_paths: bool,
    log: bool,
    output: str,
    combination_length: int,
    activation_threshold: float,
) -> None:
    """Wrap combination command in cli.

    :param graph: graph
    :param fmt: format
    :param sources: source nodes
    :param targets: target nodes
    :param lmax: max length of path
    :param simple_paths: use simple paths or cycles
    :param log: debug mode
    :param output: output directory
    :param combination_length: number of combinations
    :param activation_threshold: filtering
    """
    _setup_logging(log, lmax)

    directed_graph, source_nodes, activate_targets, inhibit_targets, process_id = _prepare_graph(
        graph=graph,
        fmt=fmt,
        sources=sources,
        targets=targets,
    )

    path_mode = "simple-paths" if simple_paths else "cycles"

    """Export results for each lmax"""
    for lmax in _handle_lmax_parameter(lmax):
        click.secho(
            f'{EMOJI} Calculating paths with lmax ({lmax}) on {path_mode} mode. This might take a while... {EMOJI}',
        )

        # Track the time it takes to run
        exe_t_0 = time.time()
        # Call main function
        results = wrapper_combine(
            graph=directed_graph,
            source_nodes=source_nodes,
            activate_targets=activate_targets,
            inhibit_targets=inhibit_targets,
            activation_threshold=activation_threshold,
            combination_length=combination_length,
            lmax=lmax,
            simple_paths=simple_paths,
        )
        # Finished time
        exe_t_f = time.time()
        running_time = exe_t_f - exe_t_0
        if process_id == 0:
            click.secho(f'{EMOJI} Finished in {running_time} seconds {EMOJI}')

            # Export results
            with open(os.path.join(output, f'combine_{combination_length}_drugs_lmax_{lmax}.json'), 'w') as f:
                json.dump(results, f, indent=2)

            click.secho(f'{EMOJI} Results exported to {output} {EMOJI}')


def _pathway_enrichment_helper(
    graph: str,
    fmt: str,
    sources: str,
    targets: str,
    lmax: int,
    simple_paths: bool,
    log: bool,
    output: str,
) -> None:
    """Wrap optimize command in cli.

    :param graph: graph
    :param fmt: format
    :param sources: source nodes
    :param targets: target nodes
    :param lmax: max length of path
    :param simple_paths: use simple paths or cycles
    :param log: debug mode
    :param output: output directory
    """
    _setup_logging(log)

    logger.debug('Getting genesets')
    kegg, reactome, wikipathways = get_genesets()

    # Ensure file is valid
    check_graph_input(sources, targets, fmt, False)

    # Load graph
    directed_graph: DiGraph = load_graph(graph, fmt)

    # Get source nodes from file
    source_nodes = read_nodes(sources)
    # Get target nodes from file
    targets_nodes = read_nodes(targets)

    click.echo(
        f"{EMOJI} A total of {len(source_nodes)} sources and {len(targets_nodes)} targets will be used. {EMOJI}"
    )

    path_mode = "simple-paths" if simple_paths else "cycles"

    """Export results for each lmax"""
    for lmax in _handle_lmax_parameter(lmax):
        click.secho(
            f'{EMOJI} Pathway analysis with lmax ({lmax}) on {path_mode} mode. This might take a while... {EMOJI}',
        )

        # Track the time it takes to run
        exe_t_0 = time.time()
        # Call main function
        _ = wrapper_pathway_enrichment(
            graph=directed_graph,
            source_nodes=source_nodes,
            target_nodes=targets_nodes,
            lmax=lmax + 1,  # TODO: Fixme since enumerate_paths uses lmax + 1 (we have to now increase 1)
            simple_paths=simple_paths,
            output=output,
            genesets={**kegg, **reactome, **wikipathways},
        )
        # Finished time
        exe_t_f = time.time()
        running_time = exe_t_f - exe_t_0

        click.secho(f'{EMOJI} Finished in {running_time} seconds {EMOJI}')

        click.secho(f'{EMOJI} Results exported to {output} {EMOJI}')


def _setup_logging(log: bool, process_id: int = 0) -> None:
    """Set up logging.

    :param log: logging boolean
    :param process_id: process id
    """
    if process_id == 0:
        if log:
            logging.getLogger('drug2ways').setLevel(logging.DEBUG)
            logging.basicConfig(level=logging.DEBUG)
            click.echo('Debug mode on')
        else:
            logging.basicConfig(level=logging.INFO)
            logger.setLevel(logging.INFO)
    else:
        logging.getLogger('drug2ways').setLevel(logging.WARNING)
        logging.basicConfig(level=logging.WARNING)
