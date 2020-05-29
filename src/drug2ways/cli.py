# -*- coding: utf-8 -*-

"""Command line interface for drug2ways."""

import logging

import click

from drug2ways.cli_helper import _explore_helper, _optimize_helper, _combine_helper, _pathway_enrichment_helper
from .constants import DEFAULT_DRUG2WAYS_DIR, FORMATS, ensure_genesets

logger = logging.getLogger(__name__)


@click.group(help='drug2ways')
def main():
    """Command line interface for drug2ways."""
    logging.basicConfig(format="%(asctime)s - %(levelname)s - %(name)s - %(message)s")


input_graph_option = click.option(
    '-g', '--graph',
    help='Path to the network',
    required=True,
    type=click.Path(exists=True, dir_okay=False),
)
input_format_option = click.option(
    '-f', '--fmt',
    help='Graph fmt',
    type=click.Choice(FORMATS),
)
sources_option = click.option(
    '-s', '--sources',
    help='Path to file with source nodes',
    type=click.Path(exists=True, dir_okay=False),
)
targets_option = click.option(
    '-t', '--targets',
    help='Path to file with source nodes',
    type=click.Path(exists=True, dir_okay=False),
)
lmax_option = click.option(
    '-l', '--lmax',
    help='Maximum length of paths',
    required=True,
    type=int,
)
simple_option = click.option(
    '--simple',
    help="Count only simple paths, i.e. paths without cycles.",
    is_flag=True,
    default=False,
    type=bool,
)
output_option = click.option(
    '-o', '--output',
    help='Output directory',
    default=DEFAULT_DRUG2WAYS_DIR,
    type=click.Path(exists=True),
)
name_option = click.option(
    '-n', '--name',
    help='Name for output file',
    default='',
)
log_option = click.option('-l', '--log', is_flag=True, help='Activate debug mode')

time_option = click.option('-t', '--time', is_flag=True, help='Export time measurements')

threshold_option = click.option(
    '-a', '--activation-threshold',
    help='Activation threshold',
    type=float,
)

combination_length_option = click.option(
    '-c', '--combination-length',
    help='Combination length. Number of drugs in each combination.',
    required=True,
    type=int,
)

drug_search_option = click.option(
    '-d', '--drug-search-bel',
    help='Drug search for BEL graphs',
    is_flag=True,
)


@input_graph_option
@input_format_option
@sources_option
@targets_option
@lmax_option
@simple_option
@output_option
@name_option
@drug_search_option
@log_option
@time_option
@main.command()
def explore(
    graph: str,
    fmt: str,
    sources: str,
    targets: str,
    lmax: int,
    simple: bool,
    output: str,
    name: str,
    drug_search_bel: bool,
    log: bool = None,
    time: bool = None,
):
    """Run drug2ways for a given network."""
    _explore_helper(
        graph=graph,
        fmt=fmt,
        sources=sources,
        targets=targets,
        lmax=lmax,
        simple_paths=simple,
        log=log,
        export_time=time,
        name=name,
        output=output,
        drug_search_bel=drug_search_bel,
    )


@input_graph_option
@input_format_option
@sources_option
@targets_option
@lmax_option
@simple_option
@threshold_option
@output_option
@log_option
@main.command()
def optimize(
    graph: str,
    fmt: str,
    sources: str,
    targets: str,
    lmax: int,
    simple: bool,
    activation_threshold: float,
    output: str,
    log: bool = None,
):
    """Run drug2ways for a given network and get sources that optimize given targets."""
    _optimize_helper(
        graph=graph,
        fmt=fmt,
        sources=sources,
        targets=targets,
        lmax=lmax,
        simple_paths=simple,
        log=log,
        output=output,
        activation_threshold=activation_threshold,
    )


@input_graph_option
@input_format_option
@sources_option
@targets_option
@lmax_option
@simple_option
@threshold_option
@combination_length_option
@output_option
@log_option
@main.command()
def combine(
    graph: str,
    fmt: str,
    sources: str,
    targets: str,
    lmax: int,
    simple: bool,
    activation_threshold: float,
    combination_length: int,
    output: str,
    log: bool = None,
):
    """Run drug2ways for a given network and get sources that optimize given targets."""
    _combine_helper(
        graph=graph,
        fmt=fmt,
        sources=sources,
        targets=targets,
        lmax=lmax,
        simple_paths=simple,
        activation_threshold=activation_threshold,
        combination_length=combination_length,
        log=log,
        output=output,
    )


@input_graph_option
@input_format_option
@sources_option
@targets_option
@lmax_option
@simple_option
@output_option
@log_option
@main.command()
def pathway_analysis(
    graph: str,
    fmt: str,
    sources: str,
    targets: str,
    lmax: int,
    simple: bool,
    output: str,
    log: bool = None,
):
    """Run drug2ways pathway enrichment on the paths."""
    ensure_genesets()
    _pathway_enrichment_helper(
        graph=graph,
        fmt=fmt,
        sources=sources,
        targets=targets,
        lmax=lmax,
        simple_paths=simple,
        log=log,
        output=output,
    )
