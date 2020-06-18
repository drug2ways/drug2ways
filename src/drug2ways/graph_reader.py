# -*- coding: utf-8 -*-

"""Methods to generalize any graph for its path analysis."""

import json
import logging
from collections import Counter
from typing import Any, Dict, List, Tuple

import pandas as pd
from networkx import DiGraph, read_graphml, read_gml, node_link_graph, read_edgelist, connected_components
from pybel.io.gpickle import from_pickle
from pybel.io.lines import from_bel_script
from tqdm import tqdm

from .bel_helper import process_bel_graph
from .constants import (
    EMOJI,
    FORMATS,
    CSV,
    TSV,
    GRAPHML,
    GML,
    BEL_PICKLE,
    BEL,
    EDGE_LIST,
    NODE_LINK_JSON,
    SOURCE,
    TARGET,
    RELATION,
    FORMAT_SEPARATOR_MAPPING,
)

__all__ = [
    'load_graph',
    'read_nodes',
    'read_nodes_to_optimize',
]

logger = logging.getLogger(__name__)


def read_nodes(path: str) -> List[Any]:
    """Read node names from lines.

    :param path: path to file
    :return: list of nodes
    """
    with open(path) as f:
        content = f.readlines()
    return [
        line.strip()
        for line in content
    ]


def read_nodes_to_optimize(path: str) -> Tuple[List[Any], List[Any]]:
    """Read target node names to optimize from lines.

    :param path: path to file
    :return: Tuple containing nodes to activate and inhibit, respectively
    """
    activate = []
    inhibit = []
    with open(path) as f:
        content = f.readlines()

        for line in content:
            target = line.rstrip().split(',')
            if target[1] in ('1', 'activate'):
                activate.append(target[0])
            elif target[1] in ('-1', 'inhibit'):
                inhibit.append(target[0])
            else:
                raise IOError(
                    f'{EMOJI} Invalid format for targets to optimize in file {path}.'
                )

    return activate, inhibit


def load_graph(path: str, fmt: str) -> DiGraph:
    """Load graph from file.

    :param path: path to the graph
    :param fmt: graph fmt
    :return: directed graph
    """
    if fmt == CSV:
        return process_network(path, CSV)

    elif fmt == TSV:
        return process_network(path, TSV)

    elif fmt == GRAPHML:
        return read_graphml(path)

    elif fmt == GML:
        return read_gml(path)

    elif fmt == BEL:
        bel_graph = from_bel_script(path)
        return process_bel_graph(bel_graph)

    elif fmt == BEL_PICKLE:
        bel_graph = from_pickle(path)
        return process_bel_graph(bel_graph)

    elif fmt == EDGE_LIST:
        return read_edgelist(path)

    elif fmt == NODE_LINK_JSON:
        data = load_json_file(path)
        return node_link_graph(data)

    raise IOError(
        f'{EMOJI} The fmt used on {path} is not valid. Please ensure you use one of the following formats: ',
        f'{FORMATS}',
    )


"""Check formats of networks """


def _harmonize_relation_column(df: pd.DataFrame) -> pd.Series:
    """Convert to integers the relation column.

    :param df: network dataframe
    :return: harmonized column
    """
    return df[RELATION].astype(str).astype(int)


def _format_checker(fmt: str) -> None:
    """Check column sep."""
    if fmt not in FORMAT_SEPARATOR_MAPPING:
        raise ValueError(
            f'The selected sep {fmt} is not valid. Please ensure you use one of the following formats: ',
            f'{FORMATS}',
        )


"""Process networks"""


def _read_network_file(path: str, fmt: str) -> pd.DataFrame:
    """Read network file."""
    _format_checker(fmt)

    df = pd.read_csv(
        path,
        sep=FORMAT_SEPARATOR_MAPPING[CSV] if fmt == CSV else FORMAT_SEPARATOR_MAPPING[TSV],
        dtype=str,
        usecols=['source', 'target', 'relation'],
    )

    if SOURCE not in df.columns or TARGET not in df.columns or RELATION not in df.columns:
        raise ValueError(
            f'Ensure that your file contains columns for {SOURCE}, {TARGET} and {RELATION}',
        )

    # Ensure quality relations
    try:
        df[RELATION] = _harmonize_relation_column(df)
    except ValueError:
        # Remove shitty relations
        original_relations = dict(Counter(df.relation))
        df = df[df[RELATION].isin([1, -1, '1', '-1'])]

        # Try to harmonize again
        df[RELATION] = _harmonize_relation_column(df)

        # Report to the user things that have been removed
        logger.warning(
            f'Drug2Ways have found not valid relations on your network: {original_relations}.'
            f'This relationships have been removed. Remember that only "1" and "-1" is allowed in the relation column'
            f'Relations remaining: {dict(Counter(df[RELATION]))}.'
        )

    return df


def load_json_file(path: str) -> Dict:
    """Read json file."""
    with open(path) as f:
        return json.load(f)


def process_network(path: str, sep: str, connectivity: bool = False) -> DiGraph:
    """Return network from dataframe."""
    _format_checker(sep)

    df = _read_network_file(path, sep)

    graph = DiGraph()

    for sub_name, obj_name, relation in tqdm(df.values, total=df.shape[0], desc='Loading graph'):
        # Store edge in the graph
        graph.add_edge(
            sub_name,
            obj_name,
            relation=relation,
        )

    logger.debug(f"Report on the number of relations: {dict(Counter(df.relation))}")

    if connectivity:
        cc = [
            len(c)
            for c in sorted(connected_components(graph.to_undirected()), key=len, reverse=True)
        ]
        logger.info(f"Connected components size: {cc}")

    return graph
