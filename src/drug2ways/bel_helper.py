# -*- coding: utf-8 -*-

"""Biological Reasoner of drug2ways."""

import logging
from typing import Dict, Tuple, List, Any

from networkx import DiGraph
from pybel import BELGraph
from pybel.constants import (
    RELATION, INCREASES, DIRECTLY_DECREASES, DIRECTLY_INCREASES, DECREASES, REGULATES
)
from pybel.struct.filters import is_causal_relation
from pybel.struct.mutation import get_subgraph_by_edge_filter

logger = logging.getLogger(__name__)

RELATION_MAPPING_BEL = {
    INCREASES: 1,
    DIRECTLY_INCREASES: 1,
    DECREASES: -1,
    DIRECTLY_DECREASES: -1,
    REGULATES: 1,
}


def _change_relationships(edge: Dict) -> Tuple[bool, bool]:
    """Validate relationship."""
    if 'increases' in edge[1]['relation'] or edge[1]['relation'] == 'positive_correlation':
        return True, True
    elif 'decreases' in edge[1]['relation'] or edge[1]['relation'] == 'negative_correlation':
        return True, False

    return False, False


def process_bel_graph(bel_graph: BELGraph) -> DiGraph:
    """Convert BEL Graph to a directed graph ready for drug2ways.

    :param bel_graph: BELGraph
    :return: directed graph
    """
    # Get only causal edges
    bel_graph = get_subgraph_by_edge_filter(bel_graph, is_causal_relation)

    directed_graph = DiGraph()
    for source, target, data in bel_graph.edges(data=True):
        if data[RELATION] not in RELATION_MAPPING_BEL:
            logger.warning(f"Unknown relation {data[RELATION]}")
            continue
        directed_graph.add_edge(source.as_bel(), target.as_bel(), relation=RELATION_MAPPING_BEL[data[RELATION]])

    return directed_graph


def _is_target_node(node: str) -> bool:
    """Check if it is valid target node in BEL.

    :param node: string representing the node
    :return: boolean checking whether the node is a valid target in BEL
    """
    if node.startswith('bp') or node.startswith('path'):
        return True
    return False


def _valid_source_node(node: str) -> bool:
    """Check if it is valid source node in BEL.

    :param node: string representing the node
    :return: boolean checking whether the node is a valid target in BEL
    """
    # Check that it is an abundance
    if not node.startswith('a'):
        return False
    # check that the namespace is CHEBI and PUBMED
    if 'CHEBI' in node or 'PUBCHEM' in node:
        return True

    return False


def get_candidate_drugs(graph):
    """Return all candidate drugs on a given BEL graph."""
    return [
        node
        for node in graph.nodes()
        if _valid_source_node(node)
    ]


def get_candidate_targets(graph):
    """Return all candidate target on a given BEL graph."""
    return [
        node
        for node in graph.nodes()
        if _is_target_node(node)
    ]


def remove_contradictory_edges(
    increases: List[Any],
    decreases: List[Any],
    debug: bool = False,
) -> Tuple[List[Any], List[Any], List[Any]]:
    """Remove contradictory edges."""
    removed = []
    for node in increases:
        if debug:
            logger.warning(f'Node: {node}')
        if node in decreases:
            if debug:
                logger.warning(f'{node} is contradictory')
            removed.append(node)
        else:
            if debug:
                logger.warning(f'{node} not contradictory')

    increases = [
        i
        for i in increases
        if i not in removed
    ]
    decreases = [
        i
        for i in decreases
        if i not in removed
    ]

    return increases, decreases, removed
