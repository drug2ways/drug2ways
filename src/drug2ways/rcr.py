# -*- coding: utf-8 -*-

"""Reverse causal reasoning on drug2ways."""

import logging
from itertools import tee
from typing import Any, Dict, List

from networkx import DiGraph

logger = logging.getLogger(__name__)


def pairwise(iterable):
    """Pairwise iteration."""
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def rcr_all_paths(graph: DiGraph, all_paths: List[List], drug_dict: Dict[str, int]) -> List[Any]:
    """Conduct causal reasoning on the paths between drug and disease based on drug experimental data."""
    valid_paths = []

    for path in all_paths:

        # remove the last node since it is the disease and it doesnt have any experimental value
        path = path[:-1]

        # Skip path if not all nodes are present in experimental data
        # First node is not considered since it corresponds to the drug which doesnt have experimental value
        if not all(node in drug_dict for node in path[1:]):
            continue

        if not _is_concordant(graph, path, drug_dict):
            continue

        valid_paths.append(path)

    return valid_paths


def _is_concordant(graph: DiGraph, path: List[str], drug_dict: Dict[str, int]) -> bool:
    """Calculate if the path is concordant.

    :param graph: original directed graph
    :param path: path to evaluate
    :param drug_dict: dictionary with the fold changes from the drug experiment
    :return: boolean with the result
    """
    # Calculate the current score
    current_polarity = 1
    for source, target in pairwise(path):
        # Update polarity
        current_polarity = current_polarity * graph.edges[source, target]['polarity']

        target_score = drug_dict[target]

        if current_polarity != target_score:
            return False

    return True


def validate_paths_with_disease_data(
    paths: List[str],
    drug_dict: Dict[str, int],
    disease_dict: Dict[str, int],
) -> List[Any]:
    """Validate paths with disease data.

    :param paths: validated paths
    :param drug_dict: path to evaluate
    :param disease_dict: dictionary with the fold changes from the drug experiment
    :return: boolean with the result
    """
    # Calculate the current score
    filtered_paths = []
    for path in paths:
        for node in path[1:]:

            # Evaluate
            if drug_dict[node] == disease_dict[node]:
                continue

        filtered_paths.append(path)

    return filtered_paths


def evaluate_data_network_overlap(
    graph: DiGraph,
    data_dictionary: Dict[str, int],
):
    """Check overlap."""
    node_names = set(graph.nodes())
    data_names = set(data_dictionary.keys())

    overlap = node_names.intersection(data_names)

    if not overlap:
        raise ValueError('Your data does not match to any node in the network')

    logger.info(
        f'Number of nodes in the network: {len(node_names)}'
        f'Number of genes/proteins measured in the dataset: {len(data_names)}'
        f'Overlap between data and network: {len(overlap)}'
    )
