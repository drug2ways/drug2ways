# -*- coding: utf-8 -*-

"""Reverse causal reasoning on drug2ways."""

import logging
from itertools import tee
from typing import Any, Dict, List, Optional

from networkx import DiGraph

logger = logging.getLogger(__name__)


def pairwise(iterable):
    """Pairwise iteration."""
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def rcr_all_paths(
    graph: DiGraph,
    all_paths: List[List],
    drug_dict: Dict[str, int],
    errors_allowed: int = 0
) -> List[Any]:
    """Conduct causal reasoning on the paths between drug and disease based on drug experimental data.

    :param graph: original directed graph
    :param all_paths: all paths to evaluate
    :param drug_dict: dictionary with the fold changes from the drug experiment
    :param errors_allowed: errors allowed in the path
    """
    valid_paths = []

    for path in all_paths:

        # remove the last node since it is the disease and it doesnt have any experimental value
        path = path[:-1]

        # Skip path if not all nodes are present in experimental data
        # First node is not considered since it corresponds to the drug which doesnt have experimental value
        if not all(node in drug_dict for node in path[1:]):
            continue

        if not _is_concordant(graph, path, drug_dict, errors_allowed):
            continue

        valid_paths.append(path)

    return valid_paths


def _is_concordant(
    graph: DiGraph,
    path: List[str],
    drug_dict: Dict[str, int],
    errors_allowed: int = 0
) -> bool:
    """Calculate if the path is concordant.

    :param graph: original directed graph
    :param path: path to evaluate
    :param drug_dict: dictionary with the fold changes from the drug experiment
    :param errors_allowed: errors allowed in the path
    :return: boolean with the result
    """
    # Calculate the current score
    current_polarity = 1
    # number of errors during evaluation
    current_errors = 0

    for source, target in pairwise(path):

        # Update polarity
        current_polarity = current_polarity * graph.edges[source, target]['polarity']

        target_score = drug_dict[target]

        if current_polarity != target_score:
            # max errors allowed reached
            if current_errors == errors_allowed:
                return False
            # allow for one more error
            current_errors += 1

    return True


def validate_paths_with_disease_data(
    paths: List[List[str]],
    drug_dict: Dict[str, int],
    disease_dict: Dict[str, int],
    errors_allowed: int = 0
) -> List[Any]:
    """Validate paths with disease data.

    :param paths: validated paths
    :param drug_dict: path to evaluate
    :param disease_dict: dictionary with the fold changes from the drug experiment
    :param errors_allowed: errors allowed in the path
    :return: boolean with the result
    """
    # valid paths
    filtered_paths = []

    # check that nodes in paths have the opposite expression in drug and disease
    for path in paths:
        result = _evaluate_opposite_expression(
            path,
            disease_dict,
            drug_dict,
            errors_allowed
        )

        if not result:
            continue

        filtered_paths.append(path)

    return filtered_paths


def _evaluate_opposite_expression(
    path: List[str],
    drug_dict: Dict[str, int],
    disease_dict: Dict[str, int],
    errors_allowed: int = 0
) -> bool:
    """Evaluate opposite expression on a single path.

    :param path: path to evaluate
    :param drug_dict: path to evaluate
    :param disease_dict: dictionary with the fold changes from the drug experiment
    :param errors_allowed: errors allowed in the path
    :return: boolean with the result
    """
    current_errors = 0
    # Disease node at the end of the path has already been removed
    for node in path[1:]:

        if drug_dict[node] == disease_dict[node]:

            # max errors allowed reached
            if current_errors == errors_allowed:
                return False
            # allow for one more error
            current_errors += 1

    return True


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


def disease_rcr_all_paths(
    graph: DiGraph,
    all_paths: List[List],
    disease_dict: Dict[str, int],
    errors_allowed: int = 0
) -> List[Any]:
    """Conduct causal reasoning on the paths between drug and disease based on disease experimental data.

    :param graph: original directed graph
    :param all_paths: all paths to evaluate
    :param disease_dict: dictionary with the fold changes from the disease experiment
    :param errors_allowed: errors allowed in the path
    """
    valid_paths = []

    for path in all_paths:

        # remove the last node since it is the disease and it doesnt have any experimental value
        path = path[:-1]

        # Skip path if not all nodes are present in experimental data
        # First node is not considered since it corresponds to the drug which doesnt have experimental value
        if not all(node in disease_dict for node in path[1:]):
            continue

        if not _is_not_concordant(graph, path, disease_dict, errors_allowed):
            continue

        valid_paths.append(path)

    return valid_paths


def _is_not_concordant(
    graph: DiGraph,
    path: List[str],
    disease_dict: Dict[str, int],
    errors_allowed: int = 0
) -> bool:
    """Calculate if the path is concordant.

    :param graph: original directed graph
    :param path: path to evaluate
    :param disease_dict: dictionary with the fold changes from the disease experiment
    :param errors_allowed: errors allowed in the path
    :return: boolean with the result
    """
    # Calculate the current score
    current_polarity = 1
    # number of errors during evaluation
    current_errors = 0

    for source, target in pairwise(path):

        # Update polarity
        current_polarity = current_polarity * graph.edges[source, target]['polarity']

        target_score = disease_dict[target]

        if current_polarity == target_score:
            # max errors allowed reached
            if current_errors == errors_allowed:
                return False
            # allow for one more error
            current_errors += 1

    return True
