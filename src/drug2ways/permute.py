# -*- coding: utf-8 -*-

"""Utils for permuting networks."""

import logging
import os
from collections import Counter

import numpy as np
import pandas as pd
from xswap.permute import permute_edge_list
from xswap.preprocessing import load_str_edges, map_str_edges

from drug2ways.constants import ROOT_DIR

logger = logging.getLogger(__name__)

__all__ = [
    'permute_network',
]


def permute_network(path: str, sep: str) -> pd.DataFrame:
    """Permute network.

    :param path: directory to the network file
    :return: dataframe with permuted network
    """
    # Read network
    network_df = pd.read_csv(path, sep=sep)
    # Read relations
    relations = network_df.relation
    logger.info(Counter(relations))

    # Read edge list
    edge_list = load_str_edges(NETWORK_PATH, node_delim=sep)
    # Get mapping since the edge list contains now integers
    edge_list_integers, node_mapping, _ = map_str_edges(edge_list, bipartite=False)
    # Permute
    permuted_edges, stats = permute_edge_list(edge_list_integers)

    logger.info(stats)

    # Reverse mapping dictionary
    node_mapping = {
        v: k
        for k, v in node_mapping.items()
    }

    logger.info("Shuffling edges...")

    np.random.permutation(relations),  # Get random edges

    return pd.DataFrame(
        {
            'source': node_mapping[source],
            'target': node_mapping[target],
            'relation': relations[index],
        }
        for index, (source, target) in enumerate(permuted_edges)
        if index < len(relations)  # Skip the header if has been read by xswap
    )


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    # TODO: Choose one of the two
    NETWORK_FILES = {
        'openbiolink': "openbiolink_network.tsv",
        'custom': "custom_network.tsv",
    }

    NETWORK_NAME = NETWORK_FILES['custom']

    # Path to the network
    NETWORK_PATH = os.path.join(ROOT_DIR, "data", "networks", "data", NETWORK_NAME)

    # Permute
    shuffled_network = permute_network(NETWORK_PATH, '\t')

    logger.info(Counter(shuffled_network.relation))

    # Export
    shuffled_network.to_csv(f"shuffled_{NETWORK_NAME}", sep='\t', index=False)
