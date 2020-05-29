# -*- coding: utf-8 -*-

"""Pathway analysis methods."""

import logging
from collections import Counter, defaultdict
from typing import Dict, List, Mapping

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

from .constants import KEGG_GENESETS, REACTOME_GENESETS, WIKIPATHWAYS_GENESETS

logger = logging.getLogger(__name__)


def _prepare_json(paths):
    """Prepare json."""
    return {
        index: [
            cosa
            for cosa in path
        ]
        for index, path in enumerate(paths)
    }


def analyze_paths(
    reduced_graph,
    paths: List[List[int]],
    id2node: Mapping[int, str],
    lmax: int,
    target_id: int,
    genesets: List[Mapping[str, List[str]]],
    min_count: int = 0,
    min_proportion: int = 0,
):
    """Analyze paths.

    :param paths: paths
    :param id2node: mapping between ids and node names
    :param lmax: lmax
    :param target_id: id of the target
    :param genesets: pathway genesets
    :param min_count: minimum number of times at a given lmax
    :param min_proportion: minimum proportion of that node based on the total
    """
    results = defaultdict(Counter)

    # Iter over all paths
    for path in paths:
        # Iterate through each node in the path while tracking its position in the path
        for index, node in enumerate(path):
            results[index][id2node[node]] += 1

    polarity_dict = {1: '->', -1: '-|'}

    final_paths = set()
    for path in paths:

        reconstructed_path = []

        for index, node in enumerate(path):

            # Avoid crashing
            if index + 1 == len(path):
                continue

            polarity = polarity_dict[reduced_graph[node][path[index + 1]]['polarity']]

            if index == 0:
                reconstructed_path.append(node)
                reconstructed_path.append(polarity)
                reconstructed_path.append(path[index + 1])
            else:
                reconstructed_path.append(polarity)
                reconstructed_path.append(path[index + 1])

        final_paths.add(
            tuple(
                id2node[cosa] if cosa in id2node else cosa
                for cosa in reconstructed_path
            )
        )

    final_paths = _prepare_json(final_paths)

    df = pd.DataFrame()

    for lmax, counter in results.items():
        # Total number of nodes (incl. duplicates) on that path position
        total = sum(counter.values())

        df.loc[:, lmax] = pd.Series([
            f"{node} ({count})"
            for node, count in counter.most_common()
            if count > min_count and (count * 100) / total > min_proportion
            # Threshold on absolute count and proportion
        ])

    # TODO: Change and all 3 databases
    # enrichment_results  = pathway_enrichment(df, genesets[0])

    return df, final_paths


def _sanitize_name(node, prefix):
    """Sanitize name."""
    return node.split('(')[0].strip().replace(prefix, '')


def pathway_enrichment(df: pd.DataFrame, geneset, prefix: str = 'NCBIGENE:'):
    """Enrich pathways on each lmax."""
    enrichment_results = {}
    # Iterate over columns
    for lmax in df:
        nodes = {
            _sanitize_name(node, prefix)
            for node in df[lmax]
            if pd.notna(node)
        }

        enrichment_results[lmax] = perform_hypergeometric_test(
            genes_to_test=nodes,
            pathway_dict=geneset,
        )

    return enrichment_results


def get_genesets():
    """Get gene sets as dicts."""
    return (
        parse_gmt_file(KEGG_GENESETS),
        parse_gmt_file(REACTOME_GENESETS),
        parse_gmt_file(WIKIPATHWAYS_GENESETS),
    )


def parse_gmt_file(gmt_file: str, min_size=3, max_size=3000) -> Dict[str, List]:
    """Parse gmt file."""
    with open(gmt_file, 'r') as file:
        geneset_dict = {
            line.strip().split("\t")[0]: line.strip().split("\t")[2:]
            for line in file.readlines()
        }
    return {
        k: v for k, v in geneset_dict.items() if len(v) >= min_size and len(v) <= max_size
    }


def _prepare_hypergeometric_test(
    query_gene_set,
    pathway_gene_set,
    gene_universe,
):
    """Prepare the matrix for hypergeometric test calculations.

    :param query_gene_set: gene set to test against pathway
    :param pathway_gene_set: pathway gene set
    :param gene_universe: number of HGNC symbols
    :return: 2x2 matrix
    """
    # Cast lists to sets
    if not isinstance(query_gene_set, set):
        query_gene_set = set(query_gene_set)
    if not isinstance(pathway_gene_set, set):
        pathway_gene_set = set(pathway_gene_set)

    # Return matrix to test hyper-geometric test
    return np.array([
        [
            len(query_gene_set.intersection(pathway_gene_set)),
            len(query_gene_set.difference(pathway_gene_set)),
        ],
        [
            len(pathway_gene_set.difference(query_gene_set)),
            gene_universe - len(pathway_gene_set.union(query_gene_set)),
        ],
    ])


def perform_hypergeometric_test(
    genes_to_test,
    pathway_dict,
    gene_universe: int = 41714,
    apply_threshold=True,
    threshold=0.05,
):
    """Perform hypergeometric tests.

    :param genes_to_test: gene set to test against pathway
    :param pathway_dict: pathway name to gene set
    :param gene_universe: number of HGNC symbols
    :param apply_threshold: return only significant pathways
    :param threshold: significance threshold (by default 0.05)
    """
    rows = []
    for pathway_id, pathway_gene_set in pathway_dict.items():
        # Prepare the test table to conduct the fisher test
        test_table = _prepare_hypergeometric_test(genes_to_test, pathway_gene_set, gene_universe)
        # Calculate fisher test (returns tuple of odds ratio and p_value
        p_value = fisher_exact(test_table, alternative='greater')[1]
        rows.append((pathway_id, p_value))

    df = pd.DataFrame(rows, columns=['pathway_id', 'pval'])
    correction_test = multipletests(df.pval, method='fdr_bh')
    df['qval'] = correction_test[1]

    if apply_threshold:
        logger.debug(f'Filtering out pathways with q-values > {threshold} according to fdr_bh')
        df = df[df['qval'] < threshold]

    return df
