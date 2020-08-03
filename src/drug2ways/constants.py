# -*- coding: utf-8 -*-

"""Constants of drug2ways."""

import logging
import os
from urllib.request import urlretrieve

logger = logging.getLogger(__name__)

dir_path = os.path.dirname(os.path.realpath(__file__))
SOURCE_DIR = os.path.join(os.path.abspath(os.path.join(dir_path, os.pardir)))
ROOT_DIR = os.path.join(os.path.abspath(os.path.join(SOURCE_DIR, os.pardir)))

#: Default drug2ways directory
DEFAULT_DRUG2WAYS_DIR = os.path.join(os.path.expanduser('~'), '.drug2ways')

KEGG_GENESETS = os.path.join(DEFAULT_DRUG2WAYS_DIR, 'kegg.tsv')
REACTOME_GENESETS = os.path.join(DEFAULT_DRUG2WAYS_DIR, 'reactome.tsv')
WIKIPATHWAYS_GENESETS = os.path.join(DEFAULT_DRUG2WAYS_DIR, 'wp.tsv')

KEGG_GENESETS_URL = 'https://raw.githubusercontent.com/pathwayforte/pathway-forte/master/data/gmt_files/kegg.gmt'
REACTOME_GENESETS_URL = 'https://raw.githubusercontent.com/pathwayforte/pathway-forte/master/data/gmt_files/reactome.gmt'
WIKIPATHWAYS_GENESETS_URL = 'https://raw.githubusercontent.com/pathwayforte/pathway-forte/master/data/gmt_files/wikipathways.gmt'


def ensure_genesets():
    """Download gene sets."""
    logger.info('Downloading genesets for pathway predictions...')
    if not os.path.exists(KEGG_GENESETS):
        download_pathway(KEGG_GENESETS_URL, KEGG_GENESETS)
    if not os.path.exists(REACTOME_GENESETS):
        download_pathway(REACTOME_GENESETS_URL, REACTOME_GENESETS)
    if not os.path.exists(WIKIPATHWAYS_GENESETS):
        download_pathway(WIKIPATHWAYS_GENESETS_URL, WIKIPATHWAYS_GENESETS)


RESOURCES_DIR = os.path.join(ROOT_DIR, 'graphs')
RESULTS_DIR = os.path.join(ROOT_DIR, 'results')


def ensure_output_dirs():
    """Ensure that the output directories exists."""
    os.makedirs(DEFAULT_DRUG2WAYS_DIR, exist_ok=True)


def download_pathway(url: str, export_path: str) -> None:
    """Make a function that downloads the data for you, or uses a cached version at the given path.

    :param url: The URL of some data
    :param export_path: folder where decompressed file will be exported
    """
    urlretrieve(url, export_path)  # noqa: S310


ensure_output_dirs()

"""Available formats"""

#: csv
CSV = 'csv'
#: tsv
TSV = 'tsv'
#: graphML
GRAPHML = 'graphml'
#: bel
BEL = 'bel'
#: node link json
NODE_LINK_JSON = 'json'
#: pickle
BEL_PICKLE = 'pickle'
#: gml
GML = 'gml'
#: edge list
EDGE_LIST = '.lst'

#: drug2ways available network formats
FORMATS = [
    CSV,
    TSV,
    GRAPHML,
    BEL,
    NODE_LINK_JSON,
    BEL_PICKLE,
]

BEL_FORMATS = [
    BEL,
    BEL_PICKLE,
]

#: Separators
FORMAT_SEPARATOR_MAPPING = {
    CSV: ',',
    TSV: '\t'
}

"""Acceptable column names for the graph"""

#: Column name for source node
SOURCE = 'source'
#: Column name for target node
TARGET = 'target'
#: Column name for relation
RELATION = 'relation'

#: drug2ways emoji
EMOJI = "💊🔬"
