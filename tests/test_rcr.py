# -*- coding: utf-8 -*-

"""Test rcr methods."""

import logging
import unittest
from networkx import DiGraph
import pandas as pd

from drug2ways.rcr import rcr_all_paths, validate_paths_with_disease_data

from tests.constants import RCR_GRAPH_EXAMPLE

log = logging.getLogger(__name__)

"""Helper functions for testing rcr"""


class TestRcr(unittest.TestCase):
    """RCR method test."""

    def create_graph(self):
        """Create example graph."""
        df = pd.read_csv(RCR_GRAPH_EXAMPLE)

        graph = DiGraph()

        for sub_name, obj_name, relation in df.values:
            # Store edge in the graph
            graph.add_edge(
                sub_name,
                obj_name,
                polarity=relation,
            )
        return graph

    def test_rcr_paths(self):
        """Test rcr method."""
        '''
            |--> B --> D --| |
        A - |                | - E
            |  --> C -->     |
        '''

        drug_dict = {'b': 1, 'd': 1, 'c': 1}

        paths = [['a', 'b', 'd', 'e'], ['a', 'c', 'e']]

        directed_graph = self.create_graph()

        filtered_paths = rcr_all_paths(
            graph=directed_graph,
            all_paths=paths,
            drug_dict=drug_dict
        )

        self.assertEqual(
            filtered_paths,
            [['a', 'b', 'd'], ['a', 'c']]
        )

    def test_rcr_paths_with_errors(self):
        """Test rcr method with error allowed."""
        drug_dict = {'b': 1, 'd': 1, 'c': -1}

        paths = [['a', 'b', 'd', 'e'], ['a', 'c', 'e']]

        directed_graph = self.create_graph()

        filtered_paths = rcr_all_paths(
            graph=directed_graph,
            all_paths=paths,
            drug_dict=drug_dict,
            errors_allowed=1
        )

        self.assertEqual(
            filtered_paths,
            [['a', 'b', 'd'], ['a', 'c']]
        )

    def test_validated_paths(self):
        """Test for validated paths."""
        disease_dict = {'b': -1, 'd': 1, 'c': -1}
        drug_dict = {'b': 1, 'd': 1, 'c': 1}

        paths = [['a', 'b', 'd'], ['a', 'c']]

        valid_paths = validate_paths_with_disease_data(
            paths=paths,
            drug_dict=drug_dict,
            disease_dict=disease_dict,
        )

        self.assertEqual(
            valid_paths,
            [['a', 'c']]
        )

    def test_validate_paths_with_errors(self):
        """Test for validated paths with errors allowed."""
        disease_dict = {'b': -1, 'd': 1, 'c': -1}
        drug_dict = {'b': 1, 'd': 1, 'c': 1}

        paths = [['a', 'b', 'd'], ['a', 'c']]

        valid_paths = validate_paths_with_disease_data(
            paths=paths,
            drug_dict=drug_dict,
            disease_dict=disease_dict,
            errors_allowed=1
        )

        self.assertEqual(
            valid_paths,
            [['a', 'b', 'd'], ['a', 'c']]
        )

    def test_rcr_pipeline(self):
        """Test rcr pipeline."""
        drug_dict = {'b': 1, 'd': 1, 'c': 1}
        disease_dict = {'b': -1, 'd': 1, 'c': -1}

        paths = [['a', 'b', 'd', 'e'], ['a', 'c', 'e']]

        directed_graph = self.create_graph()

        filtered_paths = rcr_all_paths(
            graph=directed_graph,
            all_paths=paths,
            drug_dict=drug_dict
        )

        valid_paths = validate_paths_with_disease_data(
            paths=filtered_paths,
            drug_dict=drug_dict,
            disease_dict=disease_dict,
        )

        self.assertEqual(
            valid_paths,
            [['a', 'c']]
        )

    def test_rcr_pipeline_with_errors(self):
        """Test rcr pipeline with error allowed."""
        drug_dict = {'b': 1, 'd': 1, 'c': -1}
        disease_dict = {'b': -1, 'd': 1, 'c': -1}

        paths = [['a', 'b', 'd', 'e'], ['a', 'c', 'e']]

        directed_graph = self.create_graph()

        filtered_paths = rcr_all_paths(
            graph=directed_graph,
            all_paths=paths,
            drug_dict=drug_dict,
            errors_allowed=1
        )

        valid_paths = validate_paths_with_disease_data(
            paths=filtered_paths,
            drug_dict=drug_dict,
            disease_dict=disease_dict,
            errors_allowed=1
        )

        self.assertEqual(
            valid_paths,
            [['a', 'b', 'd'], ['a', 'c']]
        )
