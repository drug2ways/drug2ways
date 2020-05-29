# -*- coding: utf-8 -*-

"""Test path search algorithms."""

import logging
import unittest

from drug2ways.alternative_graph_traversal import enumerate_paths
from drug2ways.graph_processing import generate_reduced_graph
from drug2ways.graph_reader import load_graph, read_nodes
from drug2ways.wrapper import wrapper_explore
from tests.constants import GRAPH_EXAMPLE, SOURCES_EXAMPLE, TARGETS_EXAMPLE

log = logging.getLogger(__name__)

"""Helper functions for testing."""


class TestSearch(unittest.TestCase):
    """Path search test."""

    def test_path_cycles_1(self):
        """Test path search with cycles."""
        directed_graph, source_nodes, target_nodes = self._load_helper()

        # Call main function
        results, _ = wrapper_explore(
            graph=directed_graph,
            source_nodes=source_nodes,
            target_nodes=target_nodes,
            lmax=6,
            simple_paths=False,
        )

        self.assertEqual(
            [
                {'source': 'a',
                 'target': 'g',
                 'relative_activation': 0.67,
                 'relative_inhibition': 0.33,
                 'number_of_paths': 3,
                 }
            ],
            results,
        )

    def test_path_cycles_2(self):
        """Test path search without cycles."""
        directed_graph = load_graph(GRAPH_EXAMPLE, 'csv')

        # Get source nodes from file
        source_nodes = ['b']
        # Get target nodes from file
        targets_nodes = read_nodes(TARGETS_EXAMPLE)

        # Call main function
        results, _ = wrapper_explore(
            graph=directed_graph,
            source_nodes=source_nodes,
            target_nodes=targets_nodes,
            lmax=5,
            simple_paths=False,
        )

        self.assertEqual(
            [
                {'source': 'b',
                 'target': 'g',
                 'relative_activation': 1,
                 'relative_inhibition': 0,
                 'number_of_paths': 2,
                 }
            ],
            results,
        )

    def test_simple_paths_1(self):
        """Test simple paths search reaching the target."""
        directed_graph, source_nodes, target_nodes = self._load_helper()

        # Call main function
        results, _ = wrapper_explore(
            graph=directed_graph,
            source_nodes=source_nodes,
            target_nodes=target_nodes,
            lmax=3,
            simple_paths=True,
        )

        self.assertEqual(
            [
                {'source': 'a',
                 'target': 'g',
                 'relative_activation': 0.5,
                 'relative_inhibition': 0.5,
                 'number_of_paths': 2,
                 }
            ],
            results,
        )

    def test_simple_paths_2(self):
        """Test simple paths search not reaching the target."""
        directed_graph, source_nodes, target_nodes = self._load_helper()

        # Call main function
        results, _ = wrapper_explore(
            graph=directed_graph,
            source_nodes=source_nodes,
            target_nodes=target_nodes,
            lmax=2,
            simple_paths=True,
        )

        self.assertFalse(results)

    def test_enumerate_paths_1(self):
        """Test enumerate paths two paths."""
        directed_graph, source_nodes, target_nodes = self._load_helper()

        reduced_graph, node2id = generate_reduced_graph(directed_graph, target_nodes)
        _target_nodes = [node2id[target_node] for target_node in target_nodes]

        # Call main function
        results = enumerate_paths(
            graph=reduced_graph,
            source=node2id[source_nodes[0]],
            targets=_target_nodes,
            lmax=7,  # enumerate_paths works with lmax+1 so this is lmax=3
            cycle_free=False,
        )

        # Path from a to g
        self.assertEqual(
            [[0, 1, 2, 4], [0, 5, 6, 4], [0, 1, 2, 3, 1, 2, 4]],
            results,
        )

    def test_enumerate_paths_2(self):
        """Test enumerate paths one path."""
        directed_graph, source_nodes, target_nodes = self._load_helper()

        reduced_graph, node2id = generate_reduced_graph(directed_graph, target_nodes)
        _target_nodes = [node2id[target_node] for target_node in target_nodes]

        # Call main function
        results = enumerate_paths(
            graph=reduced_graph,
            source=node2id[source_nodes[0]],
            targets=_target_nodes,
            lmax=4,  # enumerate_paths works with lmax+1 so this is lmax=3
            cycle_free=True,
        )

        # Path from a to g
        self.assertEqual(
            [[0, 1, 2, 4], [0, 5, 6, 4]],
            results,
        )

    def _load_helper(self):
        directed_graph = load_graph(GRAPH_EXAMPLE, 'csv')
        # Get source nodes from file
        source_nodes = read_nodes(SOURCES_EXAMPLE)
        # Get target nodes from file
        target_nodes = read_nodes(TARGETS_EXAMPLE)
        return directed_graph, source_nodes, target_nodes
