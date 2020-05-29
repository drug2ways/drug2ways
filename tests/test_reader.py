# -*- coding: utf-8 -*-

"""Test read graphs."""

import logging
import unittest

from drug2ways.graph_reader import load_graph

from tests.constants import GRAPH_EXAMPLE

log = logging.getLogger(__name__)

"""Helper functions for reading graphs"""


class TestReader(unittest.TestCase):
    """Read search test."""

    def test_read_graph(self):
        """Test read graph."""
        directed_graph = load_graph(GRAPH_EXAMPLE, 'csv')

        edges = list(directed_graph.edges.data())

        self.assertEqual(
            edges,
            [('a', 'b', {'relation': 1}), ('a', 'c', {'relation': 1}), ('b', 'd', {'relation': 1}),
             ('d', 'e', {'relation': 1}), ('d', 'g', {'relation': 1}), ('e', 'b', {'relation': 1}),
             ('c', 'f', {'relation': 1}), ('f', 'g', {'relation': -1})]
        )
