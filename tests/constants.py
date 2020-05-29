# -*- coding: utf-8 -*-

"""Tests constants."""
import os

TEST_FOLDER = os.path.dirname(os.path.realpath(__file__))
RESOURCES_FOLDER = os.path.join(TEST_FOLDER, 'resources')

GRAPH_EXAMPLE = os.path.join(RESOURCES_FOLDER, 'graph.csv')
SOURCES_EXAMPLE = os.path.join(RESOURCES_FOLDER, 'sources.txt')
TARGETS_EXAMPLE = os.path.join(RESOURCES_FOLDER, 'targets.txt')
