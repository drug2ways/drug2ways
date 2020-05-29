# -*- coding: utf-8 -*-

"""Methods to generalize any graph for its path analysis."""

import logging
from typing import Dict, Tuple, Any, List

from networkx import DiGraph, isolates

__all__ = [
    'generate_reduced_graph',
]

logger = logging.getLogger(__name__)


def remove_isolated_nodes(graph: DiGraph):
    """Remove isolated nodes from the graph.

    :param graph: graph to be filtered
    """
    nodes = list(isolates(graph))
    graph.remove_nodes_from(nodes)


def _dict_to_graph(data: Dict[str, Any]) -> Tuple[DiGraph, Dict[str, int]]:
    """Convert dictionary representation of the graph to a directed graph.

    :param data: graph as a dictionary
    :return: directed graph
    """
    graph = DiGraph()
    node2id = {}
    for node, properties in data['node_list'].items():
        node2id[node] = properties['id']
        graph.add_node(
            int(properties['id']),
            name=node,
            isTarget=bool(properties['isTarget'])
        )

    for node, adj in data['adj_list'].items():
        source = int(node)
        increases = adj.get('increases', [])
        decreases = adj.get('decreases', [])

        for n in increases:
            graph.add_edge(source, n, polarity=1)
        for n in decreases:
            graph.add_edge(source, n, polarity=-1)

    return graph, node2id


def generate_reduced_graph(graph: DiGraph, target_nodes: List[Any]) -> Tuple[DiGraph, Dict[str, int]]:
    """Generate a reduced version of a graph.

    :param graph: directed graph
    :param target_nodes: target nodes
    :return:
    """
    remove_isolated_nodes(graph)

    node_list = {
        f'{node}': {
            'id': i,
            'isTarget': True if node in target_nodes else False,
        }
        for i, node in enumerate(graph.nodes())
    }

    adj_list = {}
    # Counters
    num_edges = 0
    count_increases = 0
    count_decreases = 0

    for i, node in enumerate(graph.nodes()):
        increases = []
        decreases = []
        for neighbor in graph.neighbors(node):
            relation_sign = graph[node][neighbor].get('relation')

            if not relation_sign:
                raise ValueError('Ensure that your graph has been loaded within the "relation" attribute')

            # Add positive relation
            if relation_sign == 1:
                increases.append(node_list[f'{neighbor}']['id'])
                count_increases += 1

            # Add negative relation
            elif relation_sign == -1:
                decreases.append(node_list[f'{neighbor}']['id'])
                count_decreases += 1

            # Raise error if it doesnt recognize the relation type
            else:
                ValueError(f"Unknown relation: {relation_sign}")

        if increases or decreases:
            adj_list[i] = {}

        if increases:
            adj_list[i]['increases'] = increases

        if decreases:
            adj_list[i]['decreases'] = decreases

        num_edges += len(increases) + len(decreases)

    num_nodes = len(node_list)

    graph_data = {
        'num_nodes': num_nodes,
        'num_edges': num_edges,
        'node_list': node_list,
        'adj_list': adj_list
    }

    logger.debug(
        f"Number of nodes:{num_nodes}\n"
        f"Number of edges: {num_edges}\n"
        f"Number of activations: {count_increases}\n"
        f"Number of inhibitions:  {count_decreases}\n"
    )

    return _dict_to_graph(graph_data)
