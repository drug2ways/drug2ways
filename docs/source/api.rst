Drug2ways API
=============
Drug2ways implements an API with multiple methods that facilitate its programmatic usage.
The majority of this methods are similar to the ones in the CLI. They are located in cli_helper.py.

.. code-block:: sh

    from drug2ways.cli_helper import wrapper_explore
    from networkx import DiGraph

    # Initialize a directed graph
    directed_graph = nx.DiGraph()
    directed_graph.add_edges_from([(1, 2), (1, 3)])

    results, time_cache = wrapper_explore(
       graph=directed_graph, # directed graph
       source_nodes=[1], # list of source nodes
       target_nodes=[3], # list of target nodes
       lmax=2, # max length of the path
       simple_paths=True, # with or without cycles (True=no cycles are allowed)
   )
