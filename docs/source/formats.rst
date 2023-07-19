Formats
=======
Drug2ways handles the following formats:

- CSV (.csv)
- TSV (.tsv)
- GraphML_ (.graphml or .xml)
- GML_ (.gml or .xml)
- BEL_ (.bel)
- Pickle (.pickle). BELGraph object from PyBEL_ 0.13.2
- Edge_ `list`__ (.lst)
- Node-Link JSON (.json)

The GraphML, GML, EdgeList and Node-Link JSON readers of drug2ways rely on NetworkX. BEL and Pickle files rely on PyBEL.
The two easiest format to work with drug2ways are a triple-based file represented as the file below (tsv/csv).
The file must contain three columns: source, polarity and target (order is not relevant) and the only condition is that
the polarity column contains 1 and -1 to indicate the direction of the polarity (increase/decrease).

+----------+---------+----------+
| source   | target  | polarity |
+==========+=========+==========+
|  Drug1   |Protein1 |    -1    |
+----------+---------+----------+
| Protein1 |Protein2 |     1    |
+----------+---------+----------+
| Protein2 |Protein3 |    -1    |
+----------+---------+----------+
| Protein3 |Disease1 |     1    |
+----------+---------+----------+

.. _Edge: https://networkx.github.io/documentation/stable/reference/readwrite/edgelist.html
__ Edge_
.. _GraphML: http://graphml.graphdrawing.org
.. _BEL: https://language.bel.bio/
.. _GML: http://docs.yworks.com/yfiles/doc/developers-guide/gml.html
.. _PyBEL: https://github.com/pybel/pybel/
