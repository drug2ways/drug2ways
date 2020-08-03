Drug2Ways |release| Documentation
=================================
Drug2ways supports generic network formats such as JSON, CSV, GraphML, or GML. Check out drug2ways's documentation here.
Ideally, the network should contain three different types of nodes representing drugs, proteins, and
indications/phenotypes. The hypothesis underlying this software is that by reasoning over a multitude of possible paths
between a given drug and indication, the drug regulates the indication in the direction of the signs of the most
frequently occurring paths (i.e., majority rule). In other words, we assume that a drug has a greater likelihood of
interacting with its target, and its target with intermediate nodes, to modulate a pathological phenotype as the number
of possible paths connecting a drug to the phenotype increases. Based on this hypothesis, this software can be applied
for different applications outlined in the next section.

Installation is as easy as getting the code from `PyPI <https://pypi.python.org/pypi/drug2ways>`_ with
:code:`python3 -m pip install drug2ways`. See the :doc:`installation <installation>` documentation.

.. seealso::

    - Documented on `Read the Docs <http://drug2ways.readthedocs.io/>`_
    - Versioned on `GitHub <https://github.com/drug2ways/drug2ways>`_
    - Tested on `Travis CI <https://travis-ci.org/drug2ways/drug2ways>`_
    - Distributed by `PyPI <https://pypi.python.org/pypi/drug2ways>`_

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   formats
   api
   cli
   algorithm
   constants
   performance


Disclaimer
----------
Drug2Ways is a scientific software that has been developed in an academic capacity, and thus comes with no warranty or
guarantee of maintenance, support, or back-up of data.
