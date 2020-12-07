<p align="center">
  <a href="https://drug2ways.readthedocs.io/en/latest">
     <img src="docs/source/meta/logo1.jpg" height="300">
  </a>
</p>

<h1 align="center">
  drug2ways
</h1>

<p align="center">
  <a href="https://travis-ci.com/drug2ways/drug2ways">
    <img src="https://travis-ci.com/drug2ways/drug2ways.svg?branch=master"
         alt="Travis CI">
  </a>

  <a href='https://opensource.org/licenses/Apache-2.0'>
    <img src='https://img.shields.io/badge/License-Apache%202.0-blue.svg' alt='License'/>
  </a>

  <a href="https://zenodo.org/badge/latestdoi/267315762">
    <img src="https://zenodo.org/badge/267315762.svg" alt="DOI">
  </a>

</p>

<p align="center">
    <b>Drug2ways</b> is a Python package for reasoning over paths on biological networks for drug discovery
</p>

<p align="center">
  <a href="#quickstart">Quickstart</a> •
  <a href="#applications">Applications</a> •
  <a href="#installation">Installation</a>
</p>


## Quickstart
Drug2ways supports generic network formats such as JSON, CSV, GraphML, or GML. Check out [drug2ways's documentation here](https://drug2ways.readthedocs.io/en/latest). Ideally, the network should contain three different types of nodes
representing drugs, proteins, and indications/phenotypes. The hypothesis underlying this software is that by reasoning
over a multitude of possible paths between a given drug and indication, the drug regulates the indication in the
direction of the signs of the most frequently occurring paths (i.e., majority rule). In other words, we assume that a
drug has a greater likelihood of interacting with its target, and its target with intermediate nodes, to modulate a
pathological phenotype as the number of possible paths connecting a drug to the phenotype increases. Based on this
hypothesis, this software can be applied for different applications outlined in the next section.

### Citation
If you use drug2ways for your research please cite our [paper](https://doi.org/10.1371/journal.pcbi.1008464): 

> Daniel Rivas-Barragan, Sarah Mubeen, Francesc Guim-Bernat,Martin Hofmann-Apitius, and Daniel Domingo-Fernández (2020).
Drug2ways: Reasoning over causal paths in biological networks for drug discovery. *PLOS Computational Biology* 16(12): e1008464;  https://doi.org/10.1371/journal.pcbi.1008464

## Applications
Drug2ways can be applied for three different applications:

**Scripts and real examples**: https://github.com/drug2ways/drug2ways/tree/master/examples

### 1. Identifying candidate drugs
The following command of the command line interface (CLI) of drug2ways enables candidate drug identification. The
minimum required input are the path to the network and its format, a path to the nodes considered as drugs and the
ones considered as conditions/phenotypes. Finally, the maximum length allowed for a given path (i.e., lmax). Type
"python -m drug2ways explore --help" to see other optional arguments.

```python
python -m drug2ways explore \
       --graph=<path-to-graph> \
       --fmt=<format> \
       --sources=<sources> \
       --targets=<targets> \
       --lmax=<lmax>
```

### 2. Optimization of drugs' effects
The following command of the CLI of drug2ways enables searching drugs that not only target a given disease but also
activate/inhibit a set of phenotypes. This method requires the same arguments as the previous explore functionality
but the target file requires an additional second column where the desired effect on the node (e.g., 'node1,activate')
is specified. See the examples directory for more information.

```python
python -m drug2ways optimize \
       --graph=<path-to-graph> \
       --fmt=<format> \
       --sources=<sources> \
       --targets=<targets> \ # Note that this file is slightly different than the other targets
       --lmax=<lmax>
```

### 3. Proposing combination therapies
The following command of the CLI of drug2ways enables the identification of candidate drugs for combination therapies.
The minimum required input are the path to the network and its format, a path to the nodes considered as drugs and the
ones considered as conditions/phenotypes. As with the optimization command, here again the target file requires an
additional second column specifying the desired effect on the node (e.g., 'node1,activate'). Furthermore, the maximum
length allowed for a given path (i.e., lmax) and the possible number of combinations of drugs must be provided. Type
"python -m drug2ways combine --help" to see other optional arguments.

```python
python -m drug2ways combine \
       --graph=<path-to-graph> \
       --fmt=<format> \
       --sources=<sources> \
       --targets=<targets> \
       --lmax=<lmax> \
       --combination-length=<number>
```

## Installation

<p align="center">
  <a href="https://drug2ways.readthedocs.io/en/latest/">
    <img src="http://readthedocs.org/projects/drug2ways/badge/?version=latest"
         alt="Documentation">
  </a>

  <img src='https://img.shields.io/pypi/pyversions/drug2ways.svg' alt='Stable Supported Python Versions'/>
  
  <a href="https://pypi.python.org/pypi/drug2ways">
    <img src="https://img.shields.io/pypi/pyversions/drug2ways.svg"
         alt="PyPi">
  </a>
</p>

The latest stable code can be installed from [PyPI](https://pypi.python.org/pypi/drug2ways) with:

```python
python -m pip install drug2ways
```

The most recent code can be installed from the source on [GitHub](https://github.com/drug2ways/drug2ways) with:

```python
python -m pip install git+https://github.com/drug2ways/drug2ways.git
```

For developers, the repository can be cloned from [GitHub](https://github.com/drug2ways/drug2ways) and installed in
editable mode with:

```python
git clone https://github.com/drug2ways/drug2ways.git
cd drug2ways
python -m pip install -e .
```

## Requirements
```python
click==7.1.1
tqdm==4.47.0
networkx>=2.1
pandas==1.0.3
networkx>=2.4
numpy
scipy
statsmodels
```
