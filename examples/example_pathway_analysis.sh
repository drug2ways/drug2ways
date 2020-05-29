#!/usr/bin/env bash

# TODO: Set here your variables 'pkg_path' and 'output_path'
pkg_path=/Users/danieldomingo/PycharmProjects/drug2ways
output_path=/Users/danieldomingo/Downloads

lmin=8
lmax=8

if [ "$#" -eq 1 ]; then
	lmax=$1
elif [ "$#" -eq 2 ]; then
	lmin=$1
	lmax=$2
fi

echo "Running pathway analysis on the custom network for all lmax from ${lmin} to ${lmax} "

mkdir -p ${output_path}

for l in `seq $lmin $lmax`; do
	echo "Running custom network with lmax = ${l}."
	python3 -m drug2ways pathway-analysis \
	--graph ${pkg_path}/data/networks/data/openbiolink_network.tsv \
	--sources ${pkg_path}/data/validation/data/source_nodes_openbiolink.tsv \
	--targets ${pkg_path}/data/validation/data/target_nodes_openbiolink.tsv \
	--fmt tsv \
	--lmax ${l} \
	--output ${output_path}

done
