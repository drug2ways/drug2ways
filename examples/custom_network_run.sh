#!/usr/bin/env bash

# TODO: Set here your variables 'pkg_path' and 'output_path'
pkg_path=/Users/danieldomingo/PycharmProjects/drug2ways
output_path=/Users/danieldomingo/Downloads

lmin=3
lmax=9

if [ "$#" -eq 1 ]; then
	lmax=$1
elif [ "$#" -eq 2 ]; then
	lmin=$1
	lmax=$2
fi

echo "Running custom network for all lmax from ${lmin} to ${lmax} "

mkdir -p ${output_path}

for l in `seq $lmin $lmax`; do
	echo "Running custom network with lmax = ${l}."
	python3 -m drug2ways explore \
	--graph ${pkg_path}/data/networks/data/custom_network.tsv \
	--sources ${pkg_path}/data/validation/data/source_nodes_custom.tsv \
	--targets ${pkg_path}/data/validation/data/target_nodes_custom.tsv \
	--fmt tsv \
	--lmax ${l} \
	--output ${output_path} \
	--name custom \
	--log
done
