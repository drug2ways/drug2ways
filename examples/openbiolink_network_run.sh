#!/usr/bin/env bash

# TODO: Download the repo https://github.com/drug2ways/results
# TODO: Set here your variables 'pkg_path' and 'output_path'
# TODO: so pkg_path points to the cloned repo in your computer
pkg_path=/Users/danieldomingo/PycharmProjects/results
output_path=/Users/danieldomingo/Downloads

lmin=4
lmax=8

if [ "$#" -eq 1 ]; then
	lmax=$1
elif [ "$#" -eq 2 ]; then
	lmin=$1
	lmax=$2
fi

echo "Running openbiolink network for all lmax from ${lmin} to ${lmax} "


mkdir -p ${output_path}

for l in `seq $lmin $lmax`; do
	echo "Running openbiolink network with lmax = ${l}."
	python3 -m drug2ways explore \
	--graph ${pkg_path}/networks/data/openbiolink_network.tsv \
	--sources ${pkg_path}/validation/data/source_nodes_openbiolink.tsv \
	--targets ${pkg_path}/validation/data/target_nodes_openbiolink.tsv \
	--fmt tsv \
	--lmax ${l} \
	--output ${output_path} \
	--name openbiolink \
  --log
done
