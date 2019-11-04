#!/usr/bin/env python 3

import os,sys
import argparse
import numpy as np
import scipy
from scipy import spatial
from scipy.spatial.distance import pdist,squareform
import pandas as pd
import dendropy
import random

# This script takes a binary gene presence/absence file, calcualtes a pairwise hamming distance matrix, returns a neihbor-joining tree in newick format and performs boostrapping
# Bootstraps must be calculated with a separate script
# Based on code from:  https://www.drawingfromdata.com/making-a-pdm-distance-matrix-with-pandas

# Setup argparser to display help if no arguments
class ArgParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

# Determine command line arguments and get path
parser = ArgParser(description='Calculate pairwise hamming distance matrix from gene prensece/absence matrix and create nj tree')
parser.add_argument('tsv', type=str, help="Binary gene presence and absence file (tab delimited)")
parser.add_argument('matrix', type=str, help="Output matrix name")
parser.add_argument('tree', type=str, help="Output file name")
parser.add_argument('-b', type=int, metavar='bootstrapping', help="Number of iterations for bootstrapping", required=False)


args = parser.parse_args()
tsvfile = os.path.abspath(args.tsv)
matrixFile = os.path.abspath(args.matrix)
treeFile = os.path.abspath(args.tree)

matrixFile_handle = args.matrix.split(".")[0]
treeFile_handle = args.tree.split(".")[0]

# Set pandas display options (for printing)
pd.options.display.max_rows = 10
pd.options.display.max_columns = 6

# Read gene presence/absence tsv, set gene as rownames
df = pd.read_csv(tsvfile, sep="\t")
df = df.set_index('Gene')
# print(df)

# transpose data frame so genes are columns and isolates are rows
df = df.transpose()

# write transposed dataframe to tsv
df.to_csv(f"{matrixFile_handle}_transposed.tsv", sep=',', encoding='utf-8')

# calculate hamming distance matrix
pdm = pd.DataFrame(
    squareform(pdist(df, metric = 'hamming')),
    columns = df.index,
    index = df.index)

# write hamming ditance matrix to file
pdm.to_csv(matrixFile, sep='\t', encoding='utf-8')

# read hamming distance matrix in as phylogenetic distance matrix
with open(matrixFile) as src:
    pdm_phy = dendropy.PhylogeneticDistanceMatrix.from_csv(
            src,
            is_first_row_column_names=True,
            is_first_column_row_names=True,
            is_allow_new_taxa=True,
            delimiter="\t")

# Calculate neighbor-joining tree from phylogenetic distance matrix and write to file
nj_tree = pdm_phy.nj_tree()
nj_tree.write(
    path=treeFile,
    schema="newick")

if args.b is not None:
    i = 0

    # Begin bootstrapping
    while i < args.b:
        # Randomly sample columns with replacement
        rcols = np.random.choice(list(df.columns.values),len(list(df.columns.values)), replace=True)
        rdf = df[rcols]
        # Randomly re-order rows (see "jumble" option in boot.phylo function from R package ape)
        rdf = rdf.sample(frac=1, replace=False)
        # Calculate hamming distance matrix from permuted matrix
        rpdm = pd.DataFrame(
            squareform(pdist(rdf, metric = 'hamming')),
            columns = rdf.index,
            index = rdf.index)
        # Write permutation results to file
        rpdm.to_csv(f"{matrixFile_handle}_permutation_{i}.tsv", sep='\t', encoding='utf-8')
        i = i + 1

    i = 0

    # Create empty list of trees
    trees = dendropy.TreeList()

    while i < args.b:
        # read permutation in as phylogenetic distance matrix
        with open(f"{matrixFile_handle}_permutation_{i}.tsv","r") as src:
            rpdm_phy = dendropy.PhylogeneticDistanceMatrix.from_csv(
                    src,
                    is_first_row_column_names=True,
                    is_first_column_row_names=True,
                    is_allow_new_taxa=True,
                    delimiter="\t")
            # calculate neighbor-joining tree from permutation and write to file
            rnj_tree = rpdm_phy.nj_tree()
            rnj_tree.write(
                path=f"{treeFile_handle}_permutation_{i}.newick",
                schema="newick")
        # append neighbor-joining tree
        trees.append(rnj_tree)
        i = i + 1

    # write all permuted trees to file
    trees.write(
        path="bootstrapped_trees.newick",
        schema="newick")
