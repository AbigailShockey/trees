#!/usr/bin/env python 3

import os,sys
import argparse
import pandas as pd
import scipy
from scipy import spatial
from scipy.spatial.distance import pdist,squareform
import numpy as np
import random
import dendropy

# This script takes a binary gene presence/absence file, calcualtes a pairwise hamming distance matrix and
# returns a neihbor-joining tree in newick format
# The binary gene presence/absence file must have isolates as colnames and genes as rownames; if not, use -t option
# Bootstrapping of the neighbor-joining can be performed using -b option
# Bootstraps must be calculated with a separate script
# Based on code from: https://www.drawingfromdata.com/making-a-pdm-distance-matrix-with-pandas

# Setup argparser to display help if no arguments
class ArgParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

# Determine command line arguments and get path
parser = ArgParser(description='Calculate pairwise hamming distance matrix from gene prensece/absence matrix and create nj tree')
parser.add_argument('tsv', type=str, help="Binary gene presence and absence file (tab delimited)")
parser.add_argument('prefix', type=str, default="output", nargs='?', help="Output prefix, defaults to 'output'")
parser.add_argument('-t', action="store_true", required=False, help="Transpose tsv")
parser.add_argument('-b', type=int, metavar='bootstrapping', required=False, help="Perform boostrapping n times")

args = parser.parse_args()

def transpose_mat(df):
    tdf = df.transpose()
    tdf.to_csv("tmat.mat", sep='\t', encoding='utf-8')
    return tdf

def calc_hamming_distance(df, outMatrix):
    hdm = pd.DataFrame(
        squareform(pdist(df, metric = 'hamming')),
        columns = df.index,
        index = df.index)
    hdm.to_csv(outMatrix, sep='\t', encoding='utf-8')

def calc_nj_tree(outMatrix, outTree):
    with open(outMatrix) as hdm:
        pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(
                hdm,
                is_first_row_column_names=True,
                is_first_column_row_names=True,
                is_allow_new_taxa=True,
                delimiter="\t")
    nj_tree = pdm.nj_tree()
    nj_tree.write(
        path=outTree,
        schema="newick")
    return nj_tree

# Set pandas display options (for printing)
pd.options.display.max_rows = 10
pd.options.display.max_columns = 5

tsvfile = os.path.abspath(args.tsv)
prefix = args.prefix
out = os.getcwd()

# Read gene presence/absence tsv
df = pd.read_csv(tsvfile, sep="\t")
# Set gene as rownames
df = df.set_index('Gene')
<<<<<<< Updated upstream
# print(df)

# transpose data frame so genes are columns and isolates are rows
df = df.transpose()

# write transposed dataframe to tsv
df.to_csv(f"{matrixFile_handle}_transposed.tsv", sep=',', encoding='utf-8')

# calculate hamming distance matrix
hdm = pd.DataFrame(
    squareform(pdist(df, metric = 'hamming')),
    columns = df.index,
    index = df.index)

# write hamming ditance matrix to file
hdm.to_csv(matrixFile, sep='\t', encoding='utf-8')

# read hamming distance matrix in as phylogenetic distance matrix
with open(matrixFile) as src:
    pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(
            src,
            is_first_row_column_names=True,
            is_first_column_row_names=True,
            is_allow_new_taxa=True,
            delimiter="\t")

# Calculate neighbor-joining tree from phylogenetic distance matrix and write to file
nj_tree = pdm.nj_tree()
nj_tree.write(
    path=treeFile,
    schema="newick")
=======

if args.t == True:
    df = transpose_mat(df)

outMatrix = os.path.join(out,f"{prefix}_matrix.tsv")
calc_hamming_distance(df, outMatrix)

outTree = os.path.join(out,f"{prefix}_nj_tree.newick")
calc_nj_tree(outMatrix, outTree)
>>>>>>> Stashed changes

if args.b is not None:
    i = 0
    while i < args.b:
        # Randomly sample columns with replacement
        rcols = np.random.choice(list(df.columns.values),len(list(df.columns.values)), replace=True)
        rdf = df[rcols]
        # Randomly re-order rows (see "jumble" option in boot.phylo function from R package ape)
        rdf = rdf.sample(frac=1, replace=False)
        # Calculate hamming distance matrix from permuted matrix
<<<<<<< Updated upstream
        rhdm = pd.DataFrame(
            squareform(pdist(rdf, metric = 'hamming')),
            columns = rdf.index,
            index = rdf.index)
        # Write permutation results to file
        rhdm.to_csv(f"{matrixFile_handle}_permutation_{i}.tsv", sep='\t', encoding='utf-8')
=======
        outMatrix = os.path.join(out,f"matrix_permutation_{i}.tsv")
        calc_hamming_distance(rdf, outMatrix)
>>>>>>> Stashed changes
        i = i + 1

    # Create empty list of trees
    trees = dendropy.TreeList()
    i = 0
    while i < args.b:
<<<<<<< Updated upstream
        # read permutation in as phylogenetic distance matrix
        with open(f"{matrixFile_handle}_permutation_{i}.tsv","r") as src:
            rpdm = dendropy.PhylogeneticDistanceMatrix.from_csv(
                    src,
                    is_first_row_column_names=True,
                    is_first_column_row_names=True,
                    is_allow_new_taxa=True,
                    delimiter="\t")
            # calculate neighbor-joining tree from permutation and write to file
            rnj_tree = rpdm.nj_tree()
            rnj_tree.write(
                path=f"{treeFile_handle}_permutation_{i}.newick",
                schema="newick")
        # append neighbor-joining tree
=======
        outTree = os.path.join(out,f"tree_permutation_{i}.newick")
        rnj_tree = calc_nj_tree(outMatrix, outTree)
        # append permuted neighbor-joining tree
>>>>>>> Stashed changes
        trees.append(rnj_tree)
        i = i + 1
    # write all permuted trees to file
    trees.write(
        path=os.path.join(out,"bootstrapped_trees.newick"),
        schema="newick")
