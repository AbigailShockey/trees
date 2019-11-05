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
# The binary gene presence/absence file must have isolates as colnames and genes as rownames; if reverse, use -t option
# Bootstrapping of the neighbor-joining tree can be performed using -b option
# Bootstraps must be calculated with a separate script
# Based on code from: https://www.drawingfromdata.com/making-a-pdm-distance-matrix-with-pandas

# setup argparser to display help if no arguments
class ArgParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

# determine command line arguments and get path
parser = ArgParser(description='Calculate pairwise hamming distance matrix from gene prensece/absence matrix and create nj tree')
parser.add_argument('tsv', type=str, help="Binary gene presence and absence file (tab delimited)")
parser.add_argument('prefix', type=str, default="output", nargs='?', help="Output prefix, defaults to 'output'")
parser.add_argument('-t', action="store_true", required=False, help="Transpose tsv")
parser.add_argument('-b', type=int, metavar='bootstrapping', required=False, help="Perform boostrapping n times")

args = parser.parse_args()
def transpose_mat(df,tsvOut):
    tdf = df.transpose()
    tdf.to_csv(tsvOut, sep='\t', encoding='utf-8')
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

# read gene presence/absence tsv, treating first column as rownames and first row as column names
df = pd.read_csv(tsvfile, sep="\t", index_col=0, header=0)

# transpose gene presence/absence tsv and write to file
if args.t == True:
    tsvOut = os.path.basename(tsvfile).split(".")[0] + "_transposed.tsv"
    df = transpose_mat(df,tsvOut)

# hdm name and path
outMatrix = os.path.join(out,f"{prefix}_matrix.tsv")
# calculate hdm
calc_hamming_distance(df, outMatrix)

# nj tree name and path
outTree = os.path.join(out,f"{prefix}_nj_tree.newick")
# calculate nj tree from hdm
calc_nj_tree(outMatrix, outTree)

if args.b is not None:
    # create empty list of trees
    trees = dendropy.TreeList()
    
    # begin bootstrapping
    for i in range(1,(args.b + 1)):
        # randomly sample tsv columns with replacement
        rcols = np.random.choice(list(df.columns.values),len(list(df.columns.values)), replace=True)
        rdf = df[rcols]
        # randomly re-order tsv rows (see "jumble" option in boot.phylo function from R package ape)
        rdf = rdf.sample(frac=1, replace=False)
        # permuted hdm name and path
        outMatrix = os.path.join(out,f"matrix_permutation_{i}.tsv")
        # calculate hdm of permuted tsv
        calc_hamming_distance(rdf, outMatrix)
        # permuted nj tree name and path
        outTree = os.path.join(out,f"tree_permutation_{i}.newick")
        # calculate nj tree from hdm of permuted tsv
        rnj_tree = calc_nj_tree(outMatrix, outTree)
        # append permuted nj tree
        trees.append(rnj_tree)

    # write all permuted nj trees to file
    trees.write(
        path=os.path.join(out,"bootstrapped_nj_trees.newick"),
        schema="newick")
