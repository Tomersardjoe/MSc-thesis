#!/usr/bin/env python3

import argparse
import pandas as pd

def filter_matrix(matrix_path, dstat_path, output_path):
    # Load presence/absence matrix
    matrix = pd.read_csv(matrix_path, sep=",", index_col=0)

    # Load D statistics table
    dstat = pd.read_csv(dstat_path, sep="\t")

    # Ensure required columns exist
    if not {"ID", "Result"}.issubset(dstat.columns):
        raise ValueError("D-stat file must contain 'ID' and 'Result' columns")

    # Keep only genes with D > 0
    keep_genes = dstat.loc[dstat["Result"] > 0, "ID"]

    # Filter matrix rows
    filtered = matrix.loc[matrix.index.intersection(keep_genes)]

    # Save filtered matrix
    filtered.to_csv(output_path, sep=",")

    print(f"D > 0 filtered matrix saved to {output_path} with {filtered.shape[0]} genes retained.")

def main():
    parser = argparse.ArgumentParser(description="Filter gene presence/absence matrix by D statistic > 0")
    parser.add_argument("matrix", help="Path to gene presence/absence matrix (SV, rows=genes, cols=genomes)")
    parser.add_argument("dstat", help="Path to coincident_nodes_all.tsv (with ID and Result columns)")
    parser.add_argument("-o", "--output", default="filtered_matrix.csv", help="Output path for filtered matrix")
    args = parser.parse_args()

    filter_matrix(args.matrix, args.dstat, args.output)

if __name__ == "__main__":
    main()
