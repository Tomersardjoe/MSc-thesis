#!/usr/bin/env python3

import argparse
import pandas as pd
import sys
import os
import re

def get_args():
    parser = argparse.ArgumentParser(description="Summarize real pangenome results")
    parser.add_argument("-i", "--infile", required=True,
                        help="Path to input CSV file (Coinfinder/Goldfinder/PanForest output)")
    parser.add_argument("-m", "--match_file", required=True,
                        help="Path to species_categories.csv or match_all.csv file")
    parser.add_argument("-o", "--outfile", required=True,
                        help="Path to summary TSV file (append mode)")
    return parser.parse_args()

def infer_run_id(infile):
    fname = os.path.basename(infile)
    m = re.search(r'(\d+)', fname)
    return m.group(1) if m else fname

def infer_tool(infile):
    path = infile.lower()
    if "coinfinder" in path:
        return "coinfinder"
    elif "goldfinder" in path:
        return "goldfinder"
    elif "panforest" in path:
        return "panforest"
    else:
        return "unknown"

def lookup_fluidity_openness(match_file, run_id):
    try:
        df = pd.read_csv(match_file)
    except Exception as e:
        sys.stderr.write(f"Warning: could not read {match_file}: {e}\n")
        return None, None

    try:
        run_id_int = int(run_id)
    except ValueError:
        return None, None

    row = df[df["species_taxid"] == run_id_int]
    if row.empty:
        return None, None

    try:
        fluidity = float(row["fluidity_calc"].iloc[0])
    except Exception:
        fluidity = None
    try:
        openness = float(row["openness"].iloc[0])
    except Exception:
        openness = None

    return fluidity, openness

def main():
    args = get_args()

    run_id = infer_run_id(args.infile)
    tool = infer_tool(args.infile).capitalize()

    # Count total significant pairs
    try:
        df = pd.read_csv(args.infile, dtype=str, low_memory=False)
    except pd.errors.EmptyDataError:
        total_pairs = 0
    else:
        if df.empty or "Level" not in df.columns:
            total_pairs = 0
        else:
            pairs = df[df["Level"] == "Pair"].copy()
            total_pairs = len(pairs)

    # Lookup fluidity/openness
    fluidity, openness = lookup_fluidity_openness(args.match_file, run_id)

    header = not os.path.exists(args.outfile) or os.path.getsize(args.outfile) == 0
    with open(args.outfile, "a") as out:
        if header:
            out.write("Method\tpangenome_id\ttotal_significant_pairs\tfluidity\topenness\n")

        out.write(
            f"{tool}\t{run_id}\t"
            f"{total_pairs}\t"
            f"{'' if fluidity is None else f'fluidity:.6f}'\t"
            f"{'' if openness is None else f'openness:.6f}'\n"
        )

if __name__ == "__main__":
    main()