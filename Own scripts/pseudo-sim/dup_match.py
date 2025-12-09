#!/usr/bin/env python3

import argparse
import pandas as pd
import sys
import os
import re

def get_args():
    parser = argparse.ArgumentParser(description="Summarize duplicate recovery")
    parser.add_argument("-i", "--infile", required=True,
                        help="Path to input CSV file (Coinfinder/Goldfinder/PanForest output)")
    parser.add_argument("-d", "--dupfile", required=True,
                        help="Path to duplicates{run_id}_REDUCED.csv file")
    parser.add_argument("-m", "--match_file", required=True,
                        help="Path to species_categories.csv or match_all.csv file")
    parser.add_argument("-o", "--outfile", required=True,
                        help="Path to summary TSV file (append mode)")
    return parser.parse_args()

def count_duplicates(dupfile):
    try:
        df = pd.read_csv(dupfile, usecols=[0], dtype=str, low_memory=False)
    except pd.errors.EmptyDataError:
        return 0
    if df.empty:
        return 0
    return df.iloc[:,0].str.endswith("_dup").sum()

def infer_run_id(infile):
    fname = os.path.basename(infile)
    m = re.search(r'(\d+)', fname)
    return m.group(1) if m else fname

def infer_tool(infile):
    path = infile.lower()
    if "coinfinder" in path:
        return "Coinfinder"
    elif "goldfinder" in path:
        return "Goldfinder"
    elif "panforest" in path:
        return "PanForest"
    else:
        return "Unknown"

def lookup_fluidity_openness(match_file, run_id):
    try:
        df = pd.read_csv(match_file)
    except Exception as e:
        sys.stderr.write(f"Warning: could not read {match_file}: {e}\n")
        return None, None, None

    try:
        run_id_int = int(run_id)
    except ValueError:
        return None, None, None

    row = df[df["species_taxid"] == run_id_int]
    if row.empty:
        return None, None, None

    try:
        fluidity = float(row["fluidity_calc"].iloc[0])
    except Exception:
        fluidity = None
    try:
        openness = float(row["openness"].iloc[0])
    except Exception:
        openness = None
    try:
        category = str(row["category"].iloc[0])
    except Exception:
        category = None

    return fluidity, openness, category

def main():
    args = get_args()

    run_id = infer_run_id(args.infile)
    tool = infer_tool(args.infile)

    total_dups = count_duplicates(args.dupfile)

    try:
        df = pd.read_csv(args.infile, dtype=str, low_memory=False)
    except pd.errors.EmptyDataError:
        found = 0
        total_pairs = 0
    else:
        if df.empty or "Level" not in df.columns:
            found = 0
            total_pairs = 0
        else:
            pairs = df[df["Level"] == "Pair"].copy()
            total_pairs = len(pairs)
            if pairs.empty or not {"Gene_1","Gene_2"}.issubset(pairs.columns):
                found = 0
            else:
                pairs["Gene_1_base"] = pairs["Gene_1"].astype(str).str.replace("_dup$", "", regex=True)
                pairs["Gene_2_base"] = pairs["Gene_2"].astype(str).str.replace("_dup$", "", regex=True)
                dup_pairs = pairs[
                    (pairs["Gene_1_base"] == pairs["Gene_2_base"]) &
                    (pairs["Gene_1"] != pairs["Gene_2"])
                ]
                found = len(dup_pairs)

    pct_total_dups = (found / total_dups * 100) if total_dups > 0 else 0.0
    pct_total_pairs = (found / total_pairs * 100) if total_pairs > 0 else 0.0
    
    TP = found
    FN = max(total_dups - found, 0)
    FP = max(total_pairs - found, 0)
    
    precision = TP / (TP + FP) if (TP + FP) > 0 else 0.0
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

    fluidity, openness, category = lookup_fluidity_openness(args.match_file, run_id)
    
    header = not os.path.exists(args.outfile) or os.path.getsize(args.outfile) == 0
    with open(args.outfile, "a") as out:
        if header:
            out.write(
                "Method\tpangenome_id\t"
                "total_significant_pairs\t"
                "duplicate_genes_found\t"
                "total_duplicate_genes\t"
                "duplicate_found_pct\t"
                "dup_as_pct_of_pairs\t"
                "precision\trecall\tf1\t"
                "fluidity\topenness\tcategory\n"
            )
            
        out.write(
            f"{tool}\t{run_id}\t"
            f"{total_pairs}\t"
            f"{TP}\t"
            f"{total_dups}\t"
            f"{pct_total_dups:.2f}\t"
            f"{pct_total_pairs:.2f}\t"
            f"{precision:.3f}\t{recall:.3f}\t{f1:.3f}\t"
            f"{'' if fluidity is None else f'{fluidity:.6f}'}\t"
            f"{'' if openness is None else f'{openness:.6f}'}\t"
            f"{'' if category is None else category}\n"
        )

if __name__ == "__main__":
    main()
