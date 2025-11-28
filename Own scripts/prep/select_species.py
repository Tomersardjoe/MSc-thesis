#!/usr/bin/env python3

import pandas as pd
import os
import shutil
import argparse

def subset_files(species_list, gpa_dir, tree_dir, outdir):
    gpa_out = os.path.join(outdir, "gpa_matches")
    tree_out = os.path.join(outdir, "tree_matches")
    os.makedirs(gpa_out, exist_ok=True)
    os.makedirs(tree_out, exist_ok=True)

    # GPA: taxid = prefix before first underscore
    for f in os.listdir(gpa_dir):
        if f.endswith(".csv") or f.endswith(".tab"):
            taxid = f.split("_")[0]
            if taxid in species_list:
                shutil.copy(os.path.join(gpa_dir, f), os.path.join(gpa_out, f))

    # Tree: taxid = suffix after first underscore, before .nwk
    for f in os.listdir(tree_dir):
        if f.endswith(".nwk"):
            parts = os.path.splitext(f)[0].split("_")
            if len(parts) > 1:
                taxid = parts[1]
                if taxid in species_list:
                    shutil.copy(os.path.join(tree_dir, f), os.path.join(tree_out, f))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Subset matched_all.csv and associated GPA/tree files")
    parser.add_argument("--taxids", required=True, help="CSV file with species_taxid column (User created file, create this file yourself first! e.g., real_pangenomes/species.csv")
    parser.add_argument("--matched", required=True, help="Path to matched_all.csv")
    parser.add_argument("--gpa", required=True, help="Path to gpa_matches_all directory")
    parser.add_argument("--tree", required=True, help="Path to tree_matches_all directory")
    parser.add_argument("--outdir", default=".", help="Output directory")
    args = parser.parse_args()

    taxid_df = pd.read_csv(args.taxids)
    species_list = taxid_df['species_taxid'].astype(str).tolist()

    matched_df = pd.read_csv(args.matched)
    matched_df['species_taxid'] = matched_df['species_taxid'].astype(str)
    subset_df = matched_df[matched_df['species_taxid'].isin(species_list)]

    os.makedirs(args.outdir, exist_ok=True)
    subset_df.to_csv(os.path.join(args.outdir, "species_categories.csv"), index=False)

    subset_files(species_list, args.gpa, args.tree, args.outdir)

    print(f"Subset written to {args.outdir}species_categories.csv")
    print(f"GPA files copied to {args.outdir}gpa_matches")
    print(f"Tree files copied to {args.outdir}tree_matches")
