#!/usr/bin/env python3

import argparse
import pandas as pd

def get_args():
    parser = argparse.ArgumentParser(description="Find gene pairs where one is the _dup of the other")
    parser.add_argument("-i", "--infile", required=True, help="Path to input CSV file")
    return parser.parse_args()

def main():
    args = get_args()
    df = pd.read_csv(args.infile)

    # Keep only rows where Level == "Pair"
    pairs = df[df["Level"] == "Pair"].copy()

    # Normalize gene IDs by stripping "_dup"
    pairs["Gene_1_base"] = pairs["Gene_1"].str.replace("_dup$", "", regex=True)
    pairs["Gene_2_base"] = pairs["Gene_2"].str.replace("_dup$", "", regex=True)

    # Keep rows where the base IDs match but the raw IDs differ
    dup_pairs = pairs[
        (pairs["Gene_1_base"] == pairs["Gene_2_base"]) &
        (pairs["Gene_1"] != pairs["Gene_2"])
    ]

    print("Found the following _dup pairs:")
    print(dup_pairs[["Gene_1", "Gene_2"]])

if __name__ == "__main__":
    main()
