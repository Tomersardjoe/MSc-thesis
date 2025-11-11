import csv
import glob
import os
import sys

def aggregate_counts(directory):
    results = []

    # Loop over all CSV files
    for filepath in glob.glob(os.path.join(directory, "*.csv")):
        with open(filepath, newline='') as f:
            reader = list(csv.reader(f))

        num_rows = len(reader)
        num_cols = len(reader[0]) if reader else 0

        num_genes = num_rows - 1
        num_taxa = num_cols - 1

        results.append({
            "file": os.path.basename(filepath),
            "num_taxa": num_taxa,
            "num_genes": num_genes
        })

    # Write aggregated results
    parent_dir = os.path.abspath(os.path.join(directory, os.pardir))
    outpath = os.path.join(parent_dir, "aggregated_counts.csv")
    with open(outpath, "w", newline='') as f:
        writer = csv.DictWriter(f, fieldnames=["file", "num_taxa", "num_genes"])
        writer.writeheader()
        writer.writerows(results)

    print(f"Aggregated results written to {outpath}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python aggregate_counts.py <directory>")
        sys.exit(1)
    aggregate_counts(sys.argv[1])
