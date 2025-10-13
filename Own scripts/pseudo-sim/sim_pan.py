import argparse
import pandas as pd
import numpy as np
import random

def main():
    parser = argparse.ArgumentParser(
        description="Insert shuffled duplicates into random positions without changing original gene order"
    )
    parser.add_argument("nodes_tsv", help="Path to nodes_all.tsv file")
    parser.add_argument("gpa_csv", help="Path to gene_presence_absence.csv file")
    parser.add_argument(
        "-o", "--output", default="gpa_with_duplicates.csv",
        help="Output CSV file (default: gpa_with_duplicates.csv)"
    )
    args = parser.parse_args()

    # Load nodes_all.tsv
    nodes_df = pd.read_csv(args.nodes_tsv, sep="\t", encoding="utf-8")
    d_values = nodes_df["Result"]

    # Compute 5th and 95th percentiles
    p5, p95 = d_values.quantile(0.05), d_values.quantile(0.95)

    # Filter and select 10% of filtered genes
    filtered_df = nodes_df[(d_values >= p5) & (d_values <= p95)].copy()
    filtered_df = filtered_df.sort_values(by="Result").reset_index(drop=True)
    num_select = max(1, int(len(filtered_df) * 0.10))
    indices = np.linspace(0, len(filtered_df) - 1, num_select, dtype=int)
    selected_genes = set(filtered_df.iloc[indices]["ID"].astype(str))

    # Load gpa file
    gpa_df = pd.read_csv(args.gpa_csv, index_col=0, encoding="utf-8")

    # Find matches and create duplicates
    matches = [gene for gene in gpa_df.index if gene in selected_genes]
    dup_rows = gpa_df.loc[matches].copy()
    dup_rows.index = [f"{idx}_dup" for idx in dup_rows.index]

    # Shuffle duplicates
    dup_rows = dup_rows.sample(frac=1).reset_index()

    # Pre-assign random insertion positions for each duplicate
    n_dups = len(dup_rows)
    insert_positions = np.random.choice(
        range(len(gpa_df) + n_dups), size=n_dups, replace=False
    )
    insert_positions.sort()

    # Reset dup_rows for iteration
    dup_rows = dup_rows.reset_index(drop=True)

    result_rows = []
    dup_idx = 0
    for i, (gene, row) in enumerate(gpa_df.reset_index().iterrows()):
        # Add the original gene row
        result_rows.append(row)

        # Insert any duplicates assigned to this position
        while dup_idx < n_dups and insert_positions[dup_idx] == len(result_rows) - 1:
            result_rows.append(dup_rows.iloc[dup_idx])
            dup_idx += 1

    # If any duplicates were assigned after the last row
    while dup_idx < n_dups:
        result_rows.append(dup_rows.iloc[dup_idx])
        dup_idx += 1

    # Build final dataframe
    combined_df = pd.DataFrame(result_rows).set_index("index")
    combined_df.to_csv(args.output)

    print(f"Original genes: {len(gpa_df)}")
    print(f"Duplicated genes inserted: {len(matches)}")
    print(f"Final rows: {combined_df.shape[0]}")
    print(f"Saved updated gpa to: {args.output}")

if __name__ == "__main__":
    main()
