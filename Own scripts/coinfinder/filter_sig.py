import pandas as pd
from pathlib import Path
import shutil
import sys

def filter_edges_nodes_components(input_dir, fdr_cutoff=0.05, overwrite=False):
    input_dir = Path(input_dir)
    pairs_files = list(input_dir.glob("*_pairs.tsv"))

    if not pairs_files:
        print("No _pairs.tsv files found.")
        return

    for pairs_file in pairs_files:
        prefix = pairs_file.stem.replace("_pairs", "")
        edges_file = input_dir / f"{prefix}_edges.tsv"
        nodes_file = input_dir / f"{prefix}_nodes.tsv"
        comps_file = input_dir / f"{prefix}_components.tsv"

        if not edges_file.exists():
            print(f"Skipping {prefix}: no matching _edges.tsv found.")
            continue

        print(f"\nProcessing: {prefix} (cutoff {fdr_cutoff})")

        # Load pairs and edges
        pairs_df = pd.read_csv(pairs_file, sep="\t")
        edges_df = pd.read_csv(edges_file, sep="\t")

        # Build sig_pairs, sig_keys, and genes_to_keep
        sig_pairs = pairs_df[pairs_df["FDR_BH"] < fdr_cutoff][["Source", "Target", "FDR_BH"]]
        sig_pairs["key"] = sig_pairs.apply(lambda row: tuple(sorted([row["Source"], row["Target"]])), axis=1)
        sig_keys = set(sig_pairs["key"])
        genes_to_keep = set(sig_pairs["Source"]).union(sig_pairs["Target"])

        # Add keys to edges for comparison
        edges_df["key"] = edges_df.apply(lambda row: tuple(sorted([row["Source"], row["Target"]])), axis=1)
        edge_keys = set(edges_df["key"])
        edges_already_filtered = edge_keys.issubset(sig_keys) and len(edge_keys) == len(sig_keys)

        # Nodes check
        nodes_already_filtered = False
        if nodes_file.exists():
            nodes_df = pd.read_csv(nodes_file, sep="\t")
            node_ids = set(nodes_df[nodes_df.columns[0]])
            nodes_already_filtered = node_ids.issubset(genes_to_keep) and len(node_ids) == len(genes_to_keep)

        # Components check
        comps_already_filtered = False
        if comps_file.exists():
            comps_df = pd.read_csv(comps_file, sep="\t", header=None, names=["Component", "Genes"])
            comp_genes = set(g for sublist in comps_df["Genes"].str.split(',') for g in sublist)
            comps_already_filtered = comp_genes.issubset(genes_to_keep) and len(comp_genes) == len(genes_to_keep)

        # Master skip
        if edges_already_filtered and nodes_already_filtered and comps_already_filtered and not overwrite:
            print(f"Skipping {prefix}: edges, nodes and components already match significant sets.")
            continue

        # Merge to append FDR_BH to edges and overwrite
        merged_edges = edges_df.merge(sig_pairs[["key", "FDR_BH"]], on="key", how="inner").drop(columns="key")
        merged_edges.to_csv(edges_file, sep="\t", index=False)
        print(f"Overwritten edges: {edges_file}")

        # Filter nodes
        if nodes_file.exists():
            # Make a copy of the ALL nodes
            copy_file = nodes_file.with_name(f"{nodes_file.stem}_all{nodes_file.suffix}")
            shutil.copy2(nodes_file, copy_file)
            print(f"File with all nodes saved as: {copy_file}")
        
            filtered_nodes = nodes_df[nodes_df[nodes_df.columns[0]].isin(genes_to_keep)]
            filtered_nodes.to_csv(nodes_file, sep="\t", index=False)
            print(f"Overwritten nodes: {nodes_file}")
        else:
            print(f"No matching _nodes.tsv found for {prefix}")

        # Filter components
        if comps_file.exists():
            filtered_rows = []
            for _, row in comps_df.iterrows():
                genes = row["Genes"].split(',')
                kept_genes = [g for g in genes if g in genes_to_keep]
                if kept_genes:
                    filtered_rows.append({"Genes": ','.join(kept_genes)})

            filtered_comps = pd.DataFrame(filtered_rows)
            filtered_comps.insert(0, "Component", range(1, len(filtered_comps) + 1))
            filtered_comps.to_csv(comps_file, sep="\t", index=False, header=False)
            print(f"Overwritten components: {comps_file}")
        else:
            print(f"No matching _components.tsv found for {prefix}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python filter_sig.py <input_directory> [fdr_cutoff] [--overwrite]")
        sys.exit(1)

    input_directory = sys.argv[1]

    # Arg parse for significance threshold and overwrite flag
    cutoff = 0.05
    overwrite_flag = False
    for arg in sys.argv[2:]:
        if arg.startswith("--"):
            if arg == "--overwrite":
                overwrite_flag = True
        else:
            cutoff = float(arg)

    filter_edges_nodes_components(input_directory, cutoff, overwrite_flag)
