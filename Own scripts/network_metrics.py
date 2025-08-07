#!/usr/bin/env python3
import os
import sys
import csv
import networkx as nx

# Resolve where this script lives (and the project root)
SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)

def load_goldfinder_network(file_path):
    G = nx.Graph()
    with open(file_path, newline='') as csvfile:
        reader = csv.reader(csvfile)
        next(reader, None)  # skip header
        for row in reader:
            if len(row) >= 2:
                G.add_edge(row[0], row[1])
    return G

def load_coinfinder_network(dir_path):
    for fname in os.listdir(dir_path):
        if fname.endswith("_pairs.tsv"):
            file_path = os.path.join(dir_path, fname)
            G = nx.Graph()
            with open(file_path, newline='') as tsvfile:
                reader = csv.reader(tsvfile, delimiter='\t')
                next(reader, None)  # skip header
                for row in reader:
                    if len(row) >= 2:
                        G.add_edge(row[0], row[1])
            return G
    raise FileNotFoundError("No *_pairs.tsv file found in Coinfinder directory.")

def calculate_avg_degree(G):
    degrees = [deg for _, deg in G.degree()]
    return round(sum(degrees) / len(degrees), 2) if degrees else 0.0

def calculate_modularity(G):
    try:
        import community  # python-louvain
        partition = community.best_partition(G)
        return round(community.modularity(partition, G), 2)
    except ImportError:
        return "NA"

def compute_goldfinder_module_size(file_path):
    cluster_sizes = []
    current_size = 0
    try:
        with open(file_path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if current_size > 0:
                        cluster_sizes.append(current_size)
                    current_size = 0
                else:
                    current_size += 1
            if current_size > 0:
                cluster_sizes.append(current_size)
        if cluster_sizes:
            return round(sum(cluster_sizes) / len(cluster_sizes), 2)
    except Exception as e:
        print(f"Error reading Goldfinder clusters: {e}", file=sys.stderr)
    return "NA"

def compute_coinfinder_module_size(file_path):
    module_sizes = []
    try:
        with open(file_path, newline='') as tsvfile:
            reader = csv.reader(tsvfile, delimiter='\t')
            for row in reader:
                if len(row) >= 2:
                    genes = [g.strip() for g in row[1].split(',') if g.strip()]
                    module_sizes.append(len(genes))
        if module_sizes:
            return round(sum(module_sizes) / len(module_sizes), 2)
    except Exception as e:
        print(f"Error reading Coinfinder components: {e}", file=sys.stderr)
    return "NA"

def main():
    if len(sys.argv) != 2:
        print(f"Usage: {os.path.basename(__file__)} <results_directory>", file=sys.stderr)
        sys.exit(1)

    input_dir = sys.argv[1]
    gold_net = os.path.join(input_dir, "cytoscape_input.csv")
    gold_clust = os.path.join(input_dir, "association_clusters.txt")

    try:
        if os.path.isfile(gold_net) and os.path.isfile(gold_clust):
            G = load_goldfinder_network(gold_net)
            avg_module_size = compute_goldfinder_module_size(gold_clust)
        else:
            G = load_coinfinder_network(input_dir)
            comp_file = next(
                (os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith("_components.tsv")),
                None
            )
            avg_module_size = compute_coinfinder_module_size(comp_file) if comp_file else "NA"

        avg_degree = calculate_avg_degree(G)
        modularity = calculate_modularity(G)

        print(f"Avg. Degree: {avg_degree}")
        print(f"Modularity: {modularity}")
        print(f"Avg. Module Size: {avg_module_size}")

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        print("Avg. Degree: NA")
        print("Modularity: NA")
        print("Avg. Module Size: NA")

if __name__ == "__main__":
    main()
