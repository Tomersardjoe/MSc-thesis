import os
import glob
import argparse
import networkx as nx
from networkx.algorithms.community.quality import modularity

def load_graph_from_pairs(pairs_file, pval_threshold=0.05):
    G = nx.Graph()
    try:
        with open(pairs_file, 'r') as f:
            next(f)  # Skip header
            print(f"✅ Reading edges from {pairs_file}")
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    print(f"⚠️ Skipping malformed line: {line.strip()}")
                    continue
                gene_a, gene_b = parts[0], parts[1]
                try:
                    pval = float(parts[2])
                except ValueError:
                    print(f"⚠️ Invalid p-value on line: {line.strip()}")
                    continue
                if pval < pval_threshold:
                    G.add_edge(gene_a, gene_b)
        print(f"🧬 Graph has {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")
    except Exception as e:
        print(f"⚠️ Error reading {pairs_file}: {e}")
    return G

def get_communities(components_file):
    communities = []
    try:
        with open(components_file, 'r') as f:
            print(f"✅ Reading communities from {components_file}")
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    genes = parts[1].split(',')
                    communities.append(genes)
                else:
                    print(f"⚠️ Skipping line without genes: {line.strip()}")
        print(f"📦 Found {len(communities)} communities")
    except Exception as e:
        print(f"⚠️ Error reading {components_file}: {e}")
    return communities

def calculate_avg_module_size(components_file):
    sizes = []
    try:
        with open(components_file, 'r') as f:
            print(f"✅ Calculating module sizes from {components_file}")
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    genes = parts[1].split(',')
                    sizes.append(len(genes))
                else:
                    print(f"⚠️ Skipping line without gene list: {line.strip()}")
        print(f"📊 Module sizes: {sizes}")
    except Exception as e:
        print(f"⚠️ Error reading {components_file}: {e}")
    return round(sum(sizes) / len(sizes), 2) if sizes else 0

def main():
    parser = argparse.ArgumentParser(description="Compute module and network metrics from Coinfinder outputs.")
    parser.add_argument("folder", help="Folder containing *_components.tsv and *_pairs.tsv files")
    args = parser.parse_args()

    print(f"📁 Searching in folder: {args.folder}")
    components_files = glob.glob(os.path.join(args.folder, "*_components.tsv"))

    if not components_files:
        print("❌ No *_components.tsv files found in specified folder.")
        return

    for comp_file in components_files:
        prefix = os.path.basename(comp_file).replace("_components.tsv", "")
        pairs_file = os.path.join(args.folder, f"{prefix}_pairs.tsv")

        print(f"\n🔹 Processing prefix: {prefix}")

        # Avg module size
        avg_size = calculate_avg_module_size(comp_file)
        print(f"📊 Avg. Module Size for {prefix}: {avg_size}")

        # Network metrics
        if os.path.exists(pairs_file):
            G = load_graph_from_pairs(pairs_file)
            if G.number_of_nodes() > 0:
                avg_degree = round(sum(dict(G.degree()).values()) / G.number_of_nodes(), 2)
                print(f"🧮 Avg. Degree: {avg_degree}")

                communities = get_communities(comp_file)
                try:
                    mod_score = round(modularity(G, communities), 4)
                    print(f"🧩 Modularity: {mod_score}")
                except Exception as e:
                    mod_score = f"⚠️ modularity error: {e}"
                    print(mod_score)
            else:
                print("⚠️ Graph is empty. Cannot compute degree or modularity.")
        else:
            print(f"⚠️ Missing pairs.tsv file: {pairs_file}")

if __name__ == "__main__":
    main()
