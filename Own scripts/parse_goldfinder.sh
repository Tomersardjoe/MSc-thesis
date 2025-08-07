#!/bin/bash

# Inputs
sim_dir="$1"
prefix="$2"

dataset="$prefix"

presence_file="simulation/${prefix}/gene_presence_absence.csv"
pairs="${sim_dir}/simultaneous_association_significant_pairs.csv"
clusters="${sim_dir}/association_clusters.txt"
network="${sim_dir}/cytoscape_input.csv"

# Checks
if [ ! -f "$presence_file" ] || [ ! -f "$pairs" ] || [ ! -f "$clusters" ] || [ ! -f "$network" ]; then
    echo "ERROR: Missing required files in $sim_dir" >&2
    exit 1
fi

# Total pangenome genes
total_pangenome_genes=$(tail -n +2 "$presence_file" | wc -l)

# Network stats
network_edges=$(grep -v '^\s*$' "$network" | tail -n +2 | wc -l)
network_nodes=$(tail -n +2 "$network" | cut -d',' -f1,2 | tr ',' '\n' | sort | uniq | wc -l)

# Association stats
gene_associations=$(tail -n +2 "$pairs" | wc -l)
possible_pairs=$(( (network_nodes * (network_nodes - 1)) / 2 ))
association_rate=$(awk -v a="$gene_associations" -v p="$possible_pairs" 'BEGIN { if (p > 0) printf("%.2f", (a/p)*100); else print "0.00" }')

# Cluster stats
module_count=$(grep -c '^>' "$clusters")
avg_genes_per_module=$(
  awk -F, '
    /^>/ {
      clusters++
      genes += $2
    }
    END {
      if (clusters > 0) {
        printf "%.2f", genes / clusters
      } else {
        print "NA"
      }
    }
  ' "$clusters"
)

# Network metrics via Python
metrics_output=$(python3 network_metrics.py "$sim_dir")
avg_degree=$(echo "$metrics_output" | grep -oP 'Avg\. Degree: \K[\d.]+')
modularity=$(echo "$metrics_output" | grep -oP 'Modularity: \K[\d.]+')

# Output line for unified CSV
echo "$dataset,$total_pangenome_genes,$network_nodes,$possible_pairs,$gene_associations,$association_rate,$module_count,$avg_genes_per_module,$avg_degree,$modularity"
