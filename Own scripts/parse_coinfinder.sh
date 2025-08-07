#!/bin/bash

# Inputs
sim_dir="$1"
prefix="$2"

dataset="$prefix"

nodes="${sim_dir}/${prefix}_nodes.tsv"
pairs="${sim_dir}/${prefix}_pairs.tsv"
components="${sim_dir}/${prefix}_components.tsv"
presence_file="simulation/${prefix}/gene_presence_absence.csv"

# Checks
if [ ! -f "$nodes" ] || [ ! -f "$pairs" ] || [ ! -f "$components" ] || [ ! -f "$presence_file" ]; then
    echo "ERROR: Missing required files in $sim_dir" >&2
    exit 1
fi

# Total pangenome genes
total_pangenome_genes=$(tail -n +2 "$presence_file" | wc -l)

# Networked genes
network_nodes=$(tail -n +2 "$nodes" | cut -f1 | sort | uniq | wc -l)

# Association stats
gene_associations=$(awk '$3 < 0.05' "$pairs" | wc -l)
possible_pairs=$(( (network_nodes * (network_nodes - 1)) / 2 ))
association_rate=$(awk -v a="$gene_associations" -v p="$possible_pairs" 'BEGIN { if (p > 0) printf("%.2f", (a/p)*100); else print "0.00" }')

# Module stats
module_count=$(cut -f2 "$components" | sort | uniq -c | wc -l)

metrics_output=$(python3 network_metrics.py "$sim_dir")
avg_genes_per_module=$(echo "$metrics_output" | grep -oP 'Avg\. Module Size.*: \K[\d.]+')
avg_degree=$(echo "$metrics_output" | grep -oP 'Avg\. Degree: \K[\d.]+')
modularity=$(echo "$metrics_output" | grep -oP 'Modularity: \K[\d.]+')

# Output line for unified CSV
echo "$dataset,$total_pangenome_genes,$network_nodes,$possible_pairs,$gene_associations,$association_rate,$module_count,$avg_genes_per_module,$avg_degree,$modularity"
