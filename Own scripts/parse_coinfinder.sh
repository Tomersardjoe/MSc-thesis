#!/bin/bash

# Inputs
sim_dir="$1"
prefix="$2"

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
total_genes=$(tail -n +2 "$nodes" | cut -f1 | sort | uniq | wc -l)

# Association stats
sig_assoc=$(awk '$3 < 0.05' "$pairs" | wc -l)
possible_pairs=$(( (total_genes * (total_genes - 1)) / 2 ))
assoc_rate=$(awk -v a="$sig_assoc" -v p="$possible_pairs" 'BEGIN { if (p > 0) printf("%.2f", (a/p)*100); else print "0.00" }')

# Module stats
module_count=$(cut -f2 "$components" | sort | uniq -c | wc -l)

metrics_output=$(python3 network_metrics.py "$sim_dir")
avg_module_size=$(echo "$metrics_output" | grep -oP 'Avg\. Module Size.*: \K[\d.]+')
avg_degree=$(echo "$metrics_output" | grep -oP 'Avg\. Degree: \K[\d.]+')
modularity=$(echo "$metrics_output" | grep -oP 'Modularity: \K[\d.]+')

# Output line for unified CSV
echo "$total_pangenome_genes,$total_genes,$possible_pairs,$sig_assoc,$assoc_rate,$module_count,$avg_module_size,$avg_degree,$modularity"
