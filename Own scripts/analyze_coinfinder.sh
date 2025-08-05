#!/bin/bash

folder="$1"

if [ ! -d "$folder" ]; then
    echo "❌ Error: '$folder' is not a valid directory."
    exit 1
fi

echo -e "\n📁 Analyzing Coinfinder outputs in: $folder\n"

for sim_dir in "$folder"/*; do
    [ -d "$sim_dir" ] || continue
    prefix=$(basename "$sim_dir")

    nodes="${sim_dir}/${prefix}_nodes.tsv"
    pairs="${sim_dir}/${prefix}_pairs.tsv"
    components="${sim_dir}/${prefix}_components.tsv"
    presence_file="simulation/${prefix}/gene_presence_absence.csv"

    echo "🔹 ${prefix}:"

    total_genes=0
    if [ -f "$presence_file" ]; then
        total_pangenome_genes=$(tail -n +2 "$presence_file" | wc -l)
        echo "   🌐 Total Pangenome Genes: $total_pangenome_genes"
    else
        echo "   ⚠️ presence file not found at $presence_file"
    fi

    if [ -f "$nodes" ]; then
        total_genes=$(tail -n +2 "$nodes" | cut -f1 | sort | uniq | wc -l)
        echo "   🧬 Total Genes (networked): $total_genes"
    else
        echo "   ⚠️ nodes.tsv not found"
    fi

    if [ -f "$pairs" ]; then
        sig_assoc=$(awk '$3 < 0.05' "$pairs" | wc -l)
        possible_pairs=$(( (total_genes * (total_genes - 1)) / 2 ))
        assoc_rate=$(awk -v a="$sig_assoc" -v p="$possible_pairs" 'BEGIN { if (p > 0) printf("%.2f", (a/p)*100); else print "0.00" }')

        echo "   🔢 Possible Gene Pairs: $possible_pairs"
        echo "   🧪 Significant Associations: $sig_assoc"
        echo "   📈 Association Rate: $assoc_rate%"
    else
        echo "   ⚠️ pairs.tsv not found"
    fi

    if [ -f "$components" ]; then
        module_count=$(cut -f2 "$components" | sort | uniq -c | wc -l)
        echo "   🧩 Module Count: $module_count"

        # Run Python script and collect metrics
        metrics_output=$(python3 network_metrics.py "$sim_dir")

        avg_module_size=$(echo "$metrics_output" | grep "Avg. Module Size" | awk -F ': ' '{print $2}')
        avg_degree=$(echo "$metrics_output" | grep "Avg. Degree" | awk -F ': ' '{print $2}')
        modularity=$(echo "$metrics_output" | grep -i "Modularity" | awk -F ': ' '{print $2}')

        if [ -n "$avg_module_size" ]; then
            echo "   📊 Avg. Genes per Module: $avg_module_size"
        else
            echo "   ⚠️ Failed to compute avg. module size"
        fi

        if [ -n "$avg_degree" ]; then
            echo "   🧮 Avg. Degree: $avg_degree"
        else
            echo "   ⚠️ Failed to compute avg. degree"
        fi

        if [ -n "$modularity" ]; then
            echo "   🧠 Network Modularity: $modularity"
        else
            echo "   ⚠️ Failed to compute network modularity"
        fi
    else
        echo "   ⚠️ components.tsv not found"
    fi

    echo ""
done
