#!/usr/bin/env bash

# Process matrices
for type in closed moderate open; do
  for i in {1..10}; do
    echo "Starting $type rep${i}..."
    python3 process_matrix.py \
      -i ../simulation/output/${type}/grid_1/rep${i}/gene_presence_absence_all.csv \
      -o "collapsed_matrix.csv" \
      -d ${type}/rep${i}/
    echo "Finished $type rep${i}"
  done
done

# Run PanForest
for type in closed moderate open; do
  for i in {1..10}; do
    echo "Starting PanForest for $type rep${i}..."
    python3 PanForest.py \
      -n 500 \
      -d 8 \
      -m ${type}/rep${i}/collapsed_matrix.csv \
      -pres 1 \
      -abs 1 \
      -o ${type}/rep${i}/
    echo "Finished PanForest for $type rep${i}"
  done
done