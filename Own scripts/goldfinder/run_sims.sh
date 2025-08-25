#!/usr/bin/env bash

for type in closed moderate open; do
  for i in {1..10}; do
    echo "Starting $type rep${i}..."
    python3 goldfinder.py \
      -i ../../simulation/output/${type}/grid_1/rep${i}/gene_presence_absence_all.csv \
      -t ../../simulation/output/${type}/grid_1/rep${i}/tree.nwk \
      -o ${type}/rep${i}/ \
      -c both
    echo "Finished $type rep${i}"
  done
done
