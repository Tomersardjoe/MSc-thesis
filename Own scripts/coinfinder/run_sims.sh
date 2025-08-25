#!/usr/bin/env bash

for type in closed moderate open; do
  for i in {1..10}; do
    echo "Starting $type rep${i}..."
    outdir="${type}/rep${i}"

    if [ -d "$outdir" ] && [ "$(ls -A "$outdir")" ]; then
      echo "Error: $outdir already exists and is not empty. Stopping."
      exit 1
    fi

    mkdir -p "$outdir"
    (
      cd "$outdir" || exit
      coinfinder \
        -i "../../../simulation/output/${type}/grid_1/rep${i}/alpha_beta.tab" \
        -p "../../../simulation/output/${type}/grid_1/rep${i}/tree.nwk" \
        -a \
        -n \
        -E
    )
    echo "Finished $type rep${i}"
  done
done
