import os
import math
import shutil
import subprocess
import csv
import itertools
from simulation import simulate

num_reps = 10
base_dir = "output"
summary_file = os.path.join(base_dir, "grid_summary.csv")

scenarios = [
      {
          "name": "closed",
          "num_taxa": 150,
          "psize": [320, 75],
          "gsize_fraction_range": [0.9],
          "lamb": 1,
          "mu": 0,
          "scale_events_lower_range": [0.5],
          "scale_events_upper_range": [1.2],
          "num_pairs_slow_range": [0],
          "num_pairs_fast_range": [24],
          "selection_range": [2.5],
          "fluidity_range": (0.03, 0.06),
          "alpha_range": (1.2, 1.80)
      },
      {
          "name": "moderate",
          "num_taxa": 150,
          "psize": [720, 220],
          "gsize_fraction_range": [0.55],
          "lamb": 1,
          "mu": 0,
          "scale_events_lower_range": [1.0],
          "scale_events_upper_range": [3.5],
          "num_pairs_slow_range": [8],
          "num_pairs_fast_range": [100],
          "selection_range": [1.5],
          "fluidity_range": (0.10, 0.18),
          "alpha_range": (0.90, 0.95)
      },
      {
        "name": "open",
        "num_taxa": 150,
        "psize": [400, 120],
        "gsize_fraction_range": [0.28],
        "lamb": 1,
        "mu": 0,
        "scale_events_lower_range": [1.0],
        "scale_events_upper_range": [5.0],
        "num_pairs_slow_range": [20],
        "num_pairs_fast_range": [60],
        "selection_range": [1.7],
        "fluidity_range": [0.20, 0.30],
        "alpha_range": [0.05, 0.85]
      }
]

def calc_fluidity(csv_path):
    result = subprocess.run(
        ["Rscript", "--vanilla", "fluidity_calc.R", csv_path],
        capture_output=True, text=True, check=True
    )
    fluidity = None
    alpha_all = None
    alpha_slow = math.nan
    alpha_fast = math.nan

    for line in result.stdout.splitlines():
        if line.startswith("Genomic fluidity:"):
            fluidity = float(line.split(":")[1].strip())
        elif line.startswith("Pangenome openness - all genes (alpha):"):
            alpha_all = float(line.split(":")[1].strip())
        elif line.startswith("Pangenome openness - slow partition (alpha):"):
            alpha_slow = float(line.split(":")[1].strip())
        elif line.startswith("Pangenome openness - fast partition (alpha):"):
            alpha_fast = float(line.split(":")[1].strip())

    return fluidity, alpha_all, alpha_slow, alpha_fast

os.makedirs(base_dir, exist_ok=True)
with open(summary_file, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow([
        "scenario", "grid_id", "rep",
        "gsize_frac", "scale_events_lower", "scale_events_upper",
        "num_pairs_slow", "num_pairs_fast", "selection",
        "fluidity", "alpha_all", "alpha_slow", "alpha_fast"
    ])

for s in scenarios:
    grid = list(itertools.product(
        s["gsize_fraction_range"],
        s["scale_events_lower_range"],
        s["scale_events_upper_range"],
        s["num_pairs_slow_range"],
        s["num_pairs_fast_range"],
        s["selection_range"]
    ))

    num_combos = len(grid)
    total_runs = num_combos * num_reps

    print(f"Scenario '{s['name']}': {num_combos} parameter combinations x {num_reps} replicates = {total_runs} total simulation runs")

    for grid_idx, (gfrac, se_low, se_up, np_slow, np_fast, sel) in enumerate(grid, start=1):
        gsize = [int(ps * gfrac) for ps in s["psize"]]
    
        if se_up <= se_low:
            continue
        if any(gm >= ps for gm, ps in zip(gsize, s["psize"])):
            continue
    
        rep = 1
        while rep <= num_reps:
            outdir = os.path.join(base_dir, s['name'], f"grid_{grid_idx}", f"rep{rep}")
            os.makedirs(outdir, exist_ok=True)
            try:
                simulate(
                    outdir,
                    num_taxa=s["num_taxa"],
                    psize=s["psize"],
                    gsize=gsize,
                    lamb=s["lamb"],
                    scale_events=[se_low, se_up],
                    num_pairs=[np_slow, np_fast],
                    selection=sel,
                    mu=s["mu"]
                )

                csv_path = os.path.join(outdir, "gene_presence_absence_all.csv")
                if not os.path.exists(csv_path):
                    print(f"    Missing output file for {s['name']} grid {grid_idx} rep {rep}")
                    continue

                fluidity, alpha_all, alpha_slow, alpha_fast = calc_fluidity(csv_path)
                
                if fluidity is None or alpha_all is None:
                    print("    Could not parse metrics; discarding and retrying...")
                    shutil.rmtree(outdir)
                    continue
                    
                min_f, max_f = s["fluidity_range"]
                min_a, max_a = s["alpha_range"]
        
                print(f"    Fluidity = {fluidity}, Alpha(all) = {alpha_all}")
                
                if (min_f <= fluidity <= max_f) and (min_a <= alpha_all <= max_a):
                    with open(summary_file, "a", newline="") as f:
                        writer = csv.writer(f)
                        writer.writerow([
                            s["name"], grid_idx, rep,
                            gfrac, se_low, se_up, np_slow, np_fast, sel,
                            fluidity, alpha_all, alpha_slow, alpha_fast
                        ])
                    rep += 1
                else:
                    print(f"    Outside target ranges {s['fluidity_range']} & {s['alpha_range']}, discarding...")
                    shutil.rmtree(outdir)
    
            except ValueError as e:
                if "Total of weights" in str(e):
                    print(f"    Skipping: all rates zero for {s['name']} grid {grid_idx} rep {rep}")
                    shutil.rmtree(outdir)
                    rep += 1
                    continue
                else:
                    raise

print(f"All grid search simulations complete. Summary saved to {summary_file}")
