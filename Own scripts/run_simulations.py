import os
import subprocess
import shutil
from simulation import simulate

# Number of replicates to run for each scenario
num_reps = 3

# Base directory where all scenario outputs will be stored
base_dir = "output"

scenarios = [
    {
        "name": "closed",
        "num_taxa": 10,
        "psize": [320, 75],
        "gsize": [300, 30],
        "lamb": 1,
        "mu": 0,
        "scale_events": [0.05, 0.15],
        "num_pairs": [0, 12],
        "selection": 1,
        "fluidity_range": (0.03, 0.06),
        "alpha_range": (1.02, float("inf"))
    },
    {
        "name": "moderate",
        "num_taxa": 250,
        "psize": [720, 160],
        "gsize": [340, 50],
        "lamb": 1,
        "mu": 0,
        "scale_events": [0.12, 0.40],
        "num_pairs": [100, 8],
        "selection": 5.8,
        "fluidity_range": (0.10, 0.18),
        "alpha_range": (0.95, 1.02)
    },
    {
        "name": "open",
        "num_taxa": 250,
        "psize": [820, 320],
        "gsize": [300, 65],
        "lamb": 1,
        "mu": 0,
        "scale_events": [0.65, 3.5],
        "num_pairs": [100, 14],
        "selection": 2.2,
        "fluidity_range": (0.22, 0.32),
        "alpha_range": (0, 0.90)
    }
]

def calc_fluidity(csv_path):
    """
    Run the R script and return genomic fluidity and alpha values.

    Returns:
        (fluidity, alpha_all, alpha_slow, alpha_fast)
    """
    result = subprocess.run(
        ["Rscript", "--vanilla", "fluidity_calc.R", csv_path],
        capture_output=True, text=True, check=True
    )

    fluidity = None
    alpha_all = None
    alpha_slow = None
    alpha_fast = None

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


for s in scenarios:
    rep = 1
    while rep <= num_reps:
        outdir = os.path.join(base_dir, s['name'], f"rep{rep}")
        os.makedirs(outdir, exist_ok=True)
        print(f"Running simulation: {s['name']} (rep {rep})")

        simulate(
            outdir,
            num_taxa=s["num_taxa"],
            psize=s["psize"],
            gsize=s["gsize"],
            lamb=s["lamb"],
            scale_events=s["scale_events"],
            num_pairs=s["num_pairs"],
            selection=s["selection"],
            mu=s["mu"]
        )

        csv_path = os.path.join(outdir, "gene_presence_absence_all.csv")
        fluidity, alpha_all, alpha_slow, alpha_fast = calc_fluidity(csv_path)

        if fluidity is None or alpha_all is None:
            print("    Could not parse metrics; discarding and retrying...")
            shutil.rmtree(outdir)
            continue

        min_f, max_f = s["fluidity_range"]
        min_a, max_a = s["alpha_range"]

        print(f"    Fluidity = {fluidity:.6f}, "
              f"Alpha(all) = {alpha_all:.6f}, "
              f"Alpha(slow) = {alpha_slow:.6f}, "
              f"Alpha(fast) = {alpha_fast:.6f}")

        if (min_f <= fluidity <= max_f) and (min_a <= alpha_all <= max_a):
            rep += 1
        else:
            print(f"    Outside target ranges {s['fluidity_range']} & {s['alpha_range']}, discarding...")
            shutil.rmtree(outdir)

print("All simulations complete.")
