#!/usr/bin/env python3
# Compile: not needed (Python script)
# Usage: python3 autotune.py

import subprocess, time, csv, itertools

# Benchmark parameters
nrepeats = [10, 20, 30, 40, 50]
stddevs = [0, 0.01, 0.1]

# Computation parameters
threads = [1, 2, 4, 8, 16]
cells = [10, 15, 25]
gpu_flags = [
    {"useGPU": True, "fillGPU": True, "pseudoCoal": True},
    {"useGPU": True, "fillGPU": True, "pseudoCoal": False},
    {"useGPU": True, "fillGPU": False, "pseudoCoal": True},
    {"useGPU": True, "fillGPU": False, "pseudoCoal": False},
    {"useGPU": False},  # CPU mode
]
algo_flags = [
    {"useLocalHist": True, "pseudoCoal": False},
    {"useLocalHist": False, "pseudoCoal": True},
    {"useLocalHist": True, "pseudoCoal": True},
    {"useLocalHist": False, "pseudoCoal": False},
]

def run_benchmark(params):
    cmd = ["aesdebye"]  # replace with your executable
    for k, v in params.items():
        if isinstance(v, bool):
            if v: cmd.append(f"--{k}")
        else:
            cmd += [f"--{k}", str(v)]
    start = time.time()
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return time.time() - start

results = []

for nr, std, nt, gpu, algo in itertools.product(nrepeats, stddevs, threads, gpu_flags, algo_flags):
    params = {"nRepeats": nr, "stdDev": std, "nThreads": nt, **gpu, **algo}

    # Skip irrelevant params
    if params.get("useGPU", False):
        params.pop("nCells", None)  # GPU ignores nCells
    else:
        for nc in cells:
            params["nCells"] = nc
            runtime = run_benchmark(params)
            if runtime > 0.2:  # ignore very short runs
                results.append({**params, "runtime": runtime})
        continue

    # GPU case
    runtime = run_benchmark(params)
    if runtime > 0.2:
        results.append({**params, "runtime": runtime})

# Save results
with open("benchmarks.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=results[0].keys())
    writer.writeheader()
    writer.writerows(results)

print(f"✅ Completed {len(results)} runs, saved to benchmarks.csv")
