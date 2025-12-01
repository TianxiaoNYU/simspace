#!/usr/bin/env python3
"""
Benchmark harness for spatial simulators: runtime + peak RAM.

Usage example:
  python bench.py --methods simspace srt sccube scmultisim scdesign3 \
                  --grid-sizes 30 50 75 100 150 200 \
                  --repeats 10 --outfile results.csv
"""
import argparse, time, csv, os, sys, subprocess, json
from pathlib import Path

def parse_args():
    p = argparse.ArgumentParser(description="Benchmark spatial simulators (runtime + peak RAM)")
    p.add_argument("--methods", nargs="+", required=True)
    p.add_argument("--grid-sizes", nargs="+", type=int, required=True)
    p.add_argument("--repeats", type=int, default=5)
    p.add_argument("--outfile", type=str, default="results.csv")
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--mode", type=str, default="ref-free", choices=["ref-free", "ref-based"])
    p.add_argument("--extra", type=str, default="{}", help="JSON string of extra kwargs to pass to runners")
    return p.parse_args()

def env_with_threads(n):
    env = os.environ.copy()
    for k in ["OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","VECLIB_MAXIMUM_THREADS","NUMEXPR_NUM_THREADS"]:
        env[k] = str(n)
    return env

def run_once(method, n, seed, mode, extra_kws, pybin=sys.executable, threads=1):
    """Runs a single simulation in a fresh subprocess and captures runtime + peak RAM."""
    cmd = [
        pybin, "runners.py",
        "--method", method,
        "--shape", str(n), str(n),
        "--seed", str(seed),
        "--mode", mode,
        "--extra", json.dumps(extra_kws)
    ]
    env = env_with_threads(threads)
    start = time.perf_counter()
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env)
    end = time.perf_counter()

    if proc.returncode != 0:
        raise RuntimeError(f"{method} failed: {proc.stderr}")
    try:
        payload = json.loads(proc.stdout.strip())
    except Exception:
        raise RuntimeError(f"Runner did not return JSON:\n{proc.stdout}\n{proc.stderr}")
    runtime = payload.get("runtime_s", end - start)
    ram = payload.get("peak_ram_mb", None)
    return runtime, ram

def main():
    args = parse_args()
    extra = json.loads(args.extra)
    outpath = Path(args.outfile)
    outpath.parent.mkdir(parents=True, exist_ok=True)

    fields = ["method","mode","grid_n","repeat","runtime_s","peak_ram_mb"]
    with open(outpath, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()

        for method in args.methods:
            for grid_n in args.grid_sizes:
                for r in range(args.repeats):
                    t, m = run_once(method, grid_n, r, args.mode, extra, threads=args.threads)
                    writer.writerow(dict(method=method, mode=args.mode,
                                         grid_n=grid_n, repeat=r,
                                         runtime_s=t, peak_ram_mb=m))
                    print(f"[{method}] {grid_n}x{grid_n}, run {r}: {t:.2f}s | {m:.1f} MB")

if __name__ == "__main__":
    main()