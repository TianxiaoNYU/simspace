#!/usr/bin/env python3
"""
Define one function per simulator. Each must:
- Take (shape, seed, mode, **extra)
- Run a single simulation
- Return a JSON with runtime_s and peak_ram_mb (MB)
"""
import argparse, json, time, sys, platform
from memory_profiler import memory_usage

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--method", required=True)
    p.add_argument("--shape", nargs=2, type=int, required=True)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--mode", type=str, default="ref-free")
    p.add_argument("--extra", type=str, default="{}")
    return p.parse_args()

def measure(func, *args, **kwargs):
    start = time.perf_counter()
    mem_trace, _ = memory_usage((func, args, kwargs),
                                interval=0.05, retval=True,
                                include_children=True, multiprocess=True)
    runtime = time.perf_counter() - start
    return runtime, max(mem_trace) if mem_trace else None

# --- Fill in your actual simulator calls below ---
def run_simspace(shape, seed, mode, **extra):
    import simspace as ss
    param = ss.util.generate_random_parameters(n_group=2, n_state=8, seed=seed)
    sim = ss.util.sim_from_params(
        param, shape=shape, num_iteration=4, n_iter=6,
        custom_neighbor=ss.spatial.generate_offsets(3, 'manhattan'),
        seed=seed
    )
    sim.create_omics(
        n_genes=220, 
        bg_ratio=0.6, 
        lr_ratio=0.2,
        bg_param = (1, 0.5), 
        marker_param = (5, 1.6), 
        spatial=False)

    return sim

def run_srt(shape, seed, mode, **extra):
    # TODO: add SRTsim call (Python or R wrapper)
    pass

def run_sccube(shape, seed, mode, **extra):
    import pandas as pd
    import scCube
    from scCube import scCube
    import pandas as pd
    import numpy as np
    import io, contextlib

    model = scCube()
    path = '/Users/zhaotianxiao/Library/CloudStorage/Dropbox/FenyoLab/Project/Spatialsim/SimSpace/examples/runtime_benchmark/data'
    sc_data = pd.read_csv(f'{path}/sccube_data_{shape[0]}.csv', index_col=0)
    sc_meta = pd.read_csv(f'{path}/sccube_meta_{shape[0]}.csv', index_col=False)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        _, generate_sc_meta_new = model.generate_pattern_random(
            generate_sc_data=sc_data,
            generate_sc_meta=sc_meta,
            set_seed=True,
            seed=seed,
            spatial_cell_type=None,
            spatial_dim=2,
            spatial_size=50,
            delta=10,
            lamda=0.75,)

    return generate_sc_meta_new


def run_scmultisim(shape, seed, mode, **extra):
    import subprocess
    nrow, ncol = shape
    cmd = ["Rscript", "run_scMultiSim.R", str(nrow), str(ncol), str(seed)]
    # Silence R’s stdout (it prints a lot)
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        print(proc.stderr, file=sys.stderr)
        raise RuntimeError(f"scMultiSim failed:\n{proc.stderr}")
    # Optionally print limited info
    return None

def run_scdesign3(shape, seed, mode, **extra):
    # TODO: add scDesign3 spatial simulation call
    pass

DISPATCH = {
    "simspace": run_simspace,
    "srt": run_srt,
    "sccube": run_sccube,
    "scmultisim": run_scmultisim,
    "scdesign3": run_scdesign3,
}

def main():
    args = parse_args()
    shape = (args.shape[0], args.shape[1])
    func = DISPATCH.get(args.method)
    if func is None:
        raise SystemExit(f"Unknown method: {args.method}")
    runtime, peak = measure(func, shape, args.seed, args.mode)
    print(json.dumps({
        "runtime_s": runtime,
        "peak_ram_mb": peak,
        "py_version": sys.version.split()[0],
        "platform": platform.platform()
    }))

if __name__ == "__main__":
    main()