
import argparse
import glob
import hashlib
import json
import os
import subprocess
import sys
import time

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
CODE_FILES = ["__init__.py", "state.py", "phase.py", "align.py", "workflow.py",
              "annot.py", "tiles.py", "vcfio.py", "fastpath.py", "build_states.py",
              "sample_tiles.py", "split.py", "bench_io.py", "workflow_ext.py"]

OUTPUTS = {
    "state_inventory": "haplotype_state_inventory.csv",
    "diplotype_contrast": "diplotype_contrast_support.csv",
    "phase_pairs": "phase_pair_support.csv",
    "repr_info": "representation_information.csv",
    "tile_status": "tile_status.csv",
}

def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def code_hashes():
    out = {}
    for f in CODE_FILES:
        p = os.path.join(HERE, f)
        if os.path.exists(p):
            out[f] = sha256_file(p)
    return out

def build_states(args):
    t0 = time.time()
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.shard_dir, exist_ok=True)
    os.makedirs(args.tmpdir, exist_ok=True)
    procs = []
    logs = []
    for s in range(args.workers):
        log = os.path.join(args.outdir, "worker%02d.log" % s)
        logs.append(log)
        cmd = ["nice", "-n", "19", "ionice", "-c3", sys.executable, "-m", "hcf.build_states",
               "--tilelist", args.tilelist, "--vardir", args.vardir, "--split", args.split,
               "--outdir", args.shard_dir, "--tmpdir", args.tmpdir,
               "--shard", str(s), "--nshards", str(args.workers)]
        if args.limit:
            cmd += ["--limit", str(args.limit)]
        fh = open(log, "w")
        env = dict(os.environ)
        for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                  "NUMEXPR_NUM_THREADS"):
            env[v] = "1"
        procs.append((subprocess.Popen(cmd, stdout=fh, stderr=subprocess.STDOUT, env=env), fh))
    rcs = []
    for p, fh in procs:
        rcs.append(p.wait())
        fh.close()

    merged = {}
    for key, fname in OUTPUTS.items():
        parts = []
        for p in sorted(glob.glob(os.path.join(args.shard_dir, "%s.shard*.csv" % key))):
            if os.path.getsize(p) > 0:
                parts.append(pd.read_csv(p))
        if parts:
            df = pd.concat(parts, ignore_index=True)
            sort_cols = [c for c in ("chrom", "tile_start", "subwindow_id") if c in df.columns]
            if sort_cols:
                df = df.sort_values(sort_cols).reset_index(drop=True)
            df.to_csv(os.path.join(args.outdir, fname), index=False)
            merged[fname] = int(len(df))
        else:
            merged[fname] = 0

    done = {
        "stage": "build-states", "run_id": "HCF-20260922-v1/run_20260922T052403Z",
        "status": "OK" if all(r == 0 for r in rcs) else "PARTIAL",
        "worker_return_codes": rcs, "workers": args.workers,
        "config": args.config, "config_sha256": sha256_file(args.config) if args.config else None,
        "tilelist": args.tilelist, "tilelist_sha256": sha256_file(args.tilelist),
        "split_manifest": args.split, "split_sha256": sha256_file(args.split),
        "code_sha256": code_hashes(),
        "row_counts": merged,
        "wall_s": round(time.time() - t0, 1),
        "recomputed": True,
        "note": ("hcf/workflow.py was NOT modified; its build-states still returns "
                 "NOT_IMPLEMENTED. This extension entry point implements the stage."),
    }
    with open(os.path.join(args.outdir, "done.build_states.json"), "w") as fh:
        json.dump(done, fh, indent=2, ensure_ascii=False)
    print(json.dumps(done, indent=2, ensure_ascii=False))
    return 0 if done["status"] == "OK" else 3

def main(argv=None):
    ap = argparse.ArgumentParser("hcf.workflow_ext")
    sub = ap.add_subparsers(dest="stage")
    b = sub.add_parser("build-states")
    b.add_argument("--config", default="")
    b.add_argument("--tilelist", required=True)
    b.add_argument("--vardir", required=True)
    b.add_argument("--split", required=True)
    b.add_argument("--outdir", required=True)
    b.add_argument("--shard-dir", required=True)
    b.add_argument("--tmpdir", required=True)
    b.add_argument("--workers", type=int, default=4)
    b.add_argument("--limit", type=int, default=0)
    args = ap.parse_args(argv)
    if args.stage == "build-states":
        return build_states(args)
    ap.print_help()
    return 2

if __name__ == "__main__":
    sys.exit(main())
