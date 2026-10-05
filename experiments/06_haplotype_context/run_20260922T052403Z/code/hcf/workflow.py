import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

IMPLEMENTED = {"inventory", "test", "phase-audit", "io-bench"}
ALL_COMMANDS = [
    "inventory", "test", "phase-audit", "build-states", "fit-development",
    "calibrate", "freeze-evaluation", "evaluate", "controls", "report", "run",
    "io-bench",
]

def _sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def cmd_inventory(args):
    print(json.dumps({
        "status": "SEE_STAGE1",
        "note": "22/22-chromosome header/index/GT-sample inventory was produced "
                "in stage-1 registration (hcf_inventory.py / input_manifest.json). "
                "This stage's `inventory` subcommand exists for CLI-contract "
                "completeness and re-validates the frozen config hash only.",
        "config": args.config,
        "config_sha256": _sha256_file(args.config) if Path(args.config).exists() else None,
    }))
    return 0

def cmd_test(args):
    from tests.test_core import ALL_TESTS
    results = []
    n_pass = 0
    t0 = time.time()
    for fn in ALL_TESTS:
        name = fn.__name__
        try:
            fn()
            results.append({"test": name, "status": "PASS"})
            n_pass += 1
        except Exception as e:
            results.append({"test": name, "status": "FAIL", "error": f"{type(e).__name__}: {e}"})
    elapsed = time.time() - t0
    summary = {
        "status": "ALL_PASS" if n_pass == len(ALL_TESTS) else "SOME_FAILED",
        "n_total": len(ALL_TESTS), "n_pass": n_pass, "n_fail": len(ALL_TESTS) - n_pass,
        "elapsed_s": round(elapsed, 4), "results": results,
    }
    print(json.dumps(summary, indent=1))
    return 0 if n_pass == len(ALL_TESTS) else 1

def cmd_phase_audit(args):
    from hcf.phase_audit import run_phase_audit
    return run_phase_audit(args.config, args.out_dir, args.bcftools, args.chroms,
                            n_records=args.n_records, n_samples=args.n_samples, seed=args.seed)

def cmd_io_bench(args):
    from hcf.io_bench import run_io_bench
    return run_io_bench(args.config, args.out_dir, args.bcftools, args.n_tiles,
                         n_people=args.n_people, seed=args.seed)

def _not_implemented(name):
    def _cmd(args):
        print(json.dumps({
            "status": "NOT_IMPLEMENTED",
            "command": name,
            "reason": "Stage-2 scope is HC-D0 phase audit + 12 unit tests + IO "
                      "benchmark only (brief SS6.3-6.4, SS7, SS9.2, SS15.1, SS18, "
                      "SS20). HC-D1 state construction, phenotype reading, and "
                      "model fitting are explicitly out of scope this stage.",
        }))
        return 2
    return _cmd

def main(argv=None):
    p = argparse.ArgumentParser(prog="python -m hcf.workflow")
    sub = p.add_subparsers(dest="command", required=True)

    p_inv = sub.add_parser("inventory"); p_inv.add_argument("--config", required=True)
    p_inv.set_defaults(func=cmd_inventory)

    p_test = sub.add_parser("test"); p_test.add_argument("--config", required=True)
    p_test.set_defaults(func=cmd_test)

    p_pa = sub.add_parser("phase-audit")
    p_pa.add_argument("--config", required=True)
    p_pa.add_argument("--out-dir", required=True)
    p_pa.add_argument("--bcftools", default="bcftools")
    p_pa.add_argument("--chroms", default="1-22")
    p_pa.add_argument("--n-records", type=int, default=2000)
    p_pa.add_argument("--n-samples", type=int, default=1000)
    p_pa.add_argument("--seed", type=int, default=20260922)
    p_pa.set_defaults(func=cmd_phase_audit)

    p_io = sub.add_parser("io-bench")
    p_io.add_argument("--config", required=True)
    p_io.add_argument("--out-dir", required=True)
    p_io.add_argument("--bcftools", default="bcftools")
    p_io.add_argument("--n-tiles", type=int, default=12)
    p_io.add_argument("--n-people", type=int, default=2048)
    p_io.add_argument("--seed", type=int, default=20260922)
    p_io.set_defaults(func=cmd_io_bench)

    for name in ["build-states", "fit-development", "calibrate", "freeze-evaluation",
                 "evaluate", "controls", "report", "run"]:
        pn = sub.add_parser(name)
        pn.add_argument("--config", required=False)
        pn.add_argument("--resume", action="store_true")
        pn.set_defaults(func=_not_implemented(name))

    args = p.parse_args(argv)
    return args.func(args)

if __name__ == "__main__":
    sys.exit(main())
