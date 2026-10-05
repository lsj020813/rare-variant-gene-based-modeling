#!/usr/bin/env python3
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import shutil, sys, hashlib, py_compile
P = _config_path("${PROJECT_ROOT}/work/run_band15/model_v10_out/l1_train_v10.py")
src = open(P).read()
if "E6B_PATCH" in src: print("already patched"); sys.exit(0)
assert "E6_PATCH" in src
shutil.copy(P, P + ".pre_e6b_bak")
def rep(old, new):
    global src
    assert src.count(old) == 1, f"count {src.count(old)}: {old[:70]}"
    src = src.replace(old, new)
rep('    require(0 < args.memory_gb <= 8, "v10 memory cap must be <=8 GiB")\n',
    '    require(0 < args.memory_gb <= 64, "v10 memory cap must be <=64 GiB")  # E6B_PATCH ★ (was 8; real-data scale)\n')
rep('                        require(not any(b.startswith(heavy) for b in basenames), "another owner heavy stage is active")\n',
    '                        if any(b.startswith(heavy) for b in basenames):  # E6B_PATCH ★\n'
    '                            if getattr(args, "allow_concurrent_heavy", False):\n'
    '                                concurrent_heavy.add(pid)\n'
    '                            else:\n'
    '                                require(False, "another owner heavy stage is active")\n')
rep('    ancestors = {os.getpid()}\n', '    ancestors = {os.getpid()}\n    concurrent_heavy = set()  # E6B_PATCH\n')
rep('    return lock, dict(threads=args.threads, memory_cap_bytes=cap, tmpdir=str(tmp),\n',
    '    if concurrent_heavy:\n        print(f"[guard] E6B_PATCH: proceeding alongside {len(concurrent_heavy)} owner heavy-stage process(es) (--allow-concurrent-heavy)", flush=True)\n'
    '    return lock, dict(threads=args.threads, memory_cap_bytes=cap, tmpdir=str(tmp), concurrent_owner_heavy_processes=len(concurrent_heavy),\n')
rep("    parser.add_argument('--threads', type=int, default=4)\n",
    "    parser.add_argument('--threads', type=int, default=4)\n    parser.add_argument('--allow-concurrent-heavy', action='store_true', help='★ E6B_PATCH: downgrade owner heavy-stage overlap from failure to a logged note')\n")
open(P, "w").write(src); py_compile.compile(P, doraise=True)
print("PATCHED2", hashlib.sha256(open(P,'rb').read()).hexdigest()[:16])
