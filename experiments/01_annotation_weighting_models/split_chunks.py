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


import argparse, os, sys
from pathlib import Path
sys.path.insert(0, _config_path('${PROJECT_ROOT}/work/run_l3'))
import l3_common as C

ap = argparse.ArgumentParser()
ap.add_argument('--group-files', required=True, help='comma list of arm group files')
ap.add_argument('--genes-per-chunk', required=True, type=int)
ap.add_argument('--out', required=True)
a = ap.parse_args()
out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
report = {}
for gf in a.group_files.split(','):
    tag = Path(gf).stem
    blocks, cur, slots = [], [], 0
    for line in open(gf):
        f = line.split()
        if len(f) < 3:
            continue
        if f[1] == 'var':
            if cur:
                blocks.append(cur)
            cur = [line]; slots += len(f) - 2
        else:
            cur.append(line)
            if f[1] == 'anno':
                assert set(f[2:]) == {'all'}, f'label not all: {f[:4]}'
    if cur:
        blocks.append(cur)
    for old in out.glob(f'{tag}.part*.txt'):
        old.unlink()
    per = a.genes_per_chunk
    n = (len(blocks) + per - 1) // per
    wg = ws = 0
    for i in range(n):
        chunk = blocks[i * per:(i + 1) * per]
        with open(out / f'{tag}.part{i+1:03d}.txt', 'w') as fh:
            for b in chunk:
                fh.write(''.join(b)); wg += 1; ws += len(b[0].split()) - 2
    assert wg == len(blocks) and ws == slots, f'{tag}: union {wg}/{ws} vs {len(blocks)}/{slots}'
    report[tag] = dict(genes=len(blocks), slots=slots, chunks=n)
    print(f'{tag}: chunks {n} genes {wg} slots {ws}')
C.load_libraries(1)
C.atomic_json(out / 'split_manifest.json', dict(genes_per_chunk=a.genes_per_chunk, arms=report))
print('SPLIT_OK')
