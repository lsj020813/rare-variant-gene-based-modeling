#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
W=${PROJECT_ROOT}/work/run_trackA
python3 - <<'EOF'
import os as _env_os, re as _env_re
def _env_path(value):
    return _env_re.sub(r"\$\{([A-Z][A-Z0-9_]*)\}", lambda match: _env_os.environ[match.group(1)], value)
import glob, csv, os
W=_env_path("${PROJECT_ROOT}/work/run_trackA"); n=0; eff=0; cs_eff=0; mp_eff=0; gate_ok=0
for f in sorted(glob.glob(f'{W}/fm/*.arms.tsv')):
    rows=list(csv.DictReader(open(f),delimiter='\t')); d={r['arm']:r for r in rows}
    if 'phi' not in d: continue
    perms=[r for r in rows if r['arm'].startswith('perm')]
    if not perms: continue
    n+=1
    g=d.get('gate_uniform')
    if g is None:
        gf=f.replace('.arms.tsv','.gate.arms.tsv')
        if os.path.exists(gf): g={r['arm']:r for r in csv.DictReader(open(gf),delimiter='\t')}.get('gate_uniform')
    if g and g['gate_max_abs_dpip'] not in ('NA','') and float(g['gate_max_abs_dpip'])<1e-6: gate_ok+=1
    def f_(x): return float('nan') if x in ('NA','') else float(x)
    phi_cs=f_(d['phi']['mean_cs_size']); phi_mp=f_(d['phi']['max_pip'])
    pc=[f_(r['mean_cs_size']) for r in perms]; pm=[f_(r['max_pip']) for r in perms]
    # 백분위: φ 가 perm 분포에서 차지하는 위치 (CS 합 크기 작을수록·max PIP 클수록 우수)
    import math
    cs_pct=0.0 if math.isnan(phi_cs) else sum(1 for x in pc if (not math.isnan(x)) and x>phi_cs)/len(pc); mp_pct=sum(1 for x in pm if x<phi_mp)/len(pm)
    c=cs_pct>=0.95; m=mp_pct>=0.95; cs_eff+=c; mp_eff+=m; eff+= (c or m)
print(f'n={n} gate_ok={gate_ok} eff(cs|maxpip)={eff} ({eff/n:.2f}) cs_only={cs_eff} maxpip_only={mp_eff}' if n else 'n=0')
EOF
