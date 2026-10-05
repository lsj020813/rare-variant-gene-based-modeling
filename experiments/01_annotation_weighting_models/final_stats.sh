#!/bin/bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -o pipefail
cd ${PROJECT_ROOT}/work/run_trackB
nice -n 10 ./venv/bin/python diag_trackB_v2.py --set ctrl --B 2000 --out diag_ctrl.json > diag_ctrl.log 2>&1
nice -n 10 ./venv/bin/python diag_trackB_v2.py --set all --B 2000 --out diag_all.json > diag_all.log 2>&1
./venv/bin/python - <<'EOF'
import os as _env_os, re as _env_re
def _env_path(value):
    return _env_re.sub(r"\$\{([A-Z][A-Z0-9_]*)\}", lambda match: _env_os.environ[match.group(1)], value)
import json, re, hashlib
D=_env_path("${PROJECT_ROOT}/work/run_trackB")
g=[l for l in open(f"{D}/gpu_mon.log") if re.match(r"\d\d:\d\d:\d\d \d+ MiB",l)]
mem=[int(l.split()[1]) for l in g]; util=[int(l.split(",")[1].strip().split()[0]) for l in g]; temp=[int(l.split(",")[2]) for l in g]
w=[l for l in open(f"{D}/watchdog.log") if " n=" in l]
def fld(l,k):
    m=re.search(k+r"=([0-9.]+)",l); return float(m.group(1)) if m else None
loads=[fld(l,"load") for l in w if fld(l,"load") is not None]; rss=[fld(l,"rss_max") for l in w if fld(l,"rss_max") is not None]; cores=[fld(l,"our_cores") for l in w if fld(l,"our_cores") is not None]
t=sorted(float(l.split("\t")[3]) for l in open(f"{D}/scores_full.tsv").readlines()[1:] if l.split("\t")[1]=="OK")
st={"gpu_samples":len(g),"gpu_mem_MiB_max":max(mem),"gpu_util_mean":sum(util)/len(util),"gpu_util_p10":sorted(util)[len(util)//10],"gpu_temp_max":max(temp),
 "wd_samples":len(w),"load1_max":max(loads),"load1_p95":sorted(loads)[int(len(loads)*0.95)],"rss_max_G":max(rss),"our_cores_max":max(cores),"breach":sum(1 for l in w if "breach=1" in l),
 "n_ok":len(t),"t_var_mean":sum(t)/len(t),"t_var_median":t[len(t)//2],"t_var_p95":t[int(len(t)*0.95)],"t_var_max":t[-1],"wall_start":"15:02:49","wall_end":"17:49:19",
 "scores_full_sha256":hashlib.sha256(open(f"{D}/scores_full.tsv","rb").read()).hexdigest(),"variants_ordered_sha256":hashlib.sha256(open(f"{D}/variants_ordered.tsv","rb").read()).hexdigest()}
json.dump(st,open(f"{D}/run_stats.json","w"),indent=1); print(json.dumps(st))
EOF
echo EXIT $?; touch final_stats.done
