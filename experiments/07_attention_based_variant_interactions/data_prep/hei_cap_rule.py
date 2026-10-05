import os as _os
import math as _number_math
def _required_number(name, cast, positive=False):
    raw = _os.environ.get(name, "")
    if not raw.strip():
        raise ValueError(name + " must be set and nonblank")
    try:
        value = cast(raw)
    except (ValueError, OverflowError):
        raise ValueError(name + " has an invalid numeric value") from None
    if isinstance(value, float) and not _number_math.isfinite(value):
        raise ValueError(name + " must be finite")
    if positive and value <= 0:
        raise ValueError(name + " must be positive")
    return value
if not _os.environ.get("PROJECT_ROOT", "").strip():
    raise ValueError("PROJECT_ROOT must be set and nonblank")
N_TOTAL = _required_number("N_TOTAL", int, True)
N_TRAIN = _required_number("N_TRAIN", int, True)

import csv, collections, json
rows=list(csv.DictReader(open((_os.environ["PROJECT_ROOT"] + '/work/prs/out/hei/sel/gene_tokens.tsv')),delimiter="\t"))
TCHL=42*1024**2; out={}
for cap in (64,128,256):
    sel=[r for r in rows if int(r["rank"])<cap]; n=collections.Counter(r["gene"] for r in sel)
    uk=len({r["key19"] for r in sel}); cost=sum(v*v for v in n.values())
    out[cap]=dict(tokens_dup=len(sel),unique=uk,ram_gb=round(uk*N_TRAIN*2/1e9*N_TOTAL/N_TRAIN,1),attn_rel=round(cost/TCHL,2),ok=(uk*N_TOTAL*2/1e9<=42 and cost/TCHL<=1.5))
chosen=max([c for c,v in out.items() if v["ok"]] or [64]); out["chosen_cap"]=chosen
json.dump(out,open((_os.environ["PROJECT_ROOT"] + '/work/prs/out/hei/sel/cap_rule.json'),"w"),indent=1); print(json.dumps(out))
