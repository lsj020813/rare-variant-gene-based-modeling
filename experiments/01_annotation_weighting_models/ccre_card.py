
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import subprocess, json, collections, os
R=_config_path("${PROJECT_ROOT}/work/ref")
BT="bedtools"
OUT=f"{R}/l1"; os.makedirs(OUT, exist_ok=True)
tot=0; allc=collections.Counter()
for CH in [str(i) for i in range(1,23)]:
    if os.path.exists(f"{OUT}/chr{CH}.ccre_card.done"):
        continue
    out=subprocess.run(f"{BT} intersect -a {R}/b6_cards/chr{CH}.var.s.bed -b {R}/b6_cards/ccre.s.bed -wa -wb",
                       shell=True, capture_output=True, text=True)
    rows=collections.defaultdict(set)
    for line in out.stdout.splitlines():
        f=line.split("\t")
        rows[(f[0],f[1])].add(f[-1])
    with open(f"{OUT}/chr{CH}.ccre_card.tsv","w") as fh:
        fh.write("chrom\tpos0\tclasses\n")
        for (ch,p0),cs in sorted(rows.items(), key=lambda x:int(x[0][1])):
            fh.write(f"{ch}\t{p0}\t{','.join(sorted(cs))}\n")
            for cl in cs: allc[cl]+=1
    nvar=int(subprocess.run(f"wc -l < {R}/b6_cards/chr{CH}.var.s.bed", shell=True, capture_output=True, text=True).stdout.strip())
    lf=len(rows)/nvar
    ref=json.load(open(f"{R}/b6_cards/chr{CH}.b6.json"))
    assert abs(lf-ref["labeled_frac"])<0.005, f"GATE FAIL chr{CH}: {lf} vs {ref['labeled_frac']}"
    open(f"{OUT}/chr{CH}.ccre_card.done","w").write("ok\n")
    tot+=len(rows)
    print(f"chr{CH}: labeled {len(rows):,}/{nvar:,} ({lf:.4f}) GATE OK", flush=True)
print(f"TOTAL labeled variants: {tot:,}")
print("classes:", dict(allc))
print("L1_CCRE_CARD_COMPLETE")
