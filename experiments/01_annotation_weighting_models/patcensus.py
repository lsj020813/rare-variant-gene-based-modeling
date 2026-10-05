
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import collections, glob, os
R=_config_path("${PROJECT_ROOT}/work/ref")
combo=collections.Counter()
for ch in range(1,23):
    fp=f"{R}/l1/chr{ch}.ccre_card.tsv"
    with open(fp) as fh:
        fh.readline()
        for line in fh:
            combo[line.rstrip("\n").split("\t")[2]] += 1
print(f"distinct class combos (labeled variants): {len(combo)}")
print("top 12:", combo.most_common(12))
tot=sum(combo.values())
print(f"labeled variants total: {tot:,}")
for CH in ("21","22"):
    card={}
    with open(f"{R}/l1/chr{CH}.ccre_card.tsv") as fh:
        fh.readline()
        for line in fh:
            f=line.rstrip("\n").split("\t")
            card[f[1]]=f[2]
    occ=set(); genes=set(); nv=0
    import subprocess
    BCF="bcftools"
    out=subprocess.run(f"{BCF} query -f '%POS\t%ID\n' {R}/lift38_keyed/chr{CH}.keyed38.vcf.gz",
                       shell=True, capture_output=True, text=True)
    k2p={}
    for l in out.stdout.splitlines():
        p,kid=l.split("\t")
        k2p[kid[3:] if kid.startswith("chr") else kid]=str(int(p)-1)
    for line in open(f"{R}/groupfiles_bwg/chr{CH}.B_3kb_re2g.txt"):
        f=line.split()
        if f[1]!="var": continue
        g=f[0]; genes.add(g)
        for k in f[2:]:
            k=k[3:] if k.startswith("chr") else k
            nv+=1
            p0=k2p.get(k)
            cls=card.get(p0,"none") if p0 else "none"
            occ.add((g,cls))
    print(f"chr{CH}: genes {len(genes):,}  var-slots {nv:,}  occupied (gene,combo) {len(occ):,}  -> per gene {len(occ)/len(genes):.2f}")
print("PATTERN_CENSUS_DONE")
