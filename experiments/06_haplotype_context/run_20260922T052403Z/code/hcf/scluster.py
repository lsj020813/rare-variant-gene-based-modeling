
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import argparse
import json
import os
import subprocess
import sys
import time

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

BT = "bcftools"
ORIG = _config_path("${PROJECT_ROOT}/work/ref/orig_index/chr%s.vcf.gz")
HAPLA = None

def sh(cmd):
    p = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if p.returncode != 0:
        raise RuntimeError("cmd failed rc=%d: %s\n%s" % (p.returncode, cmd, p.stdout.decode()[-2000:]))
    return p.stdout.decode()

def read_bca(path, n_hap):
    raw = np.fromfile(path, dtype=np.uint8)
    hdr = len(raw) % n_hap
    if hdr not in (0, 1, 2, 3, 4, 8):
        raise RuntimeError("unexpected .bca size %d for %d haplotypes" % (len(raw), n_hap))
    arr = raw[hdr:].reshape(-1, n_hap)
    return arr

def read_ids(path):
    with open(path) as fh:
        return [l.strip() for l in fh if l.strip()]

def make_bcf(chrom, start, posfile, idsfile, out):
    sh("%s view -r %s:%d-%d -T %s -S %s -Ou %s | %s annotate -x FORMAT/GP,FORMAT/DS,INFO -Ob -o %s"
       % (BT, chrom, start + 1, start + 100000, posfile, idsfile, ORIG % chrom, BT, out))
    sh("%s index -f %s" % (BT, out))

def cluster(bcf, out, threads=1):
    sh("nice -n 19 ionice -c3 %s cluster -g %s -f 16 -s 8 -p 0.1 --max-clusters 32 "
       "--min-mac 100 --min-freq 0.000001 --medians -t %d -o %s > %s.stdout 2>&1"
       % (HAPLA, bcf, threads, out, out))

def predict(bcf, ref, out, threads=1):
    sh("nice -n 19 ionice -c3 %s predict -g %s -r %s -t %d -o %s > %s.stdout 2>&1"
       % (HAPLA, bcf, ref, threads, out, out))

def main(argv=None):
    ap = argparse.ArgumentParser("hcf.scluster")
    ap.add_argument("--tiles", required=True)
    ap.add_argument("--vardir", required=True)
    ap.add_argument("--split", required=True)
    ap.add_argument("--private", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--hapla", required=True)
    ap.add_argument("--subsamples", type=int, default=10)
    ap.add_argument("--fraction", type=float, default=0.8)
    args = ap.parse_args(argv)
    global HAPLA
    HAPLA = args.hapla

    os.makedirs(args.work, exist_ok=True)
    sp = pd.read_csv(args.split, sep="\t")
    idsA = sp.loc[sp["split"] == "A", "sample_id"].tolist()
    tl = pd.read_csv(args.tiles, sep="\t", dtype={"chrom": str})

    rows = []
    for i in range(len(tl)):
        ch, ts = str(tl["chrom"].iloc[i]), int(tl["tile_start"].iloc[i])
        tile_id = "chr%s:%d" % (ch, ts)
        wd = os.path.join(args.work, "t%s_%d" % (ch, ts))
        os.makedirs(wd, exist_ok=True)
        t0 = time.time()
        try:
            v = pd.read_csv(os.path.join(args.vardir, "chr%s.primary.tsv.gz" % ch), sep="\t")
            v = v[(v["pos"] >= ts) & (v["pos"] < ts + 100000)]
            pf = os.path.join(wd, "pos.tsv")
            v[["pos"]].assign(c=ch)[["c", "pos"]].to_csv(pf, sep="\t", header=False, index=False)
            make_bcf(ch, ts, pf, os.path.join(args.private, "idsA.txt"), os.path.join(wd, "A.bcf"))
            make_bcf(ch, ts, pf, os.path.join(args.private, "idsBC.txt"), os.path.join(wd, "BC.bcf"))
            cluster(os.path.join(wd, "A.bcf"), os.path.join(wd, "A"))
            predict(os.path.join(wd, "BC.bcf"), os.path.join(wd, "A"), os.path.join(wd, "BC"))
            fullA = read_bca(os.path.join(wd, "A.bca"), 2 * len(idsA))
            idsA_order = read_ids(os.path.join(wd, "A.ids"))
            posA = dict((s, j) for j, s in enumerate(idsA_order))
            bc = read_bca(os.path.join(wd, "BC.bca"), 2 * (len(sp) - len(idsA)))
            n_win = fullA.shape[0]
            k_full = [int(fullA[w].max()) + 1 for w in range(n_win)]
            k_bc = [int(bc[w].max()) + 1 for w in range(n_win)]

            rng = np.random.RandomState(20260922)
            for rep in range(args.subsamples):
                sub = sorted(rng.choice(len(idsA_order),
                                        size=int(args.fraction * len(idsA_order)), replace=False))
                subf = os.path.join(wd, "sub.txt")
                with open(subf, "w") as fh:
                    fh.write("\n".join(idsA_order[j] for j in sub) + "\n")
                sh("%s view -S %s -Ob -o %s/sub.bcf %s/A.bcf && %s index -f %s/sub.bcf"
                   % (BT, subf, wd, wd, BT, wd))
                cluster(os.path.join(wd, "sub.bcf"), os.path.join(wd, "sub"))
                subA = read_bca(os.path.join(wd, "sub.bca"), 2 * len(sub))
                hidx = np.concatenate([2 * np.asarray(sub), 2 * np.asarray(sub) + 1])
                hidx.sort()
                for w in range(min(n_win, subA.shape[0])):
                    a = fullA[w][hidx]
                    b = subA[w]
                    ari = float(adjusted_rand_score(a, b))
                    agree = float(np.mean(a == b))
                    rows.append(dict(
                        run_id="HCF-20260922-v1/run_20260922T052403Z",
                        protocol_version="HCF-20260922-v1", evidence_class="[실측-신규]",
                        scope_id="A", build="GRCh37", chrom=ch, tile_id=tile_id,
                        subwindow_id="%s:hapla_w%03d" % (tile_id, w),
                        representation_id="S_CLUSTER_hapla_v0.62.0",
                        perturbation="A_subsample_%.0fpct_rep%d" % (args.fraction * 100, rep),
                        split_id="A", ari=ari, assignment_agreement_raw=agree,
                        n_haplotypes_evaluated=int(len(a)),
                        n_clusters_full=int(a.max()) + 1, n_clusters_perturbed=int(b.max()) + 1,
                        phenotype_used=False,
                        note="label-id agreement is not permutation-invariant; ARI is the primary",
                        status="OK", limitation=""))
                os.remove(os.path.join(wd, "sub.bcf"))
            for w in range(n_win):
                rows.append(dict(
                    run_id="HCF-20260922-v1/run_20260922T052403Z",
                    protocol_version="HCF-20260922-v1", evidence_class="[실측-신규]",
                    scope_id="B_C", build="GRCh37", chrom=ch, tile_id=tile_id,
                    subwindow_id="%s:hapla_w%03d" % (tile_id, w),
                    representation_id="S_CLUSTER_hapla_v0.62.0",
                    perturbation="inductive_projection_A_medians_to_B_C",
                    split_id="B+C", ari=float("nan"), assignment_agreement_raw=float("nan"),
                    n_haplotypes_evaluated=int(bc.shape[1]),
                    n_clusters_full=k_full[w], n_clusters_perturbed=k_bc[w],
                    phenotype_used=False,
                    note="A medians projected onto B/C via hapla predict; no refit on B/C",
                    status="OK", limitation=""))
            sys.stderr.write("%s done %.1fs windows=%d\n" % (tile_id, time.time() - t0, n_win))
        except Exception as exc:
            rows.append(dict(run_id="HCF-20260922-v1/run_20260922T052403Z", chrom=ch,
                             tile_id=tile_id, representation_id="S_CLUSTER_hapla_v0.62.0",
                             status="TILE_FAILED",
                             limitation="%s: %s" % (type(exc).__name__, exc)))
            sys.stderr.write("%s FAILED %s\n" % (tile_id, exc))
        sys.stderr.flush()
        for f in ("A.bcf", "A.bcf.csi", "BC.bcf", "BC.bcf.csi", "A.bca", "BC.bca",
                  "A.ids", "BC.ids", "sub.bca", "sub.ids"):
            p = os.path.join(wd, f)
            if os.path.exists(p):
                os.remove(p)
    pd.DataFrame(rows).to_csv(args.out, index=False)
    print(json.dumps({"rows": len(rows), "out": args.out}))
    return 0

if __name__ == "__main__":
    sys.exit(main())
