
import gzip
import sys

import numpy as np

SPLICE_INTRONIC_BP = 8
SPLICE_EXONIC_BP = 3

def _merge(intervals):
    if not intervals:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64)
    arr = np.asarray(sorted(intervals), dtype=np.int64)
    starts = [arr[0, 0]]
    ends = [arr[0, 1]]
    for s, e in arr[1:]:
        if s <= ends[-1] + 1:
            if e > ends[-1]:
                ends[-1] = e
        else:
            starts.append(s)
            ends.append(e)
    return np.asarray(starts, dtype=np.int64), np.asarray(ends, dtype=np.int64)

def build(gtf_path, out_npz):
    cds = {}
    exon = {}
    n_lines = 0
    with gzip.open(gtf_path, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            n_lines += 1
            f = line.split("\t", 9)
            feat = f[2]
            if feat != "CDS" and feat != "exon":
                continue
            attr = f[8]
            if 'transcript_type "protein_coding"' not in attr:
                continue
            chrom = f[0]
            s = int(f[3])
            e = int(f[4])
            d = cds if feat == "CDS" else exon
            d.setdefault(chrom, []).append((s, e))

    out = {}
    chroms = sorted(set(list(cds.keys()) + list(exon.keys())))
    for ch in chroms:
        cs, ce = _merge(cds.get(ch, []))
        out["cds_%s_s" % ch] = cs
        out["cds_%s_e" % ch] = ce
        sp = []
        for s, e in exon.get(ch, []):
            sp.append((s - SPLICE_INTRONIC_BP, s + SPLICE_EXONIC_BP - 1))
            sp.append((e - SPLICE_EXONIC_BP + 1, e + SPLICE_INTRONIC_BP))
        ss, se = _merge(sp)
        out["spl_%s_s" % ch] = ss
        out["spl_%s_e" % ch] = se
    out["_chroms"] = np.asarray(chroms)
    out["_meta"] = np.asarray(
        [
            "gtf=%s" % gtf_path,
            "gtf_data_lines=%d" % n_lines,
            "transcript_type=protein_coding",
            "splice_intronic_bp=%d" % SPLICE_INTRONIC_BP,
            "splice_exonic_bp=%d" % SPLICE_EXONIC_BP,
            "assembly=GRCh38",
            "release=ANNOTATION_RELEASE_NOT_RECORDED",
        ]
    )
    np.savez_compressed(out_npz, **out)
    return {
        "out": out_npz,
        "n_chroms": len(chroms),
        "n_cds_intervals": int(sum(len(out["cds_%s_s" % c]) for c in chroms)),
        "n_splice_intervals": int(sum(len(out["spl_%s_s" % c]) for c in chroms)),
        "gtf_data_lines": n_lines,
    }

class Intervals(object):

    def __init__(self, starts, ends):
        self.s = np.asarray(starts, dtype=np.int64)
        self.e = np.asarray(ends, dtype=np.int64)

    def overlaps(self, pos, pos_end):
        pos = np.asarray(pos, dtype=np.int64)
        pos_end = np.asarray(pos_end, dtype=np.int64)
        if self.s.size == 0:
            return np.zeros(pos.shape, dtype=bool)
        idx = np.searchsorted(self.s, pos_end, side="right") - 1
        ok = idx >= 0
        res = np.zeros(pos.shape, dtype=bool)
        if not ok.any():
            return res
        res[ok] = self.e[idx[ok]] >= pos[ok]
        return res

def load(npz_path, chrom):
    z = np.load(npz_path, allow_pickle=False)
    ch = str(chrom)
    cds = Intervals(z["cds_%s_s" % ch], z["cds_%s_e" % ch]) if ("cds_%s_s" % ch) in z else Intervals([], [])
    spl = Intervals(z["spl_%s_s" % ch], z["spl_%s_e" % ch]) if ("spl_%s_s" % ch) in z else Intervals([], [])
    return cds, spl

if __name__ == "__main__":
    import json

    print(json.dumps(build(sys.argv[1], sys.argv[2]), indent=2))
