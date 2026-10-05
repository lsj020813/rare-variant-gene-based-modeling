
import gzip
import os
import subprocess
import sys

import numpy as np

TILE_BP = 100000
MAX_INDEL_LEN = 50
R2_PRIMARY_MIN = 0.8
QUARANTINE_FLAGS = ("SwappedAlleles", "ReverseComplementedAlleles")

BCFTOOLS = "bcftools"

TILEINV_COLS = [
    "tile_start",
    "n_keymap",
    "n_quarantine",
    "n_not_biallelic_snv_indel",
    "n_coding_cds",
    "n_splice_sensitive",
    "n_monomorphic",
    "n_lowq_imputed",
    "n_primary",
    "n_primary_typed",
]

def _info_dict(info):
    d = {}
    for part in info.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            d[k] = v
        else:
            d[part] = True
    return d

def scan_chrom(chrom, keymap_dir, annot_npz, out_dir):
    from hcf import annot

    cds, spl = annot.load(annot_npz, str(chrom))
    src = os.path.join(keymap_dir, "chr%s.keyed38.vcf.gz" % chrom)
    cmd = [BCFTOOLS, "query", "-f", "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO\n", src]

    tiles = {}

    def tile_row(t):
        if t not in tiles:
            tiles[t] = [t] + [0] * (len(TILEINV_COLS) - 1)
        return tiles[t]

    os.makedirs(out_dir, exist_ok=True)
    prim_path = os.path.join(out_dir, "chr%s.primary.tsv.gz" % chrom)
    n_rows = 0
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, bufsize=1 << 20)
    buf_pos, buf_ref, buf_alt, buf_af, buf_r2, buf_typed = [], [], [], [], [], []
    BATCH = 200000

    pend_idx = []
    pend_pos38 = []
    pend_end38 = []
    pend_rec = []

    def flush_annot():
        if not pend_rec:
            return
        p = np.asarray(pend_pos38, dtype=np.int64)
        e = np.asarray(pend_end38, dtype=np.int64)
        is_cds = cds.overlaps(p, e)
        is_spl = spl.overlaps(p, e)
        for i, rec in enumerate(pend_rec):
            t, pos19, ref, alt, af, r2, typed = rec
            row = tile_row(t)
            if is_cds[i]:
                row[4] += 1
                continue
            if is_spl[i]:
                row[5] += 1
                continue
            if af <= 0.0 or af >= 1.0:
                row[6] += 1
                continue
            if (not typed) and (not (r2 >= R2_PRIMARY_MIN)):
                row[7] += 1
                continue
            row[8] += 1
            if typed:
                row[9] += 1
            buf_pos.append(pos19)
            buf_ref.append(ref)
            buf_alt.append(alt)
            buf_af.append(af)
            buf_r2.append(r2)
            buf_typed.append(1 if typed else 0)
        del pend_idx[:], pend_pos38[:], pend_end38[:], pend_rec[:]

    with gzip.open(prim_path, "wt") as out:
        out.write("pos\tref\talt\taf\tr2\ttyped\n")
        for raw in proc.stdout:
            line = raw.decode("ascii")
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            n_rows += 1
            chrom38, pos38, vid, ref38, alt38, info = f[0], int(f[1]), f[2], f[3], f[4], f[5]
            kp = vid.split(":")
            if len(kp) != 4:
                continue
            key_chrom = kp[0][3:] if kp[0].startswith("chr") else kp[0]
            pos19 = int(kp[1])
            ref19, alt19 = kp[2], kp[3]
            t = (pos19 // TILE_BP) * TILE_BP
            row = tile_row(t)
            row[1] += 1
            d = _info_dict(info)
            bad = any(flag in d for flag in QUARANTINE_FLAGS)
            c38 = chrom38[3:] if chrom38.startswith("chr") else chrom38
            if c38 != key_chrom or key_chrom != str(chrom):
                bad = True
            if bad:
                row[2] += 1
                continue
            ok_allele = (
                len(ref19) <= MAX_INDEL_LEN
                and len(alt19) <= MAX_INDEL_LEN
                and "," not in alt19
                and "<" not in alt19
                and "*" not in alt19
                and set(ref19 + alt19) <= set("ACGTacgt")
            )
            if not ok_allele:
                row[3] += 1
                continue
            try:
                af = float(d.get("AF", "nan"))
            except ValueError:
                af = float("nan")
            try:
                r2 = float(d.get("R2", "nan"))
            except ValueError:
                r2 = float("nan")
            typed = "TYPED" in d
            pend_pos38.append(pos38)
            pend_end38.append(pos38 + max(len(ref38), 1) - 1)
            pend_rec.append((t, pos19, ref19, alt19, af, r2, typed))
            if len(pend_rec) >= 20000:
                flush_annot()
            if len(buf_pos) >= BATCH:
                for i in range(len(buf_pos)):
                    out.write(
                        "%d\t%s\t%s\t%.6g\t%.6g\t%d\n"
                        % (buf_pos[i], buf_ref[i], buf_alt[i], buf_af[i], buf_r2[i], buf_typed[i])
                    )
                del buf_pos[:], buf_ref[:], buf_alt[:], buf_af[:], buf_r2[:], buf_typed[:]
        flush_annot()
        for i in range(len(buf_pos)):
            out.write(
                "%d\t%s\t%s\t%.6g\t%.6g\t%d\n"
                % (buf_pos[i], buf_ref[i], buf_alt[i], buf_af[i], buf_r2[i], buf_typed[i])
            )
    proc.stdout.close()
    rc = proc.wait()
    if rc != 0:
        raise RuntimeError("bcftools query failed rc=%d for %s" % (rc, src))

    inv_path = os.path.join(out_dir, "chr%s.tileinv.tsv" % chrom)
    with open(inv_path, "w") as fh:
        fh.write("chrom\t" + "\t".join(TILEINV_COLS) + "\n")
        for t in sorted(tiles):
            fh.write("%s\t%s\n" % (chrom, "\t".join(str(v) for v in tiles[t])))
    return {"chrom": str(chrom), "keymap_records": n_rows, "n_tiles": len(tiles),
            "primary_table": prim_path, "tile_inventory": inv_path}

def subwindows(positions, markers=16, step=8, max_gap_bp=20000):
    pos = np.asarray(positions, dtype=np.int64)
    n = pos.size
    if n < markers:
        return []
    brk = np.nonzero(np.diff(pos) > max_gap_bp)[0]
    bounds = [0] + [int(b) + 1 for b in brk] + [n]
    out = []
    for a, b in zip(bounds[:-1], bounds[1:]):
        L = b - a
        if L < markers:
            continue
        s = a
        while s + markers <= b:
            out.append((s, s + markers))
            s += step
    return out
