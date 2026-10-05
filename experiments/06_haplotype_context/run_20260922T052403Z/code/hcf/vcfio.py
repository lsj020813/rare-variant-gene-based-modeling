
import os
import subprocess

import numpy as np

BCFTOOLS = "bcftools"
MISSING = -1

class GtReadError(RuntimeError):
    pass

def read_region(vcf_path, region, n_samples, regions_file=None, samples_file=None,
                keep=None, max_variants=None, use_view=False, pass_fds=()):
    fmt = "%POS\t%REF\t%ALT[\t%GT]\n"
    if use_view:
        cmd = [BCFTOOLS, "view", "-H", "-r", region]
        if samples_file:
            cmd += ["-S", samples_file]
        cmd += [vcf_path]
    else:
        cmd = [BCFTOOLS, "query"]
        if regions_file:
            cmd += ["-R", regions_file]
        else:
            cmd += ["-r", region]
        if samples_file:
            cmd += ["-S", samples_file]
        cmd += ["-f", fmt, vcf_path]

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, bufsize=1 << 22, pass_fds=tuple(pass_fds))
    pos_l, ref_l, alt_l, rows = [], [], [], []
    n_lines = 0
    sep_all_pipe = True
    try:
        for raw in proc.stdout:
            n_lines += 1
            if use_view:
                continue
            i1 = raw.index(b"\t")
            i2 = raw.index(b"\t", i1 + 1)
            i3 = raw.index(b"\t", i2 + 1)
            pos = int(raw[:i1])
            ref = raw[i1 + 1:i2].decode("ascii")
            alt = raw[i2 + 1:i3].decode("ascii")
            if keep is not None and (pos, ref, alt) not in keep:
                continue
            body = raw[i3:]
            if body.endswith(b"\n"):
                body = body[:-1]
            if len(body) != 4 * n_samples:
                raise GtReadError(
                    "GT field width mismatch at pos %d: got %d bytes, expected %d "
                    "(non-biallelic or non-3-char GT present)" % (pos, len(body), 4 * n_samples)
                )
            a = np.frombuffer(body, dtype=np.uint8).reshape(n_samples, 4)
            if not np.all(a[:, 2] == 124):
                sep_all_pipe = False
                raise GtReadError("unphased or malformed separator at pos %d" % pos)
            g = a[:, [1, 3]].astype(np.int16) - 48
            g[(g < 0) | (g > 1)] = MISSING
            pos_l.append(pos)
            ref_l.append(ref)
            alt_l.append(alt)
            rows.append(g.astype(np.int8))
            if max_variants is not None and len(rows) >= max_variants:
                break
    finally:
        try:
            proc.stdout.close()
        except Exception:
            pass
        proc.wait()

    if use_view:
        return {"n_lines": n_lines, "pos": np.zeros(0, dtype=np.int64), "ref": [], "alt": [],
                "hap": np.zeros((0, n_samples, 2), dtype=np.int8), "sep_all_pipe": None}

    hap = (np.stack(rows, axis=0) if rows else np.zeros((0, n_samples, 2), dtype=np.int8))
    return {
        "pos": np.asarray(pos_l, dtype=np.int64),
        "ref": ref_l,
        "alt": alt_l,
        "hap": hap,
        "n_lines": n_lines,
        "sep_all_pipe": sep_all_pipe,
    }

def to_person_major(hap):
    return np.ascontiguousarray(np.transpose(hap, (1, 2, 0)))
