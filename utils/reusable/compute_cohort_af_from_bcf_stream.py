#!/usr/bin/env python3
from __future__ import annotations

import math
import sys

GP_THRESHOLD = 0.90
EPS = 1e-8

def gt_alt_count(gt: str) -> int | None:
    if not gt or gt == ".":
        return None
    alleles = gt.replace("|", "/").split("/")
    if len(alleles) != 2:
        return None
    total = 0
    for allele in alleles:
        if allele in {"", "."}:
            return None
        try:
            total += 1 if int(allele) > 0 else 0
        except ValueError:
            return None
    return total

def gp90_alt_count(gp: str) -> int | None:
    if not gp or gp == ".":
        return None
    try:
        vals = [float(x) for x in gp.split(",")]
    except ValueError:
        return None
    if len(vals) != 3:
        return None
    best = max(range(3), key=lambda i: vals[i])
    if vals[best] + 1e-12 < GP_THRESHOLD:
        return None
    return best

def parse_ds(ds: str) -> float | None:
    if not ds or ds == ".":
        return None
    try:
        return float(ds)
    except ValueError:
        return None

def af_fields(ac: float, an: int) -> tuple[str, str, str]:
    if an <= 0:
        return "NA", "NA", "NA"
    af = min(max(ac / an, 0.0), 1.0)
    maf = min(af, 1.0 - af)
    maf_mb = (af * (1.0 - af) + EPS) ** (-0.5)
    return f"{af:.12g}", f"{maf:.12g}", f"{maf_mb:.12g}"

def parse_cell(cell: str) -> tuple[str, str, str]:
    parts = cell.split(":")
    if len(parts) != 3:
        return ".", ".", "."
    return parts[0], parts[1], parts[2]

def main() -> int:
    print(
        "\t".join(
            [
                "chrom",
                "pos",
                "id",
                "ref",
                "alt",
                "AC_GT",
                "AN_GT",
                "AF_GT",
                "MAF_GT",
                "MAF_MB_GT",
                "AC_GP90",
                "AN_GP90",
                "AF",
                "MAF",
                "MAF_MB",
                "AC_DS",
                "AN_DS",
                "AF_DS",
                "MAF_DS",
                "MAF_MB_DS",
                "missing_GP90",
            ]
        )
    )
    for raw in sys.stdin.buffer:
        line = raw.decode("ascii", errors="replace").rstrip("\n")
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 5:
            continue
        chrom, pos, vid, ref, alt = fields[:5]
        gt_ac = 0
        gt_an = 0
        gp_ac = 0
        gp_an = 0
        ds_ac = 0.0
        ds_an = 0
        missing_gp = 0
        for cell in fields[5:]:
            gt, gp, ds = parse_cell(cell)
            gt_count = gt_alt_count(gt)
            gp_count = gp90_alt_count(gp)
            ds_val = parse_ds(ds)
            if gt_count is not None:
                gt_ac += gt_count
                gt_an += 2
            if gp_count is not None:
                gp_ac += gp_count
                gp_an += 2
            else:
                missing_gp += 1
            if ds_val is not None and not math.isnan(ds_val):
                ds_ac += ds_val
                ds_an += 2
        gt_af, gt_maf, gt_maf_mb = af_fields(gt_ac, gt_an)
        gp_af, gp_maf, gp_maf_mb = af_fields(gp_ac, gp_an)
        ds_af, ds_maf, ds_maf_mb = af_fields(ds_ac, ds_an)
        print(
            "\t".join(
                map(
                    str,
                    [
                        chrom,
                        pos,
                        vid,
                        ref,
                        alt,
                        gt_ac,
                        gt_an,
                        gt_af,
                        gt_maf,
                        gt_maf_mb,
                        gp_ac,
                        gp_an,
                        gp_af,
                        gp_maf,
                        gp_maf_mb,
                        f"{ds_ac:.12g}",
                        ds_an,
                        ds_af,
                        ds_maf,
                        ds_maf_mb,
                        missing_gp,
                    ],
                )
            )
        )
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
