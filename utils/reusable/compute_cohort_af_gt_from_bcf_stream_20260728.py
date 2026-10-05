#!/usr/bin/env python3
from __future__ import annotations

import sys

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

def af_fields(ac: int, an: int) -> tuple[str, str, str]:
    if an <= 0:
        return "NA", "NA", "NA"
    af = min(max(ac / an, 0.0), 1.0)
    maf = min(af, 1.0 - af)
    maf_mb = (af * (1.0 - af) + EPS) ** (-0.5)
    return f"{af:.12g}", f"{maf:.12g}", f"{maf_mb:.12g}"

HEADER = [
    "chrom",
    "pos",
    "id",
    "ref",
    "alt",
    "AC_GT",
    "AN_GT",
    "AF",
    "MAF",
    "MAF_MB",
    "missing_GT",
]

def main() -> int:
    out = sys.stdout
    out.write("\t".join(HEADER) + "\n")
    for raw in sys.stdin.buffer:
        line = raw.decode("ascii", errors="replace").rstrip("\n")
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 5:
            continue
        chrom, pos, vid, ref, alt = fields[:5]
        ac = 0
        an = 0
        missing = 0
        for cell in fields[5:]:
            count = gt_alt_count(cell)
            if count is None:
                missing += 1
            else:
                ac += count
                an += 2
        af, maf, maf_mb = af_fields(ac, an)
        out.write(
            "\t".join(
                (chrom, pos, vid, ref, alt, str(ac), str(an), af, maf, maf_mb, str(missing))
            )
            + "\n"
        )
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
