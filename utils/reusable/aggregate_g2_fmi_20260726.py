#!/usr/bin/env python3
import os as _cfg_os
import math as _cfg_math

def _config_number(name, cast=int, positive=False):
    value = cast(_cfg_os.environ[name])
    if not _cfg_math.isfinite(value) or (positive and value <= 0):
        raise ValueError(name + " must be finite" + (" and positive" if positive else ""))
    return value

import math, os, sys, glob

EXPECTED_CHROMS = [str(c) for c in range(1, 23)]
EXPECTED_TOTAL_RECORDS = _config_number("EXPECTED_TOTAL_RECORDS", int, True)

def read_kv(path):
    kv = {}
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) == 2:
                kv[p[0]] = p[1]
    return kv

def main(rundir, m_max=None):
    outs = sorted(glob.glob(os.path.join(rundir, "out", "chr*.tsv")))
    outs = [p for p in outs if ".building." not in p]
    per = {}
    for p in outs:
        kv = read_kv(p)
        c = kv.get("chrom")
        if c:
            per[c] = kv

    missing = [c for c in EXPECTED_CHROMS if c not in per]
    floors = sorted({k.split("_")[1] for k in next(iter(per.values())) if k.startswith("floor_")},
                    key=float, reverse=True) if per else []

    tot_rec = sum(int(kv["total_records"]) for kv in per.values())
    lines = []
    lines.append(f"chromosomes_present\t{len(per)}/22")
    lines.append(f"missing_chromosomes\t{','.join(missing) if missing else '-'}")
    lines.append(f"total_records\t{tot_rec}")
    lines.append(f"total_records_expected\t{EXPECTED_TOTAL_RECORDS}")
    lines.append(f"total_records_match\t{tot_rec == EXPECTED_TOTAL_RECORDS}")
    lines.append(f"samples\t{next(iter(per.values()))['samples'] if per else '-'}")

    if missing:
        lines.append("STATUS\tINCOMPLETE — genome-wide FMI withheld")
        print("\n".join(lines))
        return 1

    lines.append("STATUS\tCOMPLETE")
    lines.append("")
    for band, pref in (("all_variants", "floor_"), ("rare_minor_af_lt_0.001", "rare_floor_")):
        lines.append(f"# band = {band}")
        lines.append("floor\tkept\tmaf_sum\tinfo_sum\tinformation_fraction\tFMI\tM_raw\tM_min10")
        for f in floors:
            kept = sum(int(per[c][f"{pref}{f}_kept"]) for c in EXPECTED_CHROMS)
            ms = math.fsum(float(per[c].get(f"{pref}{f}_maf_sum", 0.0)) for c in EXPECTED_CHROMS)
            isum = math.fsum(float(per[c].get(f"{pref}{f}_info_sum", 0.0)) for c in EXPECTED_CHROMS)
            if ms <= 0:
                lines.append(f"{f}\t{kept}\t0\t0\tNA\tNA\tNA\tNA")
                continue
            frac = isum / ms
            fmi = 1 - frac
            mraw = math.ceil(100 * fmi - 1e-12)
            m = max(10, mraw)
            if m_max:
                m = min(m, m_max)
            lines.append(f"{f}\t{kept}\t{ms:.17g}\t{isum:.17g}\t{frac:.17g}\t{fmi:.17g}\t{mraw}\t{m}")
        lines.append("")

    for band, pref in (("all_variants", "floor_"), ("rare_minor_af_lt_0.001", "rare_floor_")):
        lines.append(f"# per-chromosome M — band = {band}")
        lines.append("chrom\t" + "\t".join(f"FMI_{f}/M_{f}" for f in floors))
        for c in EXPECTED_CHROMS:
            cells = []
            for f in floors:
                fmi_s = per[c].get(f"{pref}{f}_FMI")
                m_s = per[c].get(f"{pref}{f}_White_M_min10")
                kept_s = per[c].get(f"{pref}{f}_kept")
                if fmi_s is None or m_s is None or kept_s is None:
                    cells.append("NA/NA")
                elif int(kept_s) == 0:
                    cells.append("empty")
                else:
                    cells.append(f"{float(fmi_s):.4f}/{m_s}")
            lines.append(f"chr{c}\t" + "\t".join(cells))
        lines.append("")

    print("\n".join(lines))
    return 0

if __name__ == "__main__":
    rd = sys.argv[1]
    mm = int(sys.argv[2]) if len(sys.argv) > 2 else None
    sys.exit(main(rd, mm))
