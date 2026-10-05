import os as _config_os
import re as _config_re

def _required(name):
    value = _config_os.environ.get(name)
    if value is None or not value.strip():
        raise RuntimeError(f'Set {name} before running this script')
    return value

def _required_int(name):
    value = int(_required(name))
    if value <= 0:
        raise ValueError(f'{name} must be a positive integer')
    return value

def _configured(value):
    defaults = {'BCFTOOLS': 'bcftools', 'PLINK2': 'plink2', 'PYTHON': 'python3', 'TABIX': 'tabix', 'BGZIP': 'bgzip', 'SAMTOOLS': 'samtools', 'BEDTOOLS': 'bedtools', 'CROSSMAP': 'CrossMap', 'UDOCKER_BIN': 'udocker'}
    def resolve(match):
        name = match.group(1) or match.group(2)
        if name in _config_os.environ:
            return _required(name)
        return defaults.get(name) or _required(name)
    return _config_re.sub(r"\$\{([A-Z][A-Z0-9_]*)\}|\$([A-Z][A-Z0-9_]*)", resolve, value)

def _source_path(base, value):
    return value if _config_os.path.isabs(value) else _config_os.path.join(base, value)
import csv, glob, os, subprocess, collections
B = _configured('${PHENO_DIR}')
BCF = _configured('${BCFTOOLS}')
gt = {s.split('_')[0] for s in subprocess.run(_configured(f'{BCF} query -l ${{PROJECT_ROOT}}/work/ref/band_vcf/chr22.band.vcf.gz'), shell=True, capture_output=True, text=True).stdout.split()}
print('=== columns matching TCHL/CHOL across ALL tables ===')
for p in sorted(glob.glob(f'{B}/**/*.tsv', recursive=True)):
    with open(p, encoding='utf-8', errors='ignore') as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
    hit = [h for h in hdr if 'TCHL' in h.upper() or 'CHOL' in h.upper()]
    if hit:
        print(f'  {os.path.relpath(p, B):<46} {hit}')
print()
for rel, col, idc in [(_configured('${COHORT_A_BIOCHEM_FILE}'), _required('COHORT_A_TCHL_COLUMN'), _required('COHORT_A_ID_COLUMN')), (_configured('${COHORT_B_BIOCHEM_FILE}'), _required('COHORT_B_TCHL_COLUMN'), _required('COHORT_B_ID_COLUMN'))]:
    p = _source_path(B, rel)
    with open(p, encoding='utf-8', errors='ignore') as fh:
        rd = csv.reader(fh, delimiter='\t')
        hdr = next(rd)
        if col not in hdr:
            print(f'  [{rel}] COLUMN {col} ABSENT — header sample: {hdr[:8]}')
            continue
        ci = hdr.index(col)
        ii = hdr.index(idc)
        n = 0
        nm = 0
        ng = 0
        vals = []
        for r in rd:
            if not r or len(r) <= ci:
                continue
            n += 1
            sid = r[ii].strip()
            v = r[ci].strip()
            if v in ('', 'NA', '.'):
                continue
            nm += 1
            if sid in gt:
                ng += 1
                try:
                    vals.append(float(v))
                except ValueError:
                    pass
        import statistics as st
        print(f'  [{rel}] rows {n:,} | non-missing {nm:,} | genotyped+measured {ng:,} | mean {st.mean(vals):.1f} sd {st.pstdev(vals):.1f} min {min(vals):.0f} max {max(vals):.0f}')
