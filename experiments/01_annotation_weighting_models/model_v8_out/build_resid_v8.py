#!/usr/bin/env python3
import os as _cfg_os
import re as _cfg_re
import math as _cfg_math

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)

def _config_number(name, cast=int, positive=False):
    value = cast(_cfg_os.environ[name])
    if not _cfg_math.isfinite(value) or (positive and value <= 0):
        raise ValueError(name + " must be finite" + (" and positive" if positive else ""))
    return value

import os as _os
N_SAMPLES = _config_number("N_SAMPLES", int, True)
import argparse
import collections
import csv
import json
import math
import os
from pathlib import Path
import sys

sys.dont_write_bytecode = True
from l1_train_v8 import (OWNER, PREREG, NUM_EPS, atomic_json, digest, log,
                         norm_iid, owned_path, read_traits, require,
                         resource_guard, safety_args, sha_file)

BASE = Path(_config_path("${PHENO_DIR}"))
EIG = Path(_config_path("${EIGENVEC_FILE}"))
COHORT_CONFIG = {"CT": "COHORT_A", "NC": "COHORT_B", "AS": "COHORT_C"}
TABLES = {cohort: {kind: prefix + "_" + kind + "_FILE"
                   for kind in ("SUBJECT", "LAB", "ANTHROPOMETRY")}
          for cohort, prefix in COHORT_CONFIG.items()}
MISS = {"", "NA", "NaN", "nan", ".", "-9", "99999", "77777", "88888", "9999", "999"}

def numeric(value):
    value = value.strip()
    if value in MISS:
        return None
    try:
        number = float(value)
    except ValueError:
        return None
    return number if math.isfinite(number) else None

class Sources:
    def __init__(self, args, sample_index):
        self.args = args
        self.sample_index = sample_index
        self.cache = {}
        self.audit = {}

    def table(self, cohort, kind):
        require(kind in TABLES[cohort], "unrecognized preregistered table")
        return Path(_config_path("${" + TABLES[cohort][kind] + "}"))

    def read_columns(self, path, columns):
        path = Path(path).resolve()
        key = (str(path), tuple(columns))
        if key in self.cache:
            return self.cache[key]
        for (cached_path, cached_columns), cached in self.cache.items():
            if cached_path == str(path) and set(columns) <= set(cached_columns):
                return cached
        result = {}
        seen = set()
        duplicate = 0
        malformed = 0
        total = 0
        matched = 0
        nonnumeric = collections.Counter()
        with open(path, newline="", errors="strict") as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = next(reader)
            require(len(set(header)) == len(header), "duplicate source header columns")
            require(set(columns) <= set(header), "preregistered source column missing")
            indexes = [header.index(c) for c in columns]
            for fields in reader:
                total += 1
                if len(fields) != len(header):
                    malformed += 1
                    continue
                iid = norm_iid(fields[0])
                if iid in seen:
                    duplicate += 1
                    continue
                seen.add(iid)
                if iid not in self.sample_index:
                    continue
                matched += 1
                values = {}
                for col, ix in zip(columns, indexes):
                    value = numeric(fields[ix])
                    if value is None:
                        nonnumeric[col] += 1
                    else:
                        values[col] = value
                result[iid] = values
        require(duplicate == 0, "duplicate normalized source IID; explicit key resolution required")
        require(malformed == 0, "malformed source rows")
        require(matched > 0, "source-to-genotype IID join matched zero rows")
        self.audit[digest(key)] = dict(file=str(path), columns=columns, rows=total,
                                      key="normalized first-column IID", duplicate_keys=duplicate,
                                      malformed_rows=malformed, matched_genotype_rows=matched,
                                      missing_or_nonnumeric=dict(nonnumeric), sha256=sha_file(path))
        self.cache[key] = result
        return result

def load_samples(path):
    with np.load(path, allow_pickle=True) as z:
        samples = z["samples"].astype(str).tolist()
    normalized = list(map(norm_iid, samples))
    require(len(samples) == N_SAMPLES, "genotype sample count differs from v72")
    require(len(normalized) == len(set(normalized)), "duplicate normalized genotype IID")
    return samples, {iid: i for i, iid in enumerate(normalized)}

def load_pcs(path, sample_index):
    pcs = {}
    seen = set()
    with open(path) as fh:
        header = fh.readline().strip().split()
        require({"IID"} | {f"PC{i}" for i in range(1, 11)} <= set(header), "eigenvec header lacks PC1..10")
        require(len(header) == len(set(header)), "duplicate eigenvec header")
        iid_col = header.index("IID")
        indexes = [header.index(f"PC{i}") for i in range(1, 11)]
        for line in fh:
            fields = line.split()
            require(len(fields) == len(header), "eigenvec row width")
            iid = norm_iid(fields[iid_col])
            require(iid not in seen, "duplicate normalized eigenvec IID")
            seen.add(iid)
            if iid in sample_index:
                values = [numeric(fields[i]) for i in indexes]
                if all(v is not None for v in values):
                    pcs[iid] = values
    require(len(pcs) == len(sample_index), "PC1..10 not complete for genotype sample set")
    return pcs

def load_covariates(sources, pcs):
    cov = {}
    cohort_counts = collections.Counter()
    for cohort in ("CT", "NC", "AS"):
        columns = [os.environ[COHORT_CONFIG[cohort] + "_AGE_COLUMN"],
                   os.environ[COHORT_CONFIG[cohort] + "_SEX_COLUMN"]]
        data = sources.read_columns(sources.table(cohort, "SUBJECT"), columns)
        for iid, fields in data.items():
            require(iid not in cov, "same normalized IID occurs in multiple cohorts")
            age = fields.get(columns[0])
            sex = fields.get(columns[1])
            require(sex is None or sex in (1, 2), "unexpected sex coding")
            cov[iid] = dict(cohort=cohort, age=age, sex=None if sex is None else int(sex == 1),
                            pcs=pcs[iid])
            cohort_counts[cohort] += 1
    require(set(cohort_counts) == {"CT", "NC", "AS"}, "cohort assignment missing")
    return cov, dict(cohort_counts)

def source_plan(row, schema, sources):
    trait = row["trait"]
    if schema == "file_column":
        path = Path(row["file"])
        if not path.is_absolute():
            path = Path(sources.args.traits).resolve().parent / path
        require(bool(row["col"]), "empty generic trait column")
        return [(None, path, [row["col"]], "direct")]
    if row["type"] == "binary":
        require(trait in ("htn", "dm", "lip") and row["CT_col"] == "(pheno_cur)",
                "wide binary source must explicitly name pheno_cur")
        return [(None, Path(sources.args.binary_dir) / f"{trait}_cur.tsv", ["y"], "binary_current")]
    plan = []
    for cohort in ("CT", "NC", "AS"):
        col = row[cohort + "_col"]
        if not col:
            continue
        path = sources.table(cohort, row["table"])
        if col == "(WEIGHT/HEIGHT²)":
            require(trait == "bmi" and cohort == "NC", "unsupported derived-column expression")
            plan.append((cohort, path, [os.environ["COHORT_B_WEIGHT_COLUMN"], os.environ["COHORT_B_HEIGHT_COLUMN"]], "bmi"))
        else:
            plan.append((cohort, path, [col], "direct"))
    require(bool(plan), "trait has no declared sources")
    return plan

def read_scales(path, traits):
    scales = {}
    if not path:
        return scales
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        require({"trait", "cohort", "scale"} <= set(reader.fieldnames or []), "unit-scale TSV schema")
        for row in reader:
            key = (row["trait"], row["cohort"])
            require(key not in scales and key[0] in traits and key[1] in ("CT", "NC", "AS"), "unit-scale key")
            scale = float(row["scale"])
            require(math.isfinite(scale) and scale > 0, "unit scale must be positive finite")
            scales[key] = scale
    return scales

def collect_trait(row, plan, sources, cov, scales):
    values = {}
    excluded = collections.Counter()
    by_cohort = collections.Counter()
    for declared_cohort, path, columns, operation in plan:
        records = sources.read_columns(path, columns)
        for iid, fields in records.items():
            if iid not in cov:
                excluded["no_cohort_assignment"] += 1
                continue
            cohort = cov[iid]["cohort"]
            require(declared_cohort is None or declared_cohort == cohort, "source/participant cohort mismatch")
            if not all(c in fields for c in columns):
                excluded["missing_source_value"] += 1
                continue
            if operation == "bmi":
                height = fields[os.environ["COHORT_B_HEIGHT_COLUMN"]]
                weight = fields[os.environ["COHORT_B_WEIGHT_COLUMN"]]
                if height <= 0 or weight <= 0:
                    excluded["nonpositive_bmi_inputs"] += 1
                    continue
                height_m = height / 100 if sources.args.height_unit == "cm" else height
                value = weight / (height_m ** 2)
            else:
                value = fields[columns[0]]
            if row["type"] == "binary":
                require(value in (0, 1), "binary input must preserve existing 0/1 definition")
            else:
                value *= scales.get((row["trait"], cohort), 1.0)
                if row["trait"] == "tchl" and not 0 < value < 600:
                    excluded["tchl_v3_range"] += 1
                    continue
            if row["transform"] == "rint(log)":
                if value <= 0:
                    excluded["nonpositive_before_log"] += 1
                    continue
                value = math.log(value)
            require(iid not in values, "duplicate phenotype IID across declared sources")
            require(math.isfinite(value), "nonfinite derived phenotype")
            values[iid] = value
            by_cohort[cohort] += 1
    require(bool(values), "no usable trait values after source gates")
    return values, dict(excluded= dict(excluded), usable_before_covariates=len(values),
                        by_cohort_before_covariates=dict(by_cohort))

def fit_residual(row, values, cov, sample_index, min_n, min_cov_retained):
    keys = [iid for iid in values if cov[iid]["age"] is not None and cov[iid]["sex"] is not None]
    keys.sort(key=sample_index.__getitem__)
    require(len(keys) >= min_n, "insufficient complete phenotype/covariate rows")
    require(len(keys) / len(values) >= min_cov_retained, "covariate join retained fraction below preregistered gate")
    y = np.array([values[k] for k in keys], dtype=np.float64)
    present = [c for c in ("CT", "NC", "AS") if any(cov[k]["cohort"] == c for k in keys)]
    baseline = "AS" if "AS" in present else present[0]
    dummy_names = [c for c in present if c != baseline]
    X = np.array([[1.0, cov[k]["age"], cov[k]["sex"]] +
                  [float(cov[k]["cohort"] == c) for c in dummy_names] + cov[k]["pcs"] for k in keys])
    names = ["intercept", "age", "sex_male"] + dummy_names + [f"PC{i}" for i in range(1, 11)]
    keep = [0]
    dropped = []
    for j in range(1, X.shape[1]):
        sd = X[:, j].std(ddof=0)
        if sd <= NUM_EPS:
            dropped.append(names[j])
        else:
            X[:, j] = (X[:, j] - X[:, j].mean()) / sd
            keep.append(j)
    X = X[:, keep]
    require(len(y) > X.shape[1], "nonpositive regression residual degrees of freedom")
    coef, _, rank, singular = np.linalg.lstsq(X, y, rcond=None)
    require(rank == X.shape[1], "covariate rank deficiency beyond constant columns")
    residual = y - X @ coef
    require(residual.std(ddof=0) > NUM_EPS, "zero regression residual variance")
    orthogonality = float(np.max(np.abs(X.T @ residual)) / (len(y) * residual.std()))
    ranks = rankdata(residual, method="average")
    transformed = ndtri((ranks - 0.375) / (len(ranks) + 0.25)).astype(np.float32)
    require(np.isfinite(transformed).all() and transformed.std() > NUM_EPS, "RINT residual gate")
    byc = collections.Counter(cov[k]["cohort"] for k in keys)
    return keys, transformed, dict(n=len(keys), loss_after_covariate_join=len(values) - len(keys),
                                   covariate_retained_fraction=len(keys) / len(values),
                                   min_cov_retained=min_cov_retained,
                                   by_cohort=dict(byc), covariates=[names[j] for j in keep],
                                   cohort_reference=baseline,
                                   dropped_constant_covariates=dropped, rank=int(rank),
                                   residual_df=len(keys) - int(rank),
                                   design_condition=float(singular[0] / singular[-1]),
                                   ols_orthogonality=orthogonality,
                                   cases=int(y.sum()) if row["type"] == "binary" else None,
                                   resid_kind="cov_ols_age_sex_cohort_PC1to10_then_RINT_supervision_only",
                                   order="optional prereg log(y) -> OLS residual -> pooled average-rank RINT",
                                   rint_formula="Phi^-1((average_rank-3/8)/(n+1/4))",
                                   rint_reprojected=False, supervision_only=True)

def write_residual(path, keys, residual, samples, sample_index):
    path = owned_path(path)
    tmp = Path(str(path) + ".tmp")
    fd = os.open(tmp, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
    with os.fdopen(fd, "w") as fh:
        fh.write("IID\tresid\n")
        for key, value in zip(keys, residual):
            fh.write(samples[sample_index[key]] + "\t" + format(float(value), ".9g") + "\n")
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)
    os.chmod(path, 0o600)

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--traits", default=str(PREREG))
    ap.add_argument("--base", default=str(BASE), help="legacy option; source files use COHORT_A/B/C_*_FILE settings")
    ap.add_argument("--eig", default=str(EIG))
    ap.add_argument("--sample-ds", default=str(OWNER / "work/ref15/annot/ds/chr22.ds.npz"))
    ap.add_argument("--binary-dir", default=str(OWNER / "work/ref/pheno_cur"))
    ap.add_argument("--out", default=str(OWNER / "work/ref/annot/resid"))
    ap.add_argument("--unit-scale-tsv", help="optional preregistered trait,cohort,scale table")
    ap.add_argument("--height-unit", choices=["cm", "m"], default="cm")
    ap.add_argument("--min-n", type=int, default=1000)
    ap.add_argument("--min-cov-retained", type=float, default=0.9)
    ap.add_argument("--overwrite", action="store_true")
    safety_args(ap)
    ap.set_defaults(device="cpu")
    args = ap.parse_args()
    require(args.device == "cpu", "residual builder requires CPU only")
    require(args.min_n > 0, "positive minimum sample count required")
    require(0 < args.min_cov_retained <= 1, "covariate retention fraction must be in (0,1]")
    rows, schema = read_traits(args.traits)
    require(all(r["type"] != "binary" or r["transform"] in ("-", "none") for r in rows),
            "binary source definitions cannot be log-transformed")
    lock, guard = resource_guard(args)
    global np, rankdata, ndtri
    import numpy as np
    from scipy.stats import rankdata
    from scipy.special import ndtri
    root = Path(args.out)
    done = root / "build_resid_v8.done"
    summary_path = root / "build_resid_v8.summary.json"
    outputs = [root / (row["trait"] + ".resid.tsv") for row in rows]
    require(args.overwrite or not any(p.exists() for p in outputs + [done, summary_path]),
            "residual outputs exist; use a new output directory or explicit --overwrite")
    if done.exists():
        done.unlink()
    samples, sample_index = load_samples(args.sample_ds)
    pcs = load_pcs(args.eig, sample_index)
    sources = Sources(args, sample_index)
    cov, cohort_counts = load_covariates(sources, pcs)
    plans = {r["trait"]: source_plan(r, schema, sources) for r in rows}
    union_columns = {}
    for plan in plans.values():
        for _, path, columns, _ in plan:
            union_columns.setdefault(path, set()).update(columns)
    for path, columns in union_columns.items():
        sources.read_columns(path, sorted(columns))
    scales = read_scales(args.unit_scale_tsv, [r["trait"] for r in rows])
    summary = dict(version="v8", code_sha256=sha_file(__file__),
                   trainer_helper_sha256=sha_file(Path(__file__).with_name("l1_train_v8.py")),
                   traits_sha256=sha_file(args.traits), trait_schema=schema, preregistered_traits=rows,
                   phenotype_selection="exact list; missing sources or failed traits stop the build",
                   sample_order_sha256=digest(samples), sample_n=len(samples),
                   eigenvec_sha256=sha_file(args.eig), eigenvec_pc_count=10,
                   cohort_counts=cohort_counts, missing_codes=sorted(MISS),
                   units_status="source names/counts alone do not verify cross-cohort units",
                   unit_scale_sha256=sha_file(args.unit_scale_tsv) if args.unit_scale_tsv else None,
                   nc_bmi_height_unit=args.height_unit,
                   residual_policy="supervision only; OLS residual then pooled RINT, including binary cov mode",
                   raw_phenotype_or_residual_values_in_summary=False, resource_guard=guard, traits={})
    full_mask = np.zeros((len(rows), len(samples)), bool)
    for it, row in enumerate(rows):
        trait = row["trait"]
        values, info = collect_trait(row, plans[trait], sources, cov, scales)
        keys, residual, fit_info = fit_residual(row, values, cov, sample_index, args.min_n, args.min_cov_retained)
        path = root / f"{trait}.resid.tsv"
        write_residual(path, keys, residual, samples, sample_index)
        full_mask[it, [sample_index[k] for k in keys]] = True
        summary["traits"][trait] = dict(info, **fit_info, residual_sha256=sha_file(path),
                                         source_plan=[dict(cohort=c, file=str(p), columns=cols, operation=op)
                                                      for c, p, cols, op in plans[trait]],
                                         unit_scales={c: scales.get((trait, c), 1.0) for c in ("CT", "NC", "AS")})
        atomic_json(root / "build_resid_v8.summary.partial.json", summary)
        log(f"built trait {trait}: n={len(keys)}, source/covariate key gates passed")
    summary["source_key_audit"] = sources.audit
    summary["complete_case_n"] = int(full_mask.all(axis=0).sum())
    summary["union_n"] = int(full_mask.any(axis=0).sum())
    summary["pairwise_overlap"] = (full_mask.astype(np.int64) @ full_mask.astype(np.int64).T).tolist()
    summary["intersection_loss"] = {r["trait"]: int(full_mask[i].sum()) - summary["complete_case_n"]
                                       for i, r in enumerate(rows)}
    atomic_json(summary_path, summary)
    atomic_json(done, dict(status="complete", summary_sha256=sha_file(summary_path),
                           traits_sha256=summary["traits_sha256"], n_traits=len(rows)))
    log("RESULT " + json.dumps(dict(summary=str(summary_path), n_traits=len(rows))))
    print("L1_DONE residual_build", flush=True)
    lock.close()

if __name__ == "__main__":
    try:
        main()
    except (ValueError, KeyError) as error:
        raise SystemExit("GATE FAIL: input parse/schema error (" + type(error).__name__ + ")") from None
