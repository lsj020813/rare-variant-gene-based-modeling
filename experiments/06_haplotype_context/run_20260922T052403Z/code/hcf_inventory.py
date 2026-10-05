#!/usr/bin/env python
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import os, sys, json, hashlib, subprocess, time, re, glob
import numpy as np
from cyvcf2 import VCF

OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)
BC = "bcftools"
REF = _config_path("${PROJECT_ROOT}/work/ref")
CHROMS = list(range(1, 23))
N_REC = 500
N_SMP = 200

def sh(cmd, timeout=3000):
    p = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=timeout)
    return p.returncode, p.stdout.decode("utf-8", "replace"), p.stderr.decode("utf-8", "replace")

def log(msg):
    print(time.strftime("%H:%M:%S"), msg, flush=True)

def stat_file(p):
    d = {"path": p, "exists": os.path.exists(p), "is_symlink": os.path.islink(p)}
    if d["is_symlink"]:
        d["realpath"] = os.path.realpath(p)
    if d["exists"]:
        st = os.stat(p)
        d["size_bytes"] = st.st_size
        d["mtime_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(st.st_mtime))
        d["readable"] = os.access(p, os.R_OK)
        d["writable_by_me"] = os.access(p, os.W_OK)
        try:
            import pwd
            d["owner"] = pwd.getpwuid(st.st_uid).pw_name
        except Exception:
            d["owner"] = st.st_uid
        d["mode"] = oct(st.st_mode & 0o777)
    return d

def bgzf_eof_ok(p):
    with open(p, "rb") as f:
        f.seek(-28, 2)
        tail = f.read(28)
    return tail.hex() == "1f8b08040000000000ff0600424302001b0003000000000000000000"

def header_info(p):
    rc, out, err = sh("%s view -h %s" % (BC, p))
    if rc != 0:
        return {"error": err.strip()[:500]}
    lines = out.rstrip("\n").split("\n")
    meta = [l for l in lines if l.startswith("##") and not l.startswith("##bcftools_view")]
    hdr_text = "\n".join(meta)
    samples = lines[-1].split("\t")[9:] if lines[-1].startswith("#CHROM") else []
    d = {
        "fileformat": [l for l in meta if l.startswith("##fileformat")],
        "source": [l for l in meta if l.startswith("##source")],
        "phasing": [l for l in meta if l.startswith("##phasing")],
        "reference": [l for l in meta if l.startswith("##reference")],
        "filedate": [l for l in meta if l.startswith("##filedate")],
        "FORMAT": [l for l in meta if l.startswith("##FORMAT")],
        "INFO": [l for l in meta if l.startswith("##INFO")],
        "FORMAT_ids": [re.search(r"ID=([^,>]+)", l).group(1) for l in meta if l.startswith("##FORMAT")],
        "INFO_ids": [re.search(r"ID=([^,>]+)", l).group(1) for l in meta if l.startswith("##INFO")],
        "n_contig_lines": sum(1 for l in meta if l.startswith("##contig")),
        "contig_ids": [re.search(r"ID=([^,>]+)", l).group(1) for l in meta if l.startswith("##contig")][:5],
        "other_meta_keys": sorted(set(l[2:].split("=")[0] for l in meta if not l.startswith(("##FORMAT","##INFO","##contig","##FILTER")))),
        "header_sha256_excl_bcftools_view": hashlib.sha256(hdr_text.encode()).hexdigest(),
        "n_samples": len(samples),
        "sample_list_md5_newline_joined": hashlib.md5(("\n".join(samples) + "\n").encode()).hexdigest() if samples else None,
        "sample_list_sha256": hashlib.sha256(("\n".join(samples) + "\n").encode()).hexdigest() if samples else None,
        "sample_ids_unique": (len(set(samples)) == len(samples)) if samples else None,
    }
    return d

def index_n(p):
    rc, out, err = sh("%s index -n %s" % (BC, p))
    try:
        return int(out.strip())
    except Exception:
        return {"error": (out + err).strip()[:300]}

def step_orig():
    res = {}
    for c in CHROMS:
        p = "%s/orig_index/chr%d.vcf.gz" % (REF, c)
        log("orig chr%d" % c)
        d = {"file": stat_file(p)}
        if d["file"]["is_symlink"]:
            d["target"] = stat_file(d["file"]["realpath"])
            d["target_tbi"] = stat_file(d["file"]["realpath"] + ".tbi")
        d["csi"] = stat_file(p + ".csi")
        d["index_n_records"] = index_n(p)
        d["bgzf_eof_ok"] = bgzf_eof_ok(p)
        d["header"] = header_info(p)
        try:
            v0 = VCF(p)
            smp_all = list(v0.samples)
            sub = smp_all[:N_SMP]
            v0.close()
            v = VCF(p, samples=sub)
            n = 0; pipe = 0; slash = 0; missing = 0; total_calls = 0
            fmt_orders = {}
            gt_ds_round_agree = 0; gt_gp_argmax_agree = 0; ds_gp_consistent = 0; n_cmp = 0
            gt_hds_like = 0
            maf_bins = {"<0.001":0, "0.001-0.01":0, "0.01-0.05":0, ">=0.05":0}
            typed = 0; r2_lt08 = 0; altmajor = 0; nonbiallelic = 0; indel = 0
            phased_flag_true = 0
            for var in v:
                n += 1
                fo = ":".join(var.FORMAT)
                fmt_orders[fo] = fmt_orders.get(fo, 0) + 1
                if len(var.ALT) != 1:
                    nonbiallelic += 1
                if len(var.REF) != 1 or any(len(a) != 1 for a in var.ALT):
                    indel += 1
                maf = var.INFO.get("MAF"); r2 = var.INFO.get("R2"); af = var.INFO.get("AF")
                if var.INFO.get("TYPED") is not None: typed += 1
                if r2 is not None and r2 < 0.8: r2_lt08 += 1
                if af is not None and af > 0.5: altmajor += 1
                if maf is not None:
                    if maf < 0.001: maf_bins["<0.001"] += 1
                    elif maf < 0.01: maf_bins["0.001-0.01"] += 1
                    elif maf < 0.05: maf_bins["0.01-0.05"] += 1
                    else: maf_bins[">=0.05"] += 1
                gts = var.genotypes
                ds = var.format("DS"); gp = var.format("GP")
                for i, g in enumerate(gts):
                    total_calls += 1
                    if g[0] < 0 or g[1] < 0:
                        missing += 1; continue
                    if g[2]: pipe += 1; phased_flag_true += 1
                    else: slash += 1
                    s = g[0] + g[1]
                    if ds is not None and gp is not None:
                        n_cmp += 1
                        dsv = float(ds[i][0]); gpv = gp[i]
                        if int(round(dsv)) == s: gt_ds_round_agree += 1
                        if int(np.argmax(gpv)) == s: gt_gp_argmax_agree += 1
                        if abs((gpv[1] + 2*gpv[2]) - dsv) <= 0.01: ds_gp_consistent += 1
                if n >= N_REC: break
            v.close()
            d["gt_sample"] = {
                "n_records": n, "n_samples": len(sub), "total_calls": total_calls,
                "phased_pipe_calls": pipe, "unphased_slash_calls": slash, "missing_calls": missing,
                "phased_fraction": (pipe / total_calls) if total_calls else None,
                "data_row_FORMAT_orders": fmt_orders,
                "gt_sum_eq_round_DS_frac": (gt_ds_round_agree / n_cmp) if n_cmp else None,
                "gt_sum_eq_argmax_GP_frac": (gt_gp_argmax_agree / n_cmp) if n_cmp else None,
                "DS_eq_GP1_plus_2GP2_within_0.01_frac": (ds_gp_consistent / n_cmp) if n_cmp else None,
                "records_typed": typed, "records_r2_lt_0.8": r2_lt08, "records_alt_af_gt_0.5": altmajor,
                "records_nonbiallelic": nonbiallelic, "records_indel_like": indel, "maf_bins": maf_bins,
                "note": "first %d records in file order x first %d samples in file order; [실측-신규][표본]" % (N_REC, N_SMP),
            }
        except Exception as e:
            d["gt_sample"] = {"error": repr(e)[:500]}
        res["chr%d" % c] = d
        json.dump(res, open(OUT + "/orig_index_inventory.json", "w"), indent=1)
    return res

def step_derived():
    res = {}
    fam = {"band_vcf": "%s/band_vcf/chr{c}.band.vcf.gz" % REF,
           "band15_vcf": "%s/band15_vcf/chr{c}.band15.vcf.gz" % REF,
           "common05": "%s/common05/chr{c}.maf05.vcf.gz" % REF}
    for name, pat in fam.items():
        res[name] = {}
        for c in CHROMS:
            p = pat.format(c=c)
            log("%s chr%d" % (name, c))
            d = {"file": stat_file(p), "csi": stat_file(p + ".csi"), "done_marker": stat_file(p + ".done")}
            if d["file"]["exists"]:
                d["index_n_records"] = index_n(p)
                d["bgzf_eof_ok"] = bgzf_eof_ok(p)
                h = header_info(p)
                d["header"] = {k: h.get(k) for k in ["FORMAT_ids", "INFO_ids", "n_samples", "sample_list_md5_newline_joined", "sample_list_sha256", "phasing", "source", "header_sha256_excl_bcftools_view"]}
                rc, out, err = sh("%s view -h %s | grep '^##bcftools_' | grep -v 'view -h' | cut -c1-400" % (BC, p))
                d["lineage_bcftools_lines"] = out.strip().split("\n")[:6]
            res[name]["chr%d" % c] = d
        json.dump(res, open(OUT + "/derived_inventory.json", "w"), indent=1)
    return res

def step_lift():
    res = {}
    orig = json.load(open(OUT + "/orig_index_inventory.json")) if os.path.exists(OUT + "/orig_index_inventory.json") else {}
    for c in CHROMS:
        log("lift chr%d" % c)
        p = "%s/lift38_keyed/chr%d.keyed38.vcf.gz" % (REF, c)
        s37 = _config_path("${PROJECT_ROOT}/work/tmp/lift38/extracts/chr%d.s37.vcf.gz") % c
        d = {"keyed38": stat_file(p), "keyed38_tbi": stat_file(p + ".tbi"), "s37_extract": stat_file(s37)}
        if d["keyed38"]["exists"]:
            d["keyed38_index_n"] = index_n(p)
            d["keyed38_bgzf_eof_ok"] = bgzf_eof_ok(p)
            rc, out, err = sh("%s view -H %s | cut -f3,8 | awk -F'\\t' 'BEGIN{n=0;idok=0;sw=0;rc=0;multi=0} {n++; if($1 ~ /^[0-9XY]+:[0-9]+:[ACGTN]+:[ACGTN]+$/) idok++; if($2 ~ /SwappedAlleles/) sw++; if($2 ~ /ReverseComplementedAlleles/) rc++} END{print n\"\\t\"idok\"\\t\"sw\"\\t\"rc}'" % (BC, p))
            try:
                n, idok, sw, rcc = [int(x) for x in out.strip().split("\t")]
                d["keyed38_stream_records"] = n
                d["keyed38_id_matches_hg19_key_pattern"] = idok
                d["keyed38_SwappedAlleles_records"] = sw
                d["keyed38_ReverseComplementedAlleles_records"] = rcc
            except Exception:
                d["keyed38_stream_error"] = (out + err)[:300]
            rc2, out2, err2 = sh("%s view -H %s | cut -f3 | sort | uniq -d | wc -l" % (BC, p))
            d["keyed38_duplicate_hg19_keys"] = int(out2.strip()) if out2.strip().isdigit() else out2[:100]
            h = header_info(p)
            d["keyed38_reference_line"] = h.get("reference")
            d["keyed38_INFO_ids"] = h.get("INFO_ids")
        if d["s37_extract"]["exists"]:
            d["s37_index_n"] = index_n(s37)
            d["s37_bgzf_eof_ok"] = bgzf_eof_ok(s37)
        on = orig.get("chr%d" % c, {}).get("index_n_records")
        d["orig_index_n"] = on
        if isinstance(on, int) and isinstance(d.get("keyed38_stream_records"), int):
            d["keyed38_coverage_vs_orig"] = d["keyed38_stream_records"] / on
        res["chr%d" % c] = d
        json.dump(res, open(OUT + "/lift38_inventory.json", "w"), indent=1)
    return res

def step_altmajor():
    s37 = _config_path("${PROJECT_ROOT}/work/tmp/lift38/extracts/chr22.s37.vcf.gz")
    orig_n = index_n("%s/orig_index/chr22.vcf.gz" % REF)
    d = {"source": s37, "orig_chr22_index_n": orig_n, "s37_index_n": index_n(s37), "s37_bgzf_eof_ok": bgzf_eof_ok(s37)}
    rc, out, err = sh("%s query -f '%%INFO/AF\\t%%INFO/MAF\\t%%INFO/R2\\t%%TYPED\\t%%REF\\t%%ALT\\n' %s" % (BC, s37))
    n=0; alt_major=0; tie=0; maf_mismatch=0; maf_gt_half=0; r2ge08=0; r2ge08_altmajor=0; typed=0; typed_altmajor=0
    bins = {"(0,0.001)":[0,0], "[0.001,0.01)":[0,0], "[0.01,0.05)":[0,0], "[0.05,0.5]":[0,0]}
    snv=0; snv_altmajor=0; indel=0; indel_altmajor=0; af_missing=0
    for line in out.split("\n"):
        if not line: continue
        f = line.split("\t")
        try:
            af = float(f[0]); maf = float(f[1])
        except Exception:
            af_missing += 1; continue
        n += 1
        try: r2 = float(f[2])
        except Exception: r2 = None
        is_typed = (f[3] not in (".", ""))
        is_snv = (len(f[4]) == 1 and len(f[5]) == 1)
        am = af > 0.5
        if af == 0.5: tie += 1
        if am: alt_major += 1
        if abs(min(af, 1-af) - maf) > 0.005: maf_mismatch += 1
        if maf > 0.5: maf_gt_half += 1
        if r2 is not None and r2 >= 0.8:
            r2ge08 += 1
            if am: r2ge08_altmajor += 1
        if is_typed:
            typed += 1
            if am: typed_altmajor += 1
        if is_snv:
            snv += 1; snv_altmajor += am
        else:
            indel += 1; indel_altmajor += am
        key = "(0,0.001)" if maf < 0.001 else "[0.001,0.01)" if maf < 0.01 else "[0.01,0.05)" if maf < 0.05 else "[0.05,0.5]"
        bins[key][0] += 1; bins[key][1] += am
    d.update({"n_records_with_AF": n, "af_missing": af_missing, "alt_major_AF_gt_0.5": alt_major, "alt_major_fraction": alt_major / n if n else None,
              "af_exactly_0.5_ties": tie, "MAF_ne_min(AF,1-AF)_tol0.005": maf_mismatch, "MAF_gt_0.5": maf_gt_half,
              "r2_ge_0.8": r2ge08, "r2_ge_0.8_alt_major": r2ge08_altmajor, "typed": typed, "typed_alt_major": typed_altmajor,
              "snv": snv, "snv_alt_major": snv_altmajor, "indel": indel, "indel_alt_major": indel_altmajor,
              "by_maf_bin_[n, alt_major]": bins,
              "rule": "REF/ALT orientation fixed for haplotype display; minor dosage = 2-DS when INFO/AF>0.5; AF==0.5 tie keeps REF/ALT as minor=ALT (DS as-is); orientation decided from upstream INFO/AF, never from phenotype",
              "grade": "[실측-신규] chr22 전체 사이트(sites-only 추출본 INFO 기준, 원본 index 레코드 수와 대조)"})
    json.dump(d, open(OUT + "/chr22_allele_orientation.json", "w"), indent=1)
    return d

def step_pheno():
    d = {}
    for t in ["tchl_v3.tsv", "lip_v3.tsv", "dm_v3.tsv", "htn_v3.tsv"]:
        p = "%s/pheno_v3/%s" % (REF, t)
        e = stat_file(p)
        if e["exists"]:
            with open(p) as f:
                hdr = f.readline().rstrip("\n").split("\t")
                nrow = sum(1 for _ in f)
            e["columns"] = hdr; e["n_data_rows"] = nrow
        d[t] = e
    d["pheno_v3_summary.json"] = json.load(open("%s/pheno_v3/pheno_v3_summary.json" % REF))
    logs = {}
    for t in ["tchl_v4", "lip_v4"]:
        p = "%s/saige_step1_v4/%s.step1.log" % (REF, t)
        if os.path.exists(p):
            rc, out, err = sh("grep -inE 'covar|samples have|phenoCol|traitType|invNormalize|LOCO|sampleIDCol|plinkFile|bedFile|nThreads|relatednessCutoff|useSparseGRM|sparseGRM|qCovar' %s | cut -c1-300 | head -60" % p)
            logs[t] = out.strip().split("\n")
    d["saige_step1_logs"] = logs
    for t in ["tchl_v4", "lip_v4"]:
        p = "%s/saige_step1_v4/%s.varianceRatio.txt" % (REF, t)
        if os.path.exists(p):
            d[t + ".varianceRatio"] = open(p).read().strip()[:200]
    json.dump(d, open(OUT + "/phenotype_schema.json", "w"), indent=1)
    return d

def step_kinship():
    roots = [_config_path("${PROJECT_ROOT}/work/ref"), _config_path("${PROJECT_ROOT}/work/gate1"), _config_path("${PROJECT_ROOT}/work/gate2"), _config_path("${PROJECT_ROOT}/work/fset"),
             _config_path("${PROJECT_ROOT}/work/ref15"), _config_path("${PROJECT_ROOT}/work/run_band15"), _config_path("${PROJECT_ROOT}/work/model_v10_out"), _config_path("${PROJECT_ROOT}/work/run_meth"),
             _config_path("${PROJECT_ROOT}/array_model_starter_20260728"), _config_path("${COHORT_DATA}"), _config_path("${DATA_ROOT}")]
    pat = r".*(kin|unrel|relat|king|grm|ibd|pihat|pedigree|\.fam|\.kin0|\.rel|sparse).*"
    found = {}
    for r in roots:
        if not os.path.isdir(r):
            found[r] = "NOT_A_DIR_OR_NO_ACCESS"; continue
        rc, out, err = sh("find %s -maxdepth 3 -iregex '%s' -not -path '*/__pycache__/*' 2>/dev/null | head -80" % (r, pat), timeout=600)
        found[r] = [l for l in out.strip().split("\n") if l]
    rc, out, err = sh("ls %s/vr_plink_v3 %s/saige_step1_v4 %s/pheno_census %s/pheno_cur 2>/dev/null" % (REF, REF, REF, REF))
    found["saige_related_dirs_listing"] = out.strip().split("\n")
    json.dump(found, open(OUT + "/kinship_manifest_search.json", "w"), indent=1)
    return found

def step_sibling():
    d = {}
    for name, pat in [("dosages", _config_path("${DOSAGE_DIR}/chr22.vcf.gz")),
                      ("dosages_v2", _config_path("${DOSAGE_V2_DIR}/chr22.vcf.gz"))]:
        e = stat_file(pat)
        if e["exists"] and e.get("readable"):
            h = header_info(pat)
            e["FORMAT_ids"] = h.get("FORMAT_ids"); e["INFO_ids"] = h.get("INFO_ids"); e["phasing"] = h.get("phasing"); e["source"] = h.get("source")
            e["n_samples"] = h.get("n_samples"); e["sample_list_md5_newline_joined"] = h.get("sample_list_md5_newline_joined")
            e["FORMAT_lines"] = h.get("FORMAT")
        d[name] = e
    rc, out, err = sh(_config_path("ls ${DOSAGE_DIR} | wc -l; ls ${DOSAGE_V2_DIR} | wc -l"))
    d["dir_file_counts"] = out.strip().split("\n")
    json.dump(d, open(OUT + "/sibling_dosage_headers.json", "w"), indent=1)
    return d

def step_annot():
    d = {}
    for p in ["%s/deductive/gencode.nochr.gtf.gz" % REF, "%s/deductive/gencode.nochr.gtf.gz.tbi" % REF, "%s/b6_cards/ccre.s.bed" % REF]:
        d[p] = stat_file(p)
    d["re2g_files"] = sorted(os.path.basename(x) for x in glob.glob("%s/re2g/*.bed.gz" % REF))
    d["remap_by_chr_files"] = len(glob.glob(_config_path("${PROJECT_ROOT}/work/fset/out/remap_by_chr/chr*.bed")))
    d["gate1_ctx_json"] = len(glob.glob(_config_path("${PROJECT_ROOT}/work/gate1/out/ctx_chr*.json")))
    d["gate2_g2_json"] = len(glob.glob(_config_path("${PROJECT_ROOT}/work/gate2/out/g2_chr*.json")))
    for m in [_config_path("${PROJECT_ROOT}/work/gate1/out/gate1_chain.done"), _config_path("${PROJECT_ROOT}/work/gate2/out/gate2.done"), _config_path("${PROJECT_ROOT}/work/fset/out/annflag.done")]:
        e = stat_file(m)
        if e["exists"]:
            e["content"] = open(m).read().strip()[:100]
        d[m] = e
    rc, out, err = sh(_config_path("ls ${PROJECT_ROOT}/work/gate1/out | head -40; echo ---; ls ${PROJECT_ROOT}/work/gate1 | head -40"))
    d["gate1_listing"] = out.strip().split("\n")
    json.dump(d, open(OUT + "/annotation_gate_inventory.json", "w"), indent=1)
    return d

if __name__ == "__main__":
    steps = sys.argv[2].split(",") if len(sys.argv) > 2 else ["pheno","kinship","sibling","annot","orig","derived","lift","altmajor"]
    for s in steps:
        log("=== step %s" % s)
        try:
            globals()["step_" + s]()
        except Exception as e:
            log("STEP %s FAILED: %r" % (s, e))
            json.dump({"error": repr(e)}, open(OUT + "/step_%s.error.json" % s, "w"))
    open(OUT + "/inventory.done", "w").write(time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()) + "\n")
    log("DONE")
