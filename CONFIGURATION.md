# Configuration of the research snapshot

Select the scripts and experimental version before configuring their inputs. Settings are
per analysis; the inventory below does not require every variable for every script.
Values from controlled data must be supplied privately by an authorized operator.
Do not commit those values, input files or derived results. `.env` is ignored; scripts
read exported environment variables rather than automatically loading an `.env` file.

## Runtime paths and external tools

`PROJECT_ROOT` is the working-data project directory used by historical analyses.
`DATA_ROOT`, `COHORT_DATA`, `GENOTYPE_DIR`, `DOSAGE_DIR`, `DOSAGE_V2_DIR`,
`PHENO_DIR`, `METADATA_DIR`, `GENOTYPE_PREFIX`, `GENOTYPE_FILE`, `EIGENVEC_FILE`
and the relevant GRM path settings replace source-specific filesystem layouts.
For example, a chromosome input is `${GENOTYPE_DIR}/chr22.vcf.gz`.

Select the environment before execution. Tools normally use `PATH`; specific launchers
also accept `PYTHON_BIN`, `BCFTOOLS_BIN`, `PLINK2_BIN` or `UDOCKER_BIN`.
Legacy environment activation uses `CONDA_INIT_SCRIPT` and `CONDA_ENV`.
Container runs require their runner settings and `SAIGE_IMAGE`; a container-state
folder is configured separately. Set `PRS_CS_ROOT` to an external official PRS-CS
checkout. The 07 launchers resolve source code from their own location and retain
working data under `${PROJECT_ROOT}/work/prs`.

Missing or blank required runtime-path settings fail explicitly. Integer dimensions
and sample counts must be positive integers; use a finite floating point value for
`EXPECTED_BASELINE_R2` and a positive finite value for `MDE_REFERENCE`.
Data-derived expected row counts, case counts, source statistics and coordinate inputs
are settings. `TSS_POSITION` and `GENE_START_B38` are genomic coordinates.
Public statistical constants, seeds, chunk sizes and experimental thresholds are retained.

## Clinical source mappings

Mappings are explicit and scoped to `COHORT_A`, `COHORT_B` and `COHORT_C`.
These labels select distinct source tables and columns; the caller must verify clinical
meaning, coding, units and missing-value rules against their authorized metadata.
No source-column correspondence has been inferred from a generic replacement name.

The phenotype builders in `run_pheno2` and `run_pheno3` use `PHENO_DIR` for relative
source paths and globs. They keep the normalized `CT`/`NC` output covariates used by
downstream scripts. Their disease column lists and case/control encodings are required
settings. The legacy `run_cur2` builder uses per-cohort disease-column maps and a
configured expected-case-count gate. Supply complete mappings for the selected version.

The v8 residual builder reads the configured `COHORT_A/B/C_*_FILE` values as direct
paths, not paths relative to `--base`. The wide trait specification uses table kinds
`SUBJECT`, `LAB` and `ANTHROPOMETRY`. Its historical source schema expects the first
column to contain participant identifiers; numeric sex codes are 1/2. It retains its
missing-value set, declared height unit for BMI, total-cholesterol range gate and OLS
residualization followed by RINT. Column settings are scoped by source cohort. Use the
file-column trait schema when the generic wide schema is inappropriate.

The 07 height preparation instead uses `PHENOTYPE_DATA_ROOT`, `HEIGHT_SOURCES_JSON`
(a list of relative path/column pairs), `COHORT_INDICATOR_COLUMNS_JSON` (two indicator
columns) and `PHENO_COVARIATE_COLUMNS` (comma-separated columns in intended order).
`HEIGHT_SCAN_GLOB` is used by its source-inspection utility. These settings are
version-specific and should not be assumed interchangeable with the phenotype builders.

## Phenotype settings

| Setting | Type | Use |
|---|---|---|
| `N_SAMPLES` | positive integer | Both phenotype builders: expected normalized genotype and PC sample count |
| `GENOTYPE_FILE` | path | Both phenotype builders: genotype VCF, PC input, output directory respectively |
| `EIGENVEC_FILE` | path | Both phenotype builders: genotype VCF, PC input, output directory respectively |
| `OUTPUT_DIR` | path | Both phenotype builders: genotype VCF, PC input, output directory respectively |
| `PHENO_DIR` | directory | Both phenotype builders: base directory used for relative configured source paths and globs |
| `CASE_VALUES_JSON` | JSON string array | Both phenotype builders: nonempty JSON arrays of source code strings; case and control sets must be disjoint |
| `CONTROL_VALUES_JSON` | JSON string array | Both phenotype builders: nonempty JSON arrays of source code strings; case and control sets must be disjoint |
| `MALE_VALUES_JSON` | JSON string array | Both phenotype builders: nonempty JSON arrays of source code strings; case and control sets must be disjoint |
| `MISSING_CODES_JSON` | JSON string array (may be empty) | Both phenotype builders: source missing-value code strings |
| `COHORT_A_SUBJECT_FILE` | path, absolute or relative to PHENO_DIR | Both builders: participant metadata table for cohort A |
| `COHORT_A_SUBJECT_ID_COLUMN` | column name | Both builders: declared participant ID, sex, age column names for cohort A |
| `COHORT_A_SEX_COLUMN` | column name | Both builders: declared participant ID, sex, age column names for cohort A |
| `COHORT_A_AGE_COLUMN` | column name | Both builders: declared participant ID, sex, age column names for cohort A |
| `COHORT_A_DISEASE_GLOB` | glob, absolute or relative to PHENO_DIR | Both builders: glob of source diagnosis tables for cohort A |
| `COHORT_A_DISEASE_ID_COLUMN` | column name | Both builders: declared diagnosis-table participant ID column for cohort A |
| `COHORT_A_HTN_COLUMNS_JSON` | nonempty JSON array of column names | Both builders: explicit diagnosis columns across configured visit tables for cohort A |
| `COHORT_A_DM_COLUMNS_JSON` | nonempty JSON array of column names | Both builders: explicit diagnosis columns across configured visit tables for cohort A |
| `COHORT_A_LIP_COLUMNS_JSON` | nonempty JSON array of column names | Both builders: explicit diagnosis columns across configured visit tables for cohort A |
| `COHORT_A_LAB_FILE` | path, absolute or relative to PHENO_DIR | v3 builder: laboratory table for cohort A |
| `COHORT_A_LAB_ID_COLUMN` | column name | v3 builder: declared laboratory participant ID and total cholesterol columns for cohort A |
| `COHORT_A_TCHL_COLUMN` | column name | v3 builder: declared laboratory participant ID and total cholesterol columns for cohort A |
| `COHORT_B_SUBJECT_FILE` | path, absolute or relative to PHENO_DIR | Both builders: participant metadata table for cohort B |
| `COHORT_B_SUBJECT_ID_COLUMN` | column name | Both builders: declared participant ID, sex, age column names for cohort B |
| `COHORT_B_SEX_COLUMN` | column name | Both builders: declared participant ID, sex, age column names for cohort B |
| `COHORT_B_AGE_COLUMN` | column name | Both builders: declared participant ID, sex, age column names for cohort B |
| `COHORT_B_DISEASE_GLOB` | glob, absolute or relative to PHENO_DIR | Both builders: glob of source diagnosis tables for cohort B |
| `COHORT_B_DISEASE_ID_COLUMN` | column name | Both builders: declared diagnosis-table participant ID column for cohort B |
| `COHORT_B_HTN_COLUMNS_JSON` | nonempty JSON array of column names | Both builders: explicit diagnosis columns across configured visit tables for cohort B |
| `COHORT_B_DM_COLUMNS_JSON` | nonempty JSON array of column names | Both builders: explicit diagnosis columns across configured visit tables for cohort B |
| `COHORT_B_LIP_COLUMNS_JSON` | nonempty JSON array of column names | Both builders: explicit diagnosis columns across configured visit tables for cohort B |
| `COHORT_B_LAB_FILE` | path, absolute or relative to PHENO_DIR | v3 builder: laboratory table for cohort B |
| `COHORT_B_LAB_ID_COLUMN` | column name | v3 builder: declared laboratory participant ID and total cholesterol columns for cohort B |
| `COHORT_B_TCHL_COLUMN` | column name | v3 builder: declared laboratory participant ID and total cholesterol columns for cohort B |
| `COHORT_C_SUBJECT_FILE` | path, absolute or relative to PHENO_DIR | Both builders: participant metadata table for cohort C |
| `COHORT_C_SUBJECT_ID_COLUMN` | column name | Both builders: declared participant ID, sex, age column names for cohort C |
| `COHORT_C_SEX_COLUMN` | column name | Both builders: declared participant ID, sex, age column names for cohort C |
| `COHORT_C_AGE_COLUMN` | column name | Both builders: declared participant ID, sex, age column names for cohort C |
| `COHORT_C_DISEASE_GLOB` | glob, absolute or relative to PHENO_DIR | Both builders: glob of source diagnosis tables for cohort C |
| `COHORT_C_DISEASE_ID_COLUMN` | column name | Both builders: declared diagnosis-table participant ID column for cohort C |
| `COHORT_C_HTN_COLUMNS_JSON` | nonempty JSON array of column names | Both builders: explicit diagnosis columns across configured visit tables for cohort C |
| `COHORT_C_DM_COLUMNS_JSON` | nonempty JSON array of column names | Both builders: explicit diagnosis columns across configured visit tables for cohort C |
| `COHORT_C_LIP_COLUMNS_JSON` | nonempty JSON array of column names | Both builders: explicit diagnosis columns across configured visit tables for cohort C |
| `COHORT_C_LAB_FILE` | path, absolute or relative to PHENO_DIR | v3 builder: laboratory table for cohort C |
| `COHORT_C_LAB_ID_COLUMN` | column name | v3 builder: declared laboratory participant ID and total cholesterol columns for cohort C |
| `COHORT_C_TCHL_COLUMN` | column name | v3 builder: declared laboratory participant ID and total cholesterol columns for cohort C |
| `PHENO_V1_HTN_CASES` | nonnegative integer | v2 builder: previous-version case count used in the original merge gate |
| `PHENO_V1_DM_CASES` | nonnegative integer | v2 builder: previous-version case count used in the original merge gate |
| `PHENO_V1_LIP_CASES` | nonnegative integer | v2 builder: previous-version case count used in the original merge gate |
| `PHENO_ADDED_COHORT_MIN_HTN_CASES` | nonnegative integer | v2 builder: required additional case count used in the original merge gate |
| `PHENO_ADDED_COHORT_MIN_DM_CASES` | nonnegative integer | v2 builder: required additional case count used in the original merge gate |
| `PHENO_ADDED_COHORT_MIN_LIP_CASES` | nonnegative integer | v2 builder: required additional case count used in the original merge gate |
| `PHENO_V2_HTN_CASES` | nonnegative integer | v3 builder: exact prior-version binary case count gate |
| `PHENO_V2_DM_CASES` | nonnegative integer | v3 builder: exact prior-version binary case count gate |
| `PHENO_V2_LIP_CASES` | nonnegative integer | v3 builder: exact prior-version binary case count gate |
| `TCHL_MIN_VALUE` | finite float | v3 builder: exclusive lower/upper source-unit cholesterol bounds; lower must be less than upper |
| `TCHL_MAX_VALUE` | finite float | v3 builder: exclusive lower/upper source-unit cholesterol bounds; lower must be less than upper |
| `MIN_TCHL_SAMPLES` | positive integer | v3 builder: minimum retained quantitative-phenotype sample count gate |
| `BCFTOOLS` | executable | Builders and bcftools-using pilot helpers: executable path or PATH command |
| `PROJECT_ROOT` | directory | Pilot Python helpers except agetest.py, and all shell launchers: generic workspace root |
| `PHENO_V2_DIR` | directory | agetest.py and run_step1v2/step1.sh: configured v2 output tables |
| `GRM_FILE` | path | step1 launchers and run_pilot/pr.sh: sparse GRM and its sample-ID file |
| `GRM_SAMPLE_IDS_FILE` | path | step1 launchers and run_pilot/pr.sh: sparse GRM and its sample-ID file |
| `PLINK_PREFIX` | path prefix | step1 launchers: PLINK BED/BIM/FAM input prefix |
| `SAIGE_STEP1_SCRIPT` | path | step1 launchers: installed SAIGE step1 script |
| `SAIGE_STEP2_SCRIPT` | path | Pilot SAIGE shell launchers: installed SAIGE step2 script |
| `SAIGE_RUNNER` | executable | SAIGE shell launchers: a single executable accepting the script path and original flags; no command-string eval |
| `REFERENCE_DIR` | directory | Pilot SAIGE launchers and step1 v3/v4: reference-data directory |
| `RUN_DIR` | directory | Pilot SAIGE launchers and step1 v3/v4: lock/log workspace |
| `PHENO_V3_DIR` | directory | step1 v3/v4: configured v3 output tables |
| `PYTHON` | executable | run_bwg.sh: Python executable for the sibling bwg.py |
| `NTHREADS` | integer | step1 v3/v4: SAIGE thread count |
| `POOL` | integer | pr.sh / pr2.sh: parallel job count |
| `NCHUNK` | integer | bchunk.sh: group-file split and parallel job count |
| `ARMS` | string | pr2.sh: whitespace-separated group-file arms |
| `NULL_MODEL_FILE` | path | readmod.R: model file to inspect |

## Literal settings inventory

This index records literal environment-variable references. Dynamic per-cohort settings
are described above and in the selected script. Some entries are local shell variables
or optional tool/resource overrides; read the script for its required/default behavior.

| Name | Referenced in |
|---|---|
| `ALLOW_HEAVY` | `utils/reusable/build_r2_only_broad_liftover_input.sh` |
| `ANALYSIS_ROOT` | `experiments/03_audits_and_null_evaluation/audit.py`<br>`experiments/03_audits_and_null_evaluation/idmap.py`<br>`experiments/03_audits_and_null_evaluation/rd.R`<br>and 8 other scripts |
| `ARM` | `pipeline/annotation/run_annot/l1_chain.sh`<br>`pipeline/phenotype_nullmodel/run_pilot/pr.sh` |
| `ARMS` | `pipeline/phenotype_nullmodel`<br>`pipeline/phenotype_nullmodel/run_pilot/pr2.sh` |
| `AV` | `experiments/02_association_and_finemapping/scripts/fm_chr.sh`<br>`experiments/02_association_and_finemapping/scripts/fm_pool.sh`<br>`experiments/05_learned_regulatory_sets/att/wd_att.sh`<br>and 6 other scripts |
| `BAGG_OUT` | `experiments/01_annotation_weighting_models/bagg2.py`<br>`experiments/01_annotation_weighting_models/bagg2b.py`<br>`experiments/01_annotation_weighting_models/bagg2b_diag.py` |
| `BAGG_SAVE_LOCK` | `experiments/01_annotation_weighting_models/bagg2b.py`<br>`experiments/01_annotation_weighting_models/bagg2b_diag.py` |
| `BCFTOOLS` | `experiments/03_audits_and_null_evaluation/audit.py`<br>`experiments/03_audits_and_null_evaluation/ca.py`<br>`experiments/03_audits_and_null_evaluation/cond1.py`<br>and 57 other scripts |
| `BCFTOOLS_BIN` | `experiments/07_attention_based_variant_interactions/baselines/prs_cs/extract_hm3.sh`<br>`experiments/07_attention_based_variant_interactions/baselines/prs_cs/extract_hm3_v2.sh`<br>`experiments/07_attention_based_variant_interactions/data_prep/build_height.py`<br>and 17 other scripts |
| `BEDTOOLS` | `pipeline/groupfiles_variants/run_b6/b6_card_join.py`<br>`pipeline/groupfiles_variants/run_g5/g5p.sh` |
| `BGZIP` | `pipeline/groupfiles_variants/run_g5/g5c2.sh` |
| `BUDGET` | `pipeline/annotation/run_annot_bbj/wd_rss.sh` |
| `CASE_VALUES_JSON` | `pipeline/phenotype_nullmodel`<br>`pipeline/phenotype_nullmodel/run_pheno2/build.py`<br>`pipeline/phenotype_nullmodel/run_pheno3/build.py` |
| `CHAIN` | `utils/reusable/build_r2_only_broad_liftover_input.sh` |
| `CHRS` | `experiments/05_learned_regulatory_sets/settest/scripts/s4_assoc.sh` |
| `CODE_ROOT` | `experiments/07_attention_based_variant_interactions/baselines/prs_cs/prscs_chr.sh`<br>`experiments/07_attention_based_variant_interactions/infra/queue_p5.sh`<br>`experiments/07_attention_based_variant_interactions/infra/queue_p5_v2.sh`<br>and 31 other scripts |
| `COHORT_A_AGE_COLUMN` | `pipeline/phenotype_nullmodel` |
| `COHORT_A_BIOCHEM_FILE` | `experiments/03_audits_and_null_evaluation/cb6.py`<br>`experiments/03_audits_and_null_evaluation/spec.py`<br>`experiments/03_audits_and_null_evaluation/tchl.py`<br>and 1 other scripts |
| `COHORT_A_DISEASE_FILE` | `experiments/03_audits_and_null_evaluation/cols.py` |
| `COHORT_A_DISEASE_GLOB` | `pipeline/phenotype_nullmodel` |
| `COHORT_A_DISEASE_ID_COLUMN` | `pipeline/phenotype_nullmodel` |
| `COHORT_A_DM_COLUMNS_JSON` | `pipeline/phenotype_nullmodel` |
| `COHORT_A_HTN_COLUMNS_JSON` | `pipeline/phenotype_nullmodel` |
| `COHORT_A_ID_COLUMN` | `experiments/03_audits_and_null_evaluation/tchl.py` |
| `COHORT_A_LAB_FILE` | `pipeline/phenotype_nullmodel` |
| `COHORT_A_LAB_ID_COLUMN` | `pipeline/phenotype_nullmodel` |
| `COHORT_A_LIP_COLUMNS_JSON` | `pipeline/phenotype_nullmodel` |
| `COHORT_A_METADATA_FILE` | `experiments/03_audits_and_null_evaluation/cb5.py`<br>`experiments/03_audits_and_null_evaluation/cb6.py`<br>`pipeline/groupfiles_variants/run_cur/htncu.py`<br>and 1 other scripts |
| `COHORT_A_SEX_COLUMN` | `pipeline/phenotype_nullmodel` |
| `COHORT_A_SUBJECT_FILE` | `experiments/03_audits_and_null_evaluation/find_missing.py`<br>`experiments/03_audits_and_null_evaluation/ov.py`<br>`experiments/03_audits_and_null_evaluation/spec.py`<br>and 1 other scripts |
| `COHORT_A_SUBJECT_ID_COLUMN` | `pipeline/phenotype_nullmodel` |
| `COHORT_A_TCHL_COLUMN` | `experiments/03_audits_and_null_evaluation/cb6.py`<br>`experiments/03_audits_and_null_evaluation/spec.py`<br>`experiments/03_audits_and_null_evaluation/tchl.py`<br>and 2 other scripts |
| `COHORT_A_TCHL_INSTRUMENT_COLUMN` | `experiments/03_audits_and_null_evaluation/cb6.py` |
| `COHORT_B_AGE_COLUMN` | `pipeline/phenotype_nullmodel` |
| `COHORT_B_BIOCHEM_FILE` | `experiments/03_audits_and_null_evaluation/spec.py`<br>`experiments/03_audits_and_null_evaluation/tchl.py`<br>`experiments/03_audits_and_null_evaluation/tchl2.py` |
| `COHORT_B_DISEASE_FILE` | `experiments/03_audits_and_null_evaluation/cols.py` |
| `COHORT_B_DISEASE_GLOB` | `pipeline/phenotype_nullmodel` |
| `COHORT_B_DISEASE_ID_COLUMN` | `pipeline/phenotype_nullmodel` |
| `COHORT_B_DM_COLUMNS_JSON` | `pipeline/phenotype_nullmodel` |
| `COHORT_B_HEIGHT_COLUMN` | `experiments/01_annotation_weighting_models/model_v8_out/build_resid_v8.py` |
| `COHORT_B_HTN_COLUMNS_JSON` | `pipeline/phenotype_nullmodel` |
| `COHORT_B_ID_COLUMN` | `experiments/03_audits_and_null_evaluation/tchl.py` |
| `COHORT_B_LAB_FILE` | `pipeline/phenotype_nullmodel` |
| `COHORT_B_LAB_ID_COLUMN` | `pipeline/phenotype_nullmodel` |
| `COHORT_B_LIP_COLUMNS_JSON` | `pipeline/phenotype_nullmodel` |
| `COHORT_B_METADATA_FILE` | `experiments/03_audits_and_null_evaluation/cb5.py`<br>`pipeline/groupfiles_variants/run_cur/htncu.py`<br>`pipeline/groupfiles_variants/run_cur/metadata_scan.py` |
| `COHORT_B_SEX_COLUMN` | `pipeline/phenotype_nullmodel` |
| `COHORT_B_SUBJECT_FILE` | `experiments/03_audits_and_null_evaluation/find_missing.py`<br>`experiments/03_audits_and_null_evaluation/ov.py`<br>`experiments/03_audits_and_null_evaluation/spec.py`<br>and 1 other scripts |
| `COHORT_B_SUBJECT_ID_COLUMN` | `pipeline/phenotype_nullmodel` |
| `COHORT_B_TCHL_COLUMN` | `experiments/03_audits_and_null_evaluation/spec.py`<br>`experiments/03_audits_and_null_evaluation/tchl.py`<br>`experiments/03_audits_and_null_evaluation/tchl2.py`<br>and 1 other scripts |
| `COHORT_B_WEIGHT_COLUMN` | `experiments/01_annotation_weighting_models/model_v8_out/build_resid_v8.py` |
| `COHORT_C_AGE_COLUMN` | `pipeline/phenotype_nullmodel` |
| `COHORT_C_BIOCHEM_FILE` | `experiments/03_audits_and_null_evaluation/cb4.py`<br>`experiments/03_audits_and_null_evaluation/spec.py`<br>`experiments/03_audits_and_null_evaluation/tchl2.py` |
| `COHORT_C_DISEASE_COLUMN_REGEX` | `pipeline/groupfiles_variants/run_cur/as_cb.py`<br>`pipeline/groupfiles_variants/run_cur/as_cb2.py`<br>`pipeline/groupfiles_variants/run_cur/as_cb4.py` |
| `COHORT_C_DISEASE_FILE` | `experiments/03_audits_and_null_evaluation/cols.py` |
| `COHORT_C_DISEASE_GLOB` | `pipeline/phenotype_nullmodel` |
| `COHORT_C_DISEASE_ID_COLUMN` | `pipeline/phenotype_nullmodel` |
| `COHORT_C_DM_COLUMNS_JSON` | `pipeline/phenotype_nullmodel` |
| `COHORT_C_HTN_COLUMNS_JSON` | `pipeline/phenotype_nullmodel` |
| `COHORT_C_LAB_FILE` | `pipeline/phenotype_nullmodel` |
| `COHORT_C_LAB_ID_COLUMN` | `pipeline/phenotype_nullmodel` |
| `COHORT_C_LIP_COLUMNS_JSON` | `pipeline/phenotype_nullmodel` |
| `COHORT_C_METADATA_FILE` | `experiments/03_audits_and_null_evaluation/cb.py`<br>`experiments/03_audits_and_null_evaluation/cb2.py`<br>`experiments/03_audits_and_null_evaluation/cb3.py`<br>and 7 other scripts |
| `COHORT_C_SEX_COLUMN` | `pipeline/phenotype_nullmodel` |
| `COHORT_C_SUBJECT_FILE` | `experiments/03_audits_and_null_evaluation/find_missing.py`<br>`experiments/03_audits_and_null_evaluation/ov.py`<br>`experiments/03_audits_and_null_evaluation/spec.py`<br>and 1 other scripts |
| `COHORT_C_SUBJECT_ID_COLUMN` | `pipeline/phenotype_nullmodel` |
| `COHORT_C_TCHL_COLUMN` | `experiments/03_audits_and_null_evaluation/cb4.py`<br>`experiments/03_audits_and_null_evaluation/spec.py`<br>`experiments/03_audits_and_null_evaluation/tchl2.py`<br>and 1 other scripts |
| `COHORT_C_TCHL_ORIGINAL_COLUMN` | `experiments/03_audits_and_null_evaluation/cb4.py`<br>`experiments/03_audits_and_null_evaluation/tchl2.py` |
| `COHORT_DATA` | `experiments/03_audits_and_null_evaluation/gate0_id_overlap.py`<br>`experiments/03_audits_and_null_evaluation/gate0b_meth_vs_exomechip.py`<br>`experiments/03_audits_and_null_evaluation/gate0c_hexa_vcf_profile.py`<br>and 3 other scripts |
| `COHORT_D_SUBJECT_FILE` | `experiments/03_audits_and_null_evaluation/ov.py` |
| `COHORT_D_TABLE_GLOB` | `experiments/03_audits_and_null_evaluation/ov.py` |
| `COHORT_INDICATOR_COLUMNS_JSON` | `experiments/07_attention_based_variant_interactions/data_prep/build_height.py` |
| `CONDA_ENV` | `experiments/01_annotation_weighting_models/run_alpha_chain2.sh`<br>`experiments/01_annotation_weighting_models/run_alpha_chain3.sh`<br>`experiments/01_annotation_weighting_models/run_gf_build.sh`<br>and 6 other scripts |
| `CONDA_INIT_SCRIPT` | `experiments/01_annotation_weighting_models/run_alpha_chain2.sh`<br>`experiments/01_annotation_weighting_models/run_alpha_chain3.sh`<br>`experiments/01_annotation_weighting_models/run_gf_build.sh`<br>and 6 other scripts |
| `CONDA_PREFIX` | `experiments/02_association_and_finemapping/scripts/fm_region.sh`<br>`experiments/02_association_and_finemapping/scripts/probe.sh`<br>`experiments/02_association_and_finemapping/scripts/salvage_tmp.sh` |
| `CONTAINER_STATE_DIR` | `experiments/01_annotation_weighting_models/run_arm_seq.sh`<br>`experiments/01_annotation_weighting_models/run_chunked.sh`<br>`experiments/01_annotation_weighting_models/run_saige_arm.sh`<br>and 21 other scripts |
| `CONTROL_VALUES_JSON` | `pipeline/phenotype_nullmodel`<br>`pipeline/phenotype_nullmodel/run_pheno2/build.py`<br>`pipeline/phenotype_nullmodel/run_pheno3/build.py` |
| `CROSSMAP` | `utils/reusable/build_r2_only_broad_liftover_input.sh` |
| `CTXGLOB` | `experiments/04_context_feasibility_gates/build_gate1_server.py` |
| `CUDA_CACHE_PATH` | `experiments/01_annotation_weighting_models/model_v10_out/l1_train_v10.py`<br>`experiments/01_annotation_weighting_models/model_v8_out/l1_train_v8.py`<br>`experiments/01_annotation_weighting_models/model_v9_out/l1_train_v9.py` |
| `DATASET_ITEM_01` | `pipeline/groupfiles_variants/run_proc2/proc2.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3b.sh`<br>and 1 other scripts |
| `DATASET_ITEM_02` | `pipeline/groupfiles_variants/run_proc2/proc2.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3b.sh`<br>and 1 other scripts |
| `DATASET_ITEM_03` | `pipeline/groupfiles_variants/run_proc2/proc2.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3b.sh`<br>and 1 other scripts |
| `DATASET_ITEM_04` | `pipeline/groupfiles_variants/run_proc2/proc2.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3b.sh`<br>and 1 other scripts |
| `DATASET_ITEM_05` | `pipeline/groupfiles_variants/run_proc2/proc2.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3b.sh`<br>and 1 other scripts |
| `DATASET_ITEM_06` | `pipeline/groupfiles_variants/run_proc2/proc2.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3b.sh`<br>and 1 other scripts |
| `DATASET_ITEM_07` | `pipeline/groupfiles_variants/run_proc2/proc2.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3b.sh`<br>and 1 other scripts |
| `DATASET_ITEM_08` | `pipeline/groupfiles_variants/run_proc2/proc2.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3b.sh`<br>and 1 other scripts |
| `DATASET_ITEM_09` | `pipeline/groupfiles_variants/run_proc2/proc2.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3b.sh`<br>and 1 other scripts |
| `DATASET_ITEM_10` | `pipeline/groupfiles_variants/run_proc2/proc2.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3b.sh`<br>and 1 other scripts |
| `DATASET_ITEM_11` | `pipeline/groupfiles_variants/run_proc2/proc2.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3b.sh`<br>and 1 other scripts |
| `DATASET_ITEM_12` | `pipeline/groupfiles_variants/run_proc2/proc2.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3b.sh`<br>and 1 other scripts |
| `DATASET_ITEM_13` | `pipeline/groupfiles_variants/run_proc2/proc2.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3b.sh`<br>and 1 other scripts |
| `DATASET_ITEM_14` | `pipeline/groupfiles_variants/run_proc2/proc2.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3b.sh`<br>and 1 other scripts |
| `DATASET_ITEM_15` | `pipeline/groupfiles_variants/run_proc2/proc2.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3b.sh`<br>and 1 other scripts |
| `DATASET_ITEM_16` | `pipeline/groupfiles_variants/run_proc2/proc2.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3b.sh`<br>and 1 other scripts |
| `DATASET_ITEM_17` | `pipeline/groupfiles_variants/run_proc2/proc2.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3b.sh`<br>and 1 other scripts |
| `DATASET_ITEM_18` | `pipeline/groupfiles_variants/run_proc2/proc2.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3b.sh`<br>and 1 other scripts |
| `DATASET_ITEM_19` | `pipeline/groupfiles_variants/run_proc2/proc2.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3b.sh`<br>and 1 other scripts |
| `DATASET_ITEM_20` | `pipeline/groupfiles_variants/run_proc2/proc2.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3b.sh`<br>and 1 other scripts |
| `DATASET_ITEM_21` | `pipeline/groupfiles_variants/run_proc2/proc2.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3b.sh`<br>and 1 other scripts |
| `DATASET_ITEM_67` | `pipeline/groupfiles_variants/run_proc2/proc2.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3.sh`<br>`pipeline/groupfiles_variants/run_proc3/proc3b.sh`<br>and 1 other scripts |
| `DATA_ROOT` | `experiments/03_audits_and_null_evaluation/verify.py`<br>`experiments/03_audits_and_null_evaluation/vrpool.py`<br>`experiments/06_haplotype_context/run_20260922T052403Z/code/hcf_inventory.py`<br>and 1 other scripts |
| `DF` | `experiments/05_learned_regulatory_sets/att/wd_b200.sh` |
| `DISEASE_COLUMN_REGEX` | `experiments/03_audits_and_null_evaluation/cols.py` |
| `DOMSTEP` | `experiments/05_learned_regulatory_sets/att/att_common_cand.py`<br>`experiments/05_learned_regulatory_sets/att/att_common_oracle.py`<br>`experiments/05_learned_regulatory_sets/att/att_common_oracle_hq.py`<br>and 5 other scripts |
| `DOSAGE_DIR` | `experiments/06_haplotype_context/run_20260922T052403Z/code/hcf_inventory.py` |
| `DOSAGE_V2_DIR` | `experiments/06_haplotype_context/run_20260922T052403Z/code/hcf_inventory.py` |
| `DUPLICATE_CANONICAL_KEY_N` | `utils/reusable/build_r2_only_broad_liftover_input.sh` |
| `EIGENVEC_FILE` | `experiments/01_annotation_weighting_models/model_v8_out/build_resid_v8.py`<br>`experiments/03_audits_and_null_evaluation/spec.py`<br>`pipeline/groupfiles_variants/run_cur2/build.py`<br>and 3 other scripts |
| `EXOME_PED_FILE` | `experiments/03_audits_and_null_evaluation/gate0b_meth_vs_exomechip.py` |
| `EXPECTED_BASELINE_R2` | `experiments/07_attention_based_variant_interactions/models/attention/train_att_hei_res.py`<br>`experiments/07_attention_based_variant_interactions/models/attention/train_att_hei_res_v2.py`<br>`experiments/07_attention_based_variant_interactions/models/attention/train_att_hei_res_v3.py`<br>and 4 other scripts |
| `EXPECTED_CASE_COUNTS_JSON` | `pipeline/groupfiles_variants/run_cur2/build.py` |
| `EXPECTED_TOTAL_RECORDS` | `utils/reusable/aggregate_g2_fmi_20260726.py` |
| `EXP_OUT` | `utils/reusable/exp/exp_f2_validator_v3_20260727.py`<br>`utils/reusable/exp/exp_f3a_cmac_distribution.py`<br>`utils/reusable/exp/exp_f3b_gene_universe.py`<br>and 1 other scripts |
| `FINAL_STATUS` | `utils/reusable/build_r2_only_broad_liftover_input.sh` |
| `FM_DEADLINE` | `experiments/02_association_and_finemapping/scripts/fm_pool.sh` |
| `FM_LAMBDA` | `experiments/02_association_and_finemapping/scripts/fm_susie.R` |
| `FM_MAXVAR_MEM` | `experiments/01_annotation_weighting_models/scripts/fm_reweight.R`<br>`experiments/01_annotation_weighting_models/scripts/fm_susie_prior.R`<br>`experiments/01_annotation_weighting_models/scripts/fm_susie_prior2.R`<br>and 5 other scripts |
| `FM_MAX_ITER` | `experiments/01_annotation_weighting_models/scripts/fm_reweight.R`<br>`experiments/01_annotation_weighting_models/scripts/fm_susie_prior.R`<br>`experiments/01_annotation_weighting_models/scripts/fm_susie_prior2.R`<br>and 1 other scripts |
| `FM_UNIFORM_PIP` | `experiments/01_annotation_weighting_models/scripts/fm_reweight.R`<br>`experiments/01_annotation_weighting_models/scripts/fm_susie_prior.R`<br>`experiments/01_annotation_weighting_models/scripts/fm_susie_prior2.R` |
| `F_BLOCKNULL` | `experiments/05_learned_regulatory_sets/f1f2/f1_tendency.py` |
| `F_JOBS` | `experiments/05_learned_regulatory_sets/f1f2/f2_modelC.py`<br>`experiments/05_learned_regulatory_sets/f1f2/f2_modelU.py` |
| `F_NQ` | `experiments/05_learned_regulatory_sets/f1f2/f1_tendency.py` |
| `F_OUT` | `experiments/05_learned_regulatory_sets/f1f2/f_common.py` |
| `F_PERDOM` | `experiments/05_learned_regulatory_sets/f1f2/f1_tendency.py` |
| `F_PRIM` | `experiments/05_learned_regulatory_sets/f1f2/f_common.py` |
| `F_WORKERS` | `experiments/05_learned_regulatory_sets/f1f2/f1_tendency.py`<br>`experiments/05_learned_regulatory_sets/f1f2/f2_modelC.py`<br>`experiments/05_learned_regulatory_sets/f1f2/f2_modelU.py` |
| `GCV_LOGMAX` | `pipeline/annotation/run_annot/l1_train.py` |
| `GCV_TRACE` | `pipeline/annotation/run_annot/l1_train.py` |
| `GENES` | `experiments/02_association_and_finemapping/LDLR2nd_v4/v4.sh` |
| `GENE_START_B38` | `experiments/03_audits_and_null_evaluation/coordchk.py` |
| `GENOTYPE_DIR` | `experiments/02_association_and_finemapping/LDLR2nd_v4/v4.sh`<br>`experiments/02_association_and_finemapping/cmp_hexa.sh`<br>`experiments/02_association_and_finemapping/cmp_hexa2.sh`<br>and 20 other scripts |
| `GENOTYPE_FILE` | `pipeline/phenotype_nullmodel`<br>`pipeline/phenotype_nullmodel/run_pheno2/build.py`<br>`pipeline/phenotype_nullmodel/run_pheno3/build.py` |
| `GENOTYPE_FILENAME_PREFIX` | `pipeline/groupfiles_variants/run_vr/calib.py`<br>`pipeline/groupfiles_variants/run_vr/calib2.py`<br>`pipeline/groupfiles_variants/run_vr/hc.sh`<br>and 8 other scripts |
| `GENOTYPE_FILENAME_SUFFIX` | `pipeline/groupfiles_variants/run_vr/calib.py`<br>`pipeline/groupfiles_variants/run_vr/calib2.py`<br>`pipeline/groupfiles_variants/run_vr/hc.sh`<br>and 8 other scripts |
| `GENOTYPE_PREFIX` | `experiments/02_association_and_finemapping/cmp_hexa.sh`<br>`experiments/02_association_and_finemapping/cmp_hexa2.sh` |
| `GENOTYPE_VCF_FILE` | `pipeline/groupfiles_variants/run_cur/census.py`<br>`pipeline/groupfiles_variants/run_cur2/build.py` |
| `GRM_FILE` | `experiments/03_audits_and_null_evaluation/verify2.py`<br>`pipeline/groupfiles_variants/run_cur3/step1cur.sh`<br>`pipeline/phenotype_nullmodel`<br>and 4 other scripts |
| `GRM_SAMPLE_IDS_FILE` | `pipeline/phenotype_nullmodel`<br>`pipeline/phenotype_nullmodel/run_pilot/pr.sh`<br>`pipeline/phenotype_nullmodel/run_step1v2/step1.sh`<br>and 2 other scripts |
| `GTAG` | `experiments/05_learned_regulatory_sets/settest/scripts/s4_assoc.sh` |
| `HEAD` | `pipeline/annotation/run_annot/l1_train.py` |
| `HEIGHT_SCAN_GLOB` | `experiments/07_attention_based_variant_interactions/data_prep/peek_height.py` |
| `HEIGHT_SOURCES_JSON` | `experiments/07_attention_based_variant_interactions/data_prep/build_height.py`<br>`experiments/07_attention_based_variant_interactions/data_prep/peek_height.py` |
| `HEI_LOGTEST` | `experiments/07_attention_based_variant_interactions/models/attention/train_att_hei_res_v2.py`<br>`experiments/07_attention_based_variant_interactions/models/attention/train_att_hei_res_v3.py`<br>`experiments/07_attention_based_variant_interactions/models/attention/train_att_hei_res_v4.py`<br>and 2 other scripts |
| `HEI_SMOKE` | `experiments/07_attention_based_variant_interactions/data_prep/hei_seed_extract_gt_merge.py`<br>`experiments/07_attention_based_variant_interactions/data_prep/hei_seed_extract_gt_part.py`<br>`experiments/07_attention_based_variant_interactions/data_prep/hei_seed_extract_merge.py`<br>and 1 other scripts |
| `HEXA_IDAT_DIR` | `experiments/03_audits_and_null_evaluation/gate0_id_overlap.py`<br>`experiments/03_audits_and_null_evaluation/gate0b_meth_vs_exomechip.py`<br>`experiments/03_audits_and_null_evaluation/gate0c_hexa_vcf_profile.py`<br>and 1 other scripts |
| `HEXA_ID_FILE` | `experiments/03_audits_and_null_evaluation/gate0_id_overlap.py` |
| `HEXA_VCF_FILE` | `experiments/03_audits_and_null_evaluation/gate0c_hexa_vcf_profile.py`<br>`experiments/03_audits_and_null_evaluation/run_gate1_chr22.sh` |
| `ID_COLUMN` | `experiments/03_audits_and_null_evaluation/audit.py`<br>`experiments/03_audits_and_null_evaluation/find_missing.py`<br>`experiments/03_audits_and_null_evaluation/idmap.py`<br>and 3 other scripts |
| `INT` | `pipeline/annotation/run_annot_bbj/wd_rss.sh` |
| `INTERVAL` | `pipeline/annotation/run_annot/wd5.sh` |
| `KARE_HIGH_DENSITY_IDAT_DIR` | `experiments/03_audits_and_null_evaluation/gate0_id_overlap.py`<br>`experiments/03_audits_and_null_evaluation/gate0b_meth_vs_exomechip.py`<br>`experiments/03_audits_and_null_evaluation/gate0c_hexa_vcf_profile.py` |
| `KARE_LOW_DENSITY_IDAT_DIR` | `experiments/03_audits_and_null_evaluation/gate0_id_overlap.py`<br>`experiments/03_audits_and_null_evaluation/gate0b_meth_vs_exomechip.py`<br>`experiments/03_audits_and_null_evaluation/gate0c_hexa_vcf_profile.py` |
| `KILLS0` | `experiments/08_height_learning_value_gate/sup_phi_v9b.sh` |
| `LABEL` | `utils/reusable/build_r2_only_broad_liftover_input.sh` |
| `LD_THREADS` | `experiments/02_association_and_finemapping/scripts/fm_ld.py` |
| `LIFTED_PATH` | `utils/reusable/build_r2_only_broad_liftover_input.sh` |
| `LIFTOVER_MAPPED_N` | `utils/reusable/build_r2_only_broad_liftover_input.sh` |
| `LIFTOVER_REJECT_STATUS` | `utils/reusable/build_r2_only_broad_liftover_input.sh` |
| `LIFTOVER_STATUS` | `utils/reusable/build_r2_only_broad_liftover_input.sh` |
| `LIFTOVER_THREADS` | `utils/reusable/build_r2_only_broad_liftover_input.sh` |
| `LIMIT` | `experiments/05_learned_regulatory_sets/f3f6/wd_f3f6.sh` |
| `MAF_HI` | `pipeline/groupfiles_variants/run_vr/smoke.sh`<br>`pipeline/groupfiles_variants/run_vr/vr.sh` |
| `MAF_LO` | `pipeline/groupfiles_variants/run_vr/smoke.sh`<br>`pipeline/groupfiles_variants/run_vr/vr.sh` |
| `MALE_VALUES_JSON` | `pipeline/phenotype_nullmodel`<br>`pipeline/phenotype_nullmodel/run_pheno2/build.py`<br>`pipeline/phenotype_nullmodel/run_pheno3/build.py` |
| `MASK_N_SAMPLES` | `utils/reusable/build_chr22_common_gene_mask_v2.py` |
| `MAXMAF` | `experiments/05_learned_regulatory_sets/settest/scripts/s4_assoc.sh` |
| `MAXV` | `experiments/05_learned_regulatory_sets/att/interact_ceiling.py`<br>`experiments/05_learned_regulatory_sets/att/run_interact.sh` |
| `MAXVAR` | `experiments/02_association_and_finemapping/scripts/fm_region.sh` |
| `MAX_JOB_MEM_GB` | `utils/reusable/run_with_resource_guard.sh` |
| `MAX_LOAD_PER_CORE` | `utils/reusable/run_with_resource_guard.sh` |
| `MAX_OTHER_USER_CPU_PCT` | `utils/reusable/run_with_resource_guard.sh` |
| `MAX_OTHER_USER_RSS_GB` | `utils/reusable/run_with_resource_guard.sh` |
| `MAX_TOTAL_OTHER_CPU_PCT` | `utils/reusable/run_with_resource_guard.sh` |
| `MAX_TOTAL_OTHER_RSS_GB` | `utils/reusable/run_with_resource_guard.sh` |
| `MDE_REFERENCE` | `experiments/01_annotation_weighting_models/l3_burden.py`<br>`experiments/01_annotation_weighting_models/l3_metrics.py`<br>`experiments/05_learned_regulatory_sets/att/att_judge.py` |
| `MEMMIN` | `experiments/05_learned_regulatory_sets/settest/scripts/wd.sh` |
| `MEM_AVAIL_MIN_GB` | `pipeline/annotation/run_annot/wd5.sh` |
| `METADATA_COLUMN_REGEX` | `pipeline/groupfiles_variants/run_cur/metadata_scan.py` |
| `METADATA_DIR` | `experiments/03_audits_and_null_evaluation/cb.py`<br>`experiments/03_audits_and_null_evaluation/cb2.py`<br>`experiments/03_audits_and_null_evaluation/cb3.py`<br>and 11 other scripts |
| `METADATA_ORIGINAL_LABEL_INDEX` | `experiments/03_audits_and_null_evaluation/cb4.py` |
| `METADATA_SEARCH_TERM` | `pipeline/groupfiles_variants/run_cur/as_cb3.py`<br>`pipeline/groupfiles_variants/run_cur/htncu.py` |
| `METADATA_SHEET` | `pipeline/groupfiles_variants/run_cur/as_cb2.py`<br>`pipeline/groupfiles_variants/run_cur/as_cb3.py`<br>`pipeline/groupfiles_variants/run_cur/as_cb4.py` |
| `METADATA_TRANSFORM_LABEL_INDEX` | `experiments/03_audits_and_null_evaluation/cb4.py` |
| `METHYLATION_FILE` | `experiments/03_audits_and_null_evaluation/gate0b_meth_vs_exomechip.py` |
| `MINV` | `experiments/05_learned_regulatory_sets/att/interact_ceiling.py`<br>`experiments/05_learned_regulatory_sets/att/run_interact.sh` |
| `MIN_AVAIL_MEM_GB` | `utils/reusable/run_with_resource_guard.sh` |
| `MIN_DATA_FREE_GB` | `utils/reusable/run_with_resource_guard.sh` |
| `MIN_TCHL_SAMPLES` | `pipeline/phenotype_nullmodel`<br>`pipeline/phenotype_nullmodel/run_pheno3/build.py` |
| `MIN_TMP_FREE_GB` | `utils/reusable/run_with_resource_guard.sh` |
| `MISSING_CODES_JSON` | `pipeline/phenotype_nullmodel`<br>`pipeline/phenotype_nullmodel/run_pheno2/build.py`<br>`pipeline/phenotype_nullmodel/run_pheno3/build.py` |
| `MKL_NUM_THREADS` | `experiments/05_learned_regulatory_sets/f1f2/f1_tendency.py` |
| `NCHUNK` | `pipeline/phenotype_nullmodel`<br>`pipeline/phenotype_nullmodel/run_pilot/bchunk.sh` |
| `NPERM` | `experiments/05_learned_regulatory_sets/att/att_common_oracle.py`<br>`experiments/05_learned_regulatory_sets/att/att_common_oracle_hq.py`<br>`experiments/05_learned_regulatory_sets/att/ceiling_cand2.py`<br>and 4 other scripts |
| `NREC` | `experiments/02_association_and_finemapping/scripts/fm_region.sh` |
| `NTHREADS` | `pipeline/phenotype_nullmodel`<br>`pipeline/phenotype_nullmodel/run_step1v3/step1.sh`<br>`pipeline/phenotype_nullmodel/run_step1v4/step1.sh` |
| `NTILES` | `experiments/06_haplotype_context/run_20260922T052403Z/code/watchdog.sh` |
| `NULL_MODEL_FILE` | `pipeline/phenotype_nullmodel`<br>`pipeline/phenotype_nullmodel/run_pilot/readmod.R` |
| `NWORK` | `experiments/05_learned_regulatory_sets/att/run_cache_hq.sh`<br>`experiments/05_learned_regulatory_sets/att/run_interact.sh`<br>`experiments/05_learned_regulatory_sets/att/run_oracle_hq.sh` |
| `N_ANNOTATION_VARIANTS` | `pipeline/annotation/run_annot/build_cache.py`<br>`pipeline/annotation/run_annot/l1_train.py`<br>`pipeline/annotation/run_annot/l1_train_v7.py` |
| `N_BBJ_VARIANTS` | `pipeline/annotation/run_annot_bbj/map_bbj.py`<br>`pipeline/annotation/run_annot_bbj/pip_bbj_A.py`<br>`pipeline/annotation/run_annot_bbj/rsq_bbj.py` |
| `N_F46` | `experiments/05_learned_regulatory_sets/f3f6/f3f6_prep.py` |
| `N_GENES` | `experiments/07_attention_based_variant_interactions/models/attention/bench_cpu_gpu.py`<br>`experiments/07_attention_based_variant_interactions/models/attention/bench_cpu_gpu2.py` |
| `N_MATRIX_ROWS` | `experiments/05_learned_regulatory_sets/f_cmatrix_v1.py` |
| `N_SAMPLES` | `experiments/01_annotation_weighting_models/bagg.py`<br>`experiments/01_annotation_weighting_models/bagg2.py`<br>`experiments/01_annotation_weighting_models/bagg2b.py`<br>and 42 other scripts |
| `N_TOTAL` | `experiments/07_attention_based_variant_interactions/data_prep/hei_cap_rule.py` |
| `N_TRAIN` | `experiments/07_attention_based_variant_interactions/data_prep/hei_cap_rule.py` |
| `OMP_NUM_THREADS` | `experiments/01_annotation_weighting_models/l1_train_v7.py`<br>`experiments/01_annotation_weighting_models/l1_train_v72.py`<br>`experiments/05_learned_regulatory_sets/f1f2/f1_tendency.py`<br>and 2 other scripts |
| `OPENBLAS_NUM_THREADS` | `experiments/05_learned_regulatory_sets/f1f2/f1_tendency.py` |
| `OUTD` | `experiments/04_context_feasibility_gates/build_gate1_server.py` |
| `OUTPUT_DIR` | `pipeline/phenotype_nullmodel`<br>`pipeline/phenotype_nullmodel/run_pheno2/build.py`<br>`pipeline/phenotype_nullmodel/run_pheno3/build.py` |
| `PART` | `pipeline/phenotype_nullmodel/run_pilot/pr.sh` |
| `PAT` | `experiments/05_learned_regulatory_sets/f3f6/wd_f3f6.sh` |
| `PCT` | `pipeline/annotation/run_annot/gwas_common.sh` |
| `PERL5LIB` | `pipeline/groupfiles_variants/run_g5/g5b.sh`<br>`pipeline/groupfiles_variants/run_g5/g5c.sh`<br>`pipeline/groupfiles_variants/run_g5/g5c2.sh` |
| `PFX` | `pipeline/groupfiles_variants/run_g5/g5c2.sh` |
| `PHENOTYPE_DATA_ROOT` | `experiments/07_attention_based_variant_interactions/data_prep/build_height.py`<br>`experiments/07_attention_based_variant_interactions/data_prep/peek_height.py` |
| `PHENO_ADDED_COHORT_MIN_DM_CASES` | `pipeline/phenotype_nullmodel` |
| `PHENO_ADDED_COHORT_MIN_HTN_CASES` | `pipeline/phenotype_nullmodel` |
| `PHENO_ADDED_COHORT_MIN_LIP_CASES` | `pipeline/phenotype_nullmodel` |
| `PHENO_CENSUS_OUT` | `pipeline/groupfiles_variants/run_cur/census.py` |
| `PHENO_COVARIATE_COLUMNS` | `experiments/07_attention_based_variant_interactions/data_prep/hei_prep_model.py`<br>`experiments/07_attention_based_variant_interactions/data_prep/mk_bench_base.py`<br>`experiments/07_attention_based_variant_interactions/data_prep/prep_model.py`<br>and 4 other scripts |
| `PHENO_CUR_OUT` | `pipeline/groupfiles_variants/run_cur2/build.py` |
| `PHENO_DIR` | `experiments/01_annotation_weighting_models/model_v8_out/build_resid_v8.py`<br>`experiments/03_audits_and_null_evaluation/audit.py`<br>`experiments/03_audits_and_null_evaluation/ca.py`<br>and 15 other scripts |
| `PHENO_V1_DM_CASES` | `pipeline/phenotype_nullmodel` |
| `PHENO_V1_HTN_CASES` | `pipeline/phenotype_nullmodel` |
| `PHENO_V1_LIP_CASES` | `pipeline/phenotype_nullmodel` |
| `PHENO_V2_DIR` | `pipeline/phenotype_nullmodel`<br>`pipeline/phenotype_nullmodel/run_pilot/agetest.py`<br>`pipeline/phenotype_nullmodel/run_step1v2/step1.sh` |
| `PHENO_V2_DM_CASES` | `pipeline/phenotype_nullmodel` |
| `PHENO_V2_HTN_CASES` | `pipeline/phenotype_nullmodel` |
| `PHENO_V2_LIP_CASES` | `pipeline/phenotype_nullmodel` |
| `PHENO_V3_DIR` | `pipeline/phenotype_nullmodel`<br>`pipeline/phenotype_nullmodel/run_step1v3/step1.sh`<br>`pipeline/phenotype_nullmodel/run_step1v4/step1.sh` |
| `PLINK2` | `experiments/03_audits_and_null_evaluation/vrpool.py`<br>`pipeline/groupfiles_variants/run_ldb/ldb1.sh`<br>`pipeline/groupfiles_variants/run_vr/cmp.sh`<br>and 12 other scripts |
| `PLINK2_BIN` | `experiments/07_attention_based_variant_interactions/baselines/prs_cs/extract_hm3.sh`<br>`experiments/07_attention_based_variant_interactions/baselines/prs_cs/extract_hm3_v2.sh`<br>`experiments/07_attention_based_variant_interactions/baselines/prs_cs/score_chr.sh`<br>and 3 other scripts |
| `PLINK_PREFIX` | `pipeline/phenotype_nullmodel`<br>`pipeline/phenotype_nullmodel/run_step1v2/step1.sh`<br>`pipeline/phenotype_nullmodel/run_step1v3/step1.sh`<br>and 1 other scripts |
| `POOL` | `experiments/03_audits_and_null_evaluation/single_typed15.sh`<br>`experiments/03_audits_and_null_evaluation/single_typed15_v2.sh`<br>`experiments/05_learned_regulatory_sets/settest/scripts/s4_assoc.sh`<br>and 6 other scripts |
| `POS_DELTAS` | `experiments/05_learned_regulatory_sets/att/pos_control.py` |
| `POS_NCAUSAL` | `experiments/05_learned_regulatory_sets/att/pos_control.py` |
| `POS_SCENARIOS` | `experiments/05_learned_regulatory_sets/att/pos_control.py` |
| `PREFILTER_OMITTED_COUNT` | `utils/reusable/build_chr22_common_gene_mask_v2.py` |
| `PROJECT_ROOT` | `experiments/01_annotation_weighting_models/acatv_axis.py`<br>`experiments/01_annotation_weighting_models/assemble.py`<br>`experiments/01_annotation_weighting_models/bagg.py`<br>and 619 other scripts |
| `PRS_CS_ROOT` | `experiments/07_attention_based_variant_interactions/baselines/prs_cs/prscs_chr.sh` |
| `PRUNED` | `pipeline/groupfiles_variants/run_vr/vr2.sh` |
| `PRUNED_PLINK_PREFIX` | `pipeline/groupfiles_variants/run_vr/cmp.sh`<br>`pipeline/groupfiles_variants/run_vr/mac.sh`<br>`pipeline/groupfiles_variants/run_vr/merge_only.sh`<br>and 6 other scripts |
| `PYTHON` | `experiments/03_audits_and_null_evaluation/run_ident_r2_chr22_v2.sh`<br>`experiments/03_audits_and_null_evaluation/run_typed_windows15.sh`<br>`pipeline/groupfiles_variants/run_b6/run_all_b6.sh`<br>and 3 other scripts |
| `PYTHONDONTWRITEBYTECODE` | `experiments/01_annotation_weighting_models/model_v10_out/l1_train_v10.py`<br>`experiments/01_annotation_weighting_models/model_v8_out/l1_train_v8.py`<br>`experiments/01_annotation_weighting_models/model_v9_out/l1_train_v9.py`<br>and 6 other scripts |
| `PYTHON_BIN` | `experiments/01_annotation_weighting_models/model_v10_out/run_model_v10.sh`<br>`experiments/07_attention_based_variant_interactions/baselines/prs_cs/prscs_chr.sh`<br>`experiments/07_attention_based_variant_interactions/infra/queue_p5.sh`<br>and 27 other scripts |
| `R2_MIN` | `experiments/05_learned_regulatory_sets/att/cache_hq.py`<br>`experiments/05_learned_regulatory_sets/att/run_cache_hq.sh` |
| `REF37` | `utils/reusable/build_r2_only_broad_liftover_input.sh` |
| `REF38` | `utils/reusable/build_r2_only_broad_liftover_input.sh` |
| `REFERENCE_DIR` | `pipeline/phenotype_nullmodel`<br>`pipeline/phenotype_nullmodel/run_pilot/bchunk.sh`<br>`pipeline/phenotype_nullmodel/run_pilot/pr.sh`<br>and 3 other scripts |
| `REJECT_N` | `utils/reusable/build_r2_only_broad_liftover_input.sh` |
| `REPORT` | `utils/reusable/run_with_resource_guard.sh` |
| `ROOT` | `utils/reusable/run_with_resource_guard.sh` |
| `RSS` | `experiments/05_learned_regulatory_sets/att/night_monitor.sh`<br>`experiments/05_learned_regulatory_sets/att/wd_fit.sh`<br>`experiments/05_learned_regulatory_sets/att/wd_gen.sh`<br>and 3 other scripts |
| `RSS_BUDGET_GB` | `pipeline/annotation/run_annot/wd5.sh` |
| `RUN_DIR` | `pipeline/phenotype_nullmodel`<br>`pipeline/phenotype_nullmodel/run_pilot/bchunk.sh`<br>`pipeline/phenotype_nullmodel/run_pilot/pr.sh`<br>and 3 other scripts |
| `SAIGE_IMAGE` | `experiments/01_annotation_weighting_models/run_arm_seq.sh`<br>`experiments/01_annotation_weighting_models/run_chunked.sh`<br>`experiments/02_association_and_finemapping/LDLR2nd_v4/v4.sh`<br>and 24 other scripts |
| `SAIGE_PHENOTYPE_DIR` | `experiments/03_audits_and_null_evaluation/prior.py` |
| `SAIGE_RUNNER` | `pipeline/phenotype_nullmodel`<br>`pipeline/phenotype_nullmodel/run_pilot/bchunk.sh`<br>`pipeline/phenotype_nullmodel/run_pilot/pr.sh`<br>and 4 other scripts |
| `SAIGE_STEP1_SCRIPT` | `pipeline/phenotype_nullmodel`<br>`pipeline/phenotype_nullmodel/run_step1v2/step1.sh`<br>`pipeline/phenotype_nullmodel/run_step1v3/step1.sh`<br>and 1 other scripts |
| `SAIGE_STEP2_SCRIPT` | `pipeline/phenotype_nullmodel`<br>`pipeline/phenotype_nullmodel/run_pilot/bchunk.sh`<br>`pipeline/phenotype_nullmodel/run_pilot/pr.sh`<br>and 1 other scripts |
| `SAMTOOLS` | `pipeline/groupfiles_variants/run_g5/g5c2.sh` |
| `SOURCE_DUPLICATE_KEY_N` | `utils/reusable/build_r2_only_broad_liftover_input.sh` |
| `SOURCE_DUP_STATUS` | `utils/reusable/build_r2_only_broad_liftover_input.sh` |
| `SOURCE_INPUT_N` | `utils/reusable/build_r2_only_broad_liftover_input.sh` |
| `SOURCE_INPUT_STATUS` | `utils/reusable/build_r2_only_broad_liftover_input.sh` |
| `SOURCE_NORMALIZED_N` | `utils/reusable/build_r2_only_broad_liftover_input.sh` |
| `SOURCE_NORMALIZED_STATUS` | `utils/reusable/build_r2_only_broad_liftover_input.sh` |
| `SOURCE_NORM_PATH` | `utils/reusable/build_r2_only_broad_liftover_input.sh` |
| `SOURCE_REF_MISMATCH_N` | `utils/reusable/build_r2_only_broad_liftover_input.sh` |
| `SOURCE_REF_STATUS` | `utils/reusable/build_r2_only_broad_liftover_input.sh` |
| `SPARSE_GRM_PATH` | `experiments/07_attention_based_variant_interactions/data_prep/prep_split.py`<br>`experiments/07_attention_based_variant_interactions/infra/step1_hei.sh`<br>`experiments/07_attention_based_variant_interactions/infra/step1_hei_v2.sh` |
| `ST` | `pipeline/groupfiles_variants/run_ldb/ldb1.sh` |
| `TABIX` | `pipeline/groupfiles_variants/run_b6/b6_card_join.py`<br>`pipeline/groupfiles_variants/run_g5/g5c2.sh` |
| `TAG` | `experiments/01_annotation_weighting_models/run_alpha_chain.sh`<br>`experiments/01_annotation_weighting_models/run_alpha_chain2.sh`<br>`experiments/01_annotation_weighting_models/run_alpha_chain3.sh`<br>and 7 other scripts |
| `TARGET_DUP_STATUS` | `utils/reusable/build_r2_only_broad_liftover_input.sh` |
| `TARGET_NORMALIZED_N` | `utils/reusable/build_r2_only_broad_liftover_input.sh` |
| `TARGET_NORMALIZED_STATUS` | `utils/reusable/build_r2_only_broad_liftover_input.sh` |
| `TARGET_NORM_PATH` | `utils/reusable/build_r2_only_broad_liftover_input.sh` |
| `TARGET_REF_MISMATCH_N` | `utils/reusable/build_r2_only_broad_liftover_input.sh` |
| `TARGET_REF_STATUS` | `utils/reusable/build_r2_only_broad_liftover_input.sh` |
| `TCHL_MAX_VALUE` | `pipeline/phenotype_nullmodel`<br>`pipeline/phenotype_nullmodel/run_pheno3/build.py` |
| `TCHL_MIN_VALUE` | `pipeline/phenotype_nullmodel`<br>`pipeline/phenotype_nullmodel/run_pheno3/build.py` |
| `TORCH_HOME` | `experiments/01_annotation_weighting_models/model_v10_out/l1_train_v10.py`<br>`experiments/01_annotation_weighting_models/model_v8_out/l1_train_v8.py`<br>`experiments/01_annotation_weighting_models/model_v9_out/l1_train_v9.py` |
| `TRAIT` | `experiments/05_learned_regulatory_sets/att/att_fit_trait.py`<br>`experiments/05_learned_regulatory_sets/att/ceiling_w.py`<br>`pipeline/phenotype_nullmodel/run_pilot/pr.sh` |
| `TRAITS` | `experiments/05_learned_regulatory_sets/settest/scripts/s4_assoc.sh` |
| `TRAIT_COLUMN_KEYWORDS_JSON` | `experiments/03_audits_and_null_evaluation/audit.py` |
| `TSS_POSITION` | `experiments/03_audits_and_null_evaluation/ap3.py`<br>`experiments/03_audits_and_null_evaluation/coordchk.py` |
| `UDOCKER_BIN` | `experiments/01_annotation_weighting_models/run_arm_seq.sh`<br>`experiments/01_annotation_weighting_models/run_chunked.sh`<br>`experiments/02_association_and_finemapping/LDLR2nd_v4/v4.sh`<br>and 24 other scripts |
| `UL` | `pipeline/annotation/run_annot/run_pool.sh` |
| `USER` | `utils/reusable/run_with_resource_guard.sh` |
| `VARIANCE_RATIO_PLINK_BASENAME` | `experiments/03_audits_and_null_evaluation/vrchk.py` |
| `WA` | `experiments/02_association_and_finemapping/scripts/fm_chr.sh` |
| `WAITED` | `experiments/01_annotation_weighting_models/patch_wait.py`<br>`experiments/01_annotation_weighting_models/run_arm_seq.sh` |
| `WANT` | `pipeline/groupfiles_variants/run_vr/hc.sh` |
| `WD_HEARTBEAT` | `experiments/01_annotation_weighting_models/bagg2.py`<br>`experiments/01_annotation_weighting_models/bagg2b.py`<br>`experiments/01_annotation_weighting_models/bagg2b_diag.py` |
| `WD_MAX_AGE` | `experiments/01_annotation_weighting_models/bagg2.py`<br>`experiments/01_annotation_weighting_models/bagg2b.py`<br>`experiments/01_annotation_weighting_models/bagg2b_diag.py` |
| `WGRID` | `experiments/05_learned_regulatory_sets/att/ceiling_cand.py`<br>`experiments/05_learned_regulatory_sets/att/ceiling_cand2.py`<br>`experiments/05_learned_regulatory_sets/att/ceiling_w.py` |
| `XDG_CACHE_HOME` | `experiments/01_annotation_weighting_models/model_v10_out/l1_train_v10.py`<br>`experiments/01_annotation_weighting_models/model_v8_out/l1_train_v8.py`<br>`experiments/01_annotation_weighting_models/model_v9_out/l1_train_v9.py` |
| `XP` | `pipeline/groupfiles_variants/run_g5/g5p.sh` |
