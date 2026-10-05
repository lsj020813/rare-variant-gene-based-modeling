# experiments/07_attention_based_variant_interactions

111 scripts. PRS-CS source is not bundled; set `PRS_CS_ROOT` to an external official checkout.

Shell launchers resolve code from their location and keep working data in `${PROJECT_ROOT}/work/prs`. Tools use PATH, with optional `PYTHON_BIN`, `BCFTOOLS_BIN`, `PLINK2_BIN`, and `UDOCKER_BIN` overrides.

Phenotype preparation reads `PHENOTYPE_DATA_ROOT`, `HEIGHT_SOURCES_JSON` (a JSON list of relative path/column pairs), `COHORT_INDICATOR_COLUMNS_JSON` (two indicator columns), and `PHENO_COVARIATE_COLUMNS` (comma separated covariates in the intended order). The height inspection utility also reads `HEIGHT_SCAN_GLOB`. Null-model scripts require `SPARSE_GRM_PATH` and `SAIGE_IMAGE`.

## baselines/prs_cs

- `baselines/prs_cs/dl_ldref.sh`
- `baselines/prs_cs/extract_hm3.sh`
- `baselines/prs_cs/extract_hm3_v2.sh`
- `baselines/prs_cs/mk_bim.py`
- `baselines/prs_cs/prep_bbj.py`
- `baselines/prs_cs/prep_bbj_hei.py`
- `baselines/prs_cs/prep_bbj_ht.py`
- `baselines/prs_cs/prscs_chr.sh`
- `baselines/prs_cs/qc_prs.py`
- `baselines/prs_cs/score_chr.sh`

## models/linear

- `models/linear/hei_linear.py`
- `models/linear/hei_linear_v2.py`
- `models/linear/hei_seed_linear.py`
- `models/linear/hei_seed_p3oof.py`
- `models/linear/p3_maf_strata.py`
- `models/linear/train_linear.py`
- `models/linear/train_linear_v2.py`
- `models/linear/train_linear_v3.py`

## models/attention

- `models/attention/bench_cpu_gpu.py`
- `models/attention/bench_cpu_gpu2.py`
- `models/attention/bench_scale.sh`
- `models/attention/test_fastpath.py`
- `models/attention/test_fastpath_gpu.py`
- `models/attention/train_att.py`
- `models/attention/train_att_hei.py`
- `models/attention/train_att_hei_res.py`
- `models/attention/train_att_hei_res_v2.py`
- `models/attention/train_att_hei_res_v3.py`
- `models/attention/train_att_hei_res_v4.py`
- `models/attention/train_att_hei_res_v5.py`
- `models/attention/train_att_hei_res_v6.py`
- `models/attention/train_att_hei_seed.py`
- `models/attention/train_att_hei_v2.py`
- `models/attention/train_att_v3.py`
- `models/attention/train_att_v4.py`
- `models/attention/train_att_v5.py`
- `models/attention/train_att_v6.py`

## data_prep

- `data_prep/build_height.py`
- `data_prep/build_locus.py`
- `data_prep/f_features_hei.py`
- `data_prep/f_features_hei_seed.py`
- `data_prep/f_features_prs.py`
- `data_prep/hei_cap_calc.py`
- `data_prep/hei_cap_rule.py`
- `data_prep/hei_clump_select.py`
- `data_prep/hei_clump_select_v2.py`
- `data_prep/hei_clump_select_v3.py`
- `data_prep/hei_extract.py`
- `data_prep/hei_extract_gt_merge.py`
- `data_prep/hei_extract_gt_part.py`
- `data_prep/hei_extract_merge.py`
- `data_prep/hei_extract_part.py`
- `data_prep/hei_locus_counts.py`
- `data_prep/hei_prep_model.py`
- `data_prep/hei_rint.py`
- `data_prep/hei_seed_extract_gt_merge.py`
- `data_prep/hei_seed_extract_gt_part.py`
- `data_prep/hei_seed_extract_merge.py`
- `data_prep/hei_seed_extract_part.py`
- `data_prep/hei_seed_smoke_check.py`
- `data_prep/hei_select.py`
- `data_prep/hei_select_seed.py`
- `data_prep/hei_select_seed_r3.py`
- `data_prep/hei_select_seed_r4.py`
- `data_prep/hei_split.py`
- `data_prep/hei_tissue_scan.py`
- `data_prep/hei_tissue_scan_v2.py`
- `data_prep/hei_tissue_scan_v3.py`
- `data_prep/hei_token_count.py`
- `data_prep/mk_bench_base.py`
- `data_prep/peek_height.py`
- `data_prep/peek_meta.py`
- `data_prep/prep_model.py`
- `data_prep/prep_model_v2.py`
- `data_prep/prep_split.py`
- `data_prep/prs_common.py`

## infra

- `infra/queue_p5.sh`
- `infra/queue_p5_v2.sh`
- `infra/run_build_loci.sh`
- `infra/run_chain1.sh`
- `infra/run_extract_hm3.sh`
- `infra/run_extract_hm3_b.sh`
- `infra/run_extract_hm3_v2.sh`
- `infra/run_feat.sh`
- `infra/run_feat_hei.sh`
- `infra/run_feat_hei_seed.sh`
- `infra/run_hei_att.sh`
- `infra/run_hei_att_v2.sh`
- `infra/run_hei_att_v3.sh`
- `infra/run_hei_att_v4.sh`
- `infra/run_hei_chain.sh`
- `infra/run_hei_chain_v2.sh`
- `infra/run_hei_chain_v3.sh`
- `infra/run_hei_extract.sh`
- `infra/run_hei_extract_gt.sh`
- `infra/run_hei_extract_v2.sh`
- `infra/run_hei_prscs.sh`
- `infra/run_hei_res_att.sh`
- `infra/run_hei_res_att_v2.sh`
- `infra/run_hei_res_att_v3.sh`
- `infra/run_hei_res_att_v4.sh`
- `infra/run_hei_res_att_v5.sh`
- `infra/run_hei_res_att_v6.sh`
- `infra/run_hei_seed_att.sh`
- `infra/run_hei_seed_extract.sh`
- `infra/run_hei_seed_lin.sh`
- `infra/run_prscs_bbj.sh`
- `infra/step1_hei.sh`
- `infra/step1_hei_v2.sh`
- `infra/wd_prs.sh`
- `infra/wd_prs_big.sh`
