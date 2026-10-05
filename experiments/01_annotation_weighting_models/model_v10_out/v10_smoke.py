import argparse
import ast
import json
from pathlib import Path
from types import SimpleNamespace

def run_smoke(v, args, guard):
    np = v.np
    cfg, manifest = v.fixture_args(args, args.smoke_suite)
    cfg.require_convergence = True
    out = v.owned_path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    v.require(not (out / 'SMOKE_SUITE_PASS').exists(), 'smoke output already completed; use another directory')
    v.atomic_json(out / 'SMOKE_SETTINGS.json', dict(manifest=manifest, effective_settings=vars(cfg), guard=guard, code_sha256=v.sha_file(v.__file__)))
    data = v.Data(cfg)
    checks = {}

    def check(name, condition, details=None):
        checks[name] = dict(pass_=bool(condition), details=details)
        v.atomic_json(out / 'SMOKE_CHECKS.partial.json', checks)
        v.require(condition, 'smoke ' + name)

    def must_fail(name, action):
        try:
            action()
        except RuntimeError as error:
            check(name, str(error).startswith('GATE FAIL:'))
        else:
            check(name, False)

    check('fixture_dimensions', len(data.keys) == 20000 and len(data.traits) == 3 and
          len(set(data.regions)) == 20 and len(data.cs_blocks) == 40,
          dict(variants=len(data.keys), traits=len(data.traits), regions=len(set(data.regions)), cs=len(data.cs_blocks)))
    check('unlabelled_not_zero_imputed', len(data.y) == 54000 and data.audit['annotation_unlabelled'] == 2000)
    check('hard_resource_limits', guard['threads'] <= 4 and guard['memory_cap_bytes'] <= 8 * 1024 ** 3)

    tiny = SimpleNamespace(ti=np.array([0, 0, 0, 0]), regions=np.array(['r'] * 4),
                           cs=np.array([1, 1, -1, -1]), y=np.array([.2, .02, .5, .001]))
    w, wm = v.weights(tiny, np.arange(4))
    check('B2_weights_hand_calculation', np.allclose(w, [1 / 3, 1 / 3, 2 / 3, 2 / 3]) and
          np.allclose(wm['class_final_mass'], [1, 1]))
    yy = np.array([.2, .5, .5, .001])
    changed, _ = v.weights(tiny, np.arange(4), yy)
    check('B2_reweights_shuffled_continuous_labels', np.allclose(changed, [.25, .25, .5, 1]))

    check('CE_independent_formula', np.allclose(v.ce_values(np.array([.25, .8]), np.array([.3, .9])),
          [-.3 * np.log(.25) - .7 * np.log(.75), -.9 * np.log(.8) - .1 * np.log(.2)]))
    check('CE_endpoint_finite', np.isfinite(v.ce_values(np.array([0., 1.]), np.array([1., 0.]))).all())

    splits = v.split_data(cfg, data)
    check('sealed_chr_excluded', not np.isin(data.chr, ['21', '22']).any())
    check('two_axis_holdout', not np.isin(data.groups[splits['train']], cfg.holdout_trait_group).any() and
          not set(data.chr[splits['train']]) & set(data.chr[splits['chromosome']]) and
          not set(data.chr[splits['inner_train']]) & set(data.chr[splits['inner_valid']]))
    sealed = out / 'sealed_contract.tsv'
    v.atomic_tsv(sealed, ['chr', 'pip'], [dict(chr='21', pip='DO_NOT_PARSE')])
    must_fail('sealed_row_rejected_before_PIP_parse', lambda: list(v.tsv_rows(sealed, ['chr', 'pip'])))
    must_fail('write_outside_workspace_rejected', lambda: v.owned_path('/tmp/v10_forbidden'))
    must_fail('key_position_mismatch_rejected', lambda: v.variant_key('1:100:A:C', '1', 101))
    must_fail('float_key_position_rejected', lambda: v.variant_key('1:100.0:A:C', '1', 100))

    t = 0
    name = data.traits[t]['trait']
    rows = list(v.tsv_rows(Path(cfg.bbj_marginal) / (name + '.tsv.gz'),
                           ['variant_hg19', 'region', 'beta_marginal', 'se']))
    first_region = rows[0]['region']
    rows = [r for r in rows if r['region'] == first_region]
    z = np.array([abs(float(r['beta_marginal']) / float(r['se'])) for r in rows])
    expected = np.array([(np.sum(z < value) + (np.sum(z == value) + 1) / 2) / len(z) for value in z])
    key_to_rank = {r['variant_hg19']: rr for r, rr in zip(rows, expected)}
    idx = np.flatnonzero((data.ti == t) & (data.regions == first_region))
    check('C2_independent_region_and_rank', np.allclose(data.C2[idx, 2], len(z)) and
          np.allclose(data.C2[idx, 3], [key_to_rank[data.keys[data.vi[i]]] for i in idx]) and
          np.allclose(data.C2[idx, :2], np.column_stack([data.maf[data.vi[idx]], data.rsq[data.vi[idx]]])))

    c1tiny = SimpleNamespace(keys=np.array(['a', 'b', 'c', 'd']),
                             col=lambda _: np.array([1., 1., 3., np.nan]))
    ecdf, _ = v.staar_pair_weights(c1tiny, 'synthetic', np.arange(4))
    check('C1_ECDF_ties_missing_reference', np.allclose(ecdf, [.375, .375, .75, .5]))
    probabilities = np.array([[.2, .6], [.4, .8]])
    direct = .5 - np.arctan(np.tan((.5 - probabilities) * np.pi).mean(axis=0)) / np.pi
    check('C1_Cauchy_independent_formula', np.allclose(v.cauchy_combine(probabilities), direct, atol=1e-14))
    normal_scores = np.array([[1., 2.], [2., 3.]])
    acat_p = v.cauchy_combine(2 * v.norm.sf(abs(normal_scores)))
    check('ACAT_equivalent_z_consistency', np.allclose(2 * v.norm.sf(v.acat_equivalent_z(normal_scores)), acat_p))

    plan = v.PermutationPlan(data, splits['chromosome'], cfg.null_perm, cfg.test_seed)
    shuffled = next(plan.labels(data.y))
    original = data.y[splits['chromosome']]
    check('permutation_preserves_each_trait_region', all(np.array_equal(np.sort(shuffled[b]), np.sort(original[b])) for b in plan.blocks))
    check('same_permutation_plan_reproducible', np.array_equal(shuffled, next(plan.labels(data.y))))
    check('permutation_not_identity', not np.array_equal(shuffled, original))

    candidate = dict(strong=dict(score=.4), weak=dict(score=1.))
    chosen, _ = v.select_one_se(candidate, dict(strong=np.zeros(4), weak=np.array([-1., 1., -1., 1.])))
    check('one_SE_uses_paired_SD_not_SD_over_sqrt_B', chosen == 'strong')
    rank, counts, eligible = v.cs_spearman(np.array([[1., 2., 3.], [1., 1., 1.]]),
                                         np.array([3., 2., 1.]), [np.arange(3)])
    check('CS_rank_and_constant_undefined', np.allclose(rank[0], -1) and np.isnan(rank[1]) and counts == [1, 0] and eligible == 1)

    design = v.make_design(cfg, data, splits['inner_train'])
    trainer = v.Trainer(cfg, data, design, splits['inner_train'])
    rng = np.random.default_rng(cfg.inner_seed)
    a = trainer.initial + rng.normal(scale=.02, size=len(trainer.initial))
    direction = rng.normal(size=len(a)); direction /= np.linalg.norm(direction)
    lam, h = .01, 1e-5
    _, gradient = trainer.objective(a, lam)
    numerical = (trainer.objective(a + h * direction, lam)[0] - trainer.objective(a - h * direction, lam)[0]) / (2 * h)
    check('soft_PIP_penalty_gradient', np.isclose(numerical, gradient @ direction, rtol=1e-5, atol=1e-7))
    modified = data.X.copy()
    held_variants = np.setdiff1d(np.arange(len(data.keys)), np.unique(data.vi[splits['inner_train']]))
    modified[held_variants, data.cols.index('cons')] += 100.
    alt = v.Design(cfg, modified, data.cols, np.unique(data.vi[splits['inner_train']]), cfg.phi_columns)
    check('heldout_annotations_do_not_fit_design', design.signature == alt.signature and np.array_equal(design.P, alt.P))
    legacy = v.Design(cfg, data.X, data.cols, np.unique(data.vi[splits['inner_train']]), 'v9-eighteen')
    check('explicit_feature_modes', len(design.selected) == 34 and len(legacy.selected) == 18)
    group_cfg = argparse.Namespace(**vars(cfg))
    group_cfg.group_phi = True
    group_cfg.unseen_phi_group = '지질'
    group_cfg.export_phi_group = '지질'
    group_trainer = v.Trainer(group_cfg, data, design, splits['inner_train'])
    ga = group_trainer.initial.copy()
    ga[0] = 1.; ga[group_trainer.q] = -1.
    gm = v.export_model(group_cfg, data, group_trainer, ga)
    gi = splits['joint']
    check('group_phi_transfer_export', np.allclose(group_trainer.predict(ga, gi),
          v.PhiPredictor(gm).predict(data.X[data.vi[gi]], data.cols), atol=1e-12))
    group_trainer.args = argparse.Namespace(**vars(group_cfg))
    group_trainer.args.unseen_phi_group = None
    must_fail('group_phi_unsealed_transfer_rejected', lambda: group_trainer.predict(ga, gi))
    c1_cfg = argparse.Namespace(**vars(cfg)); c1_cfg.c1_columns = 'all34'
    pp, _ = v.c1_predictions(c1_cfg, data, splits['train'], splits['chromosome'])
    check('C1_all34_option', len(pp) == 35 and all(np.isfinite(p).all() for p in pp.values()))
    del group_trainer, pp
    del alt, legacy, trainer, design

    result, exported, final_trainer, a = v.run_pipeline(cfg, data, out)
    check('selected_and_final_fits_converged', all(result['final_fit'][arm]['success'] for arm in ('C2', 'phi')) and
          all(c['optimizer']['success'] for arm in ('C2', 'phi') for c in result['selection'][arm]['candidates'].values()))
    idx = splits['chromosome']
    true_signal = np.load(Path(args.smoke_suite) / 'oracle_signal.npy', allow_pickle=False)
    fitted_signal = final_trainer.basis_for(idx) @ a[:final_trainer.ncoef]
    correlation = float(np.corrcoef(fitted_signal, true_signal[data.vi[idx]])[0, 1])
    check('oracle_phi_recovery', correlation > manifest['smoke_gates']['oracle_correlation_min'], dict(correlation=correlation, minimum=.9))
    predictor = v.PhiPredictor(exported)
    known_trait = data.traits[final_trainer.trait_ids[0]]['trait']
    known_t = final_trainer.trait_ids[0]
    selected = idx[data.ti[idx] == known_t]
    expected = final_trainer.predict(a, selected)
    restored = predictor.predict(data.X[data.vi[selected]], data.cols, trait=known_trait)
    check('export_reconstruction', np.allclose(expected, restored, rtol=1e-12, atol=1e-12))
    row = {c: data.X[data.vi[selected[0]], j] for j, c in enumerate(data.cols)}
    row = {c: x if np.isfinite(x) else 'NA' for c, x in row.items()}
    check('apply_phi_row_interface', np.isclose(v.apply_phi(row, exported, trait=known_trait), expected[0]))
    unbound = dict(exported, default_intercept=None)
    must_fail('unsealed_transfer_intercept_rejected', lambda: v.PhiPredictor(unbound).predict(data.X[:1], data.cols))
    check('band_score_file', (out / 'phi_scores_1-5.tsv.gz').is_file() and result['band_scores']['1-5']['n'] == 5000)
    lc = result['learning_curves']
    check('learning_curves_partial_and_complete', len(lc['rows']) == 7 and
          (out / 'learning_curve.partial.json').is_file() and (out / 'learning_curve.json').is_file() and
          all((out / ('curve_' + str(i)) / 'POINT_DONE').is_file() for i in range(7)))
    check('learning_curves_fixed_design_and_permutations',
          len({row['phi_design_sha256'] for row in lc['rows']}) == 1 and
          len({row['C2_design_sha256'] for row in lc['rows']}) == 1 and
          len({row['evaluation']['null_plan']['sha256'] for row in lc['rows']}) == 1)

    reps = manifest['smoke_gates']['null_replicates']
    primary_idx = splits['chromosome']
    frozen = final_trainer.predict(a, primary_idx)
    c2_trainer, c2_a = final_trainer.comparator
    c2_pred = c2_trainer.predict(c2_a, primary_idx)
    zscores = []
    for r in range(reps):
        seeds = np.random.SeedSequence([cfg.test_seed, r]).generate_state(2)
        perm = v.PermutationPlan(data, primary_idx, 1, int(seeds[0]))
        yy = data.y.copy(); yy[primary_idx] = next(perm.labels(data.y))
        pp = dict(flat=np.full(len(primary_idx), data.y[splits['train']].mean()), C2=c2_pred, phi=frozen)
        metrics, _ = v.evaluate(cfg, data, primary_idx, pp, int(seeds[1]), yy, only_bands=[v.PRIMARY_BAND])
        zscores.append(metrics['bands'][v.PRIMARY_BAND]['models']['phi']['delta_ce_c2']['z'])
    finite = [z for z in zscores if z is not None and np.isfinite(z)]
    check('null_shuffle_z_centering', len(finite) == reps and abs(float(np.mean(finite))) <= manifest['smoke_gates']['null_z_mean_tolerance'],
          dict(mean=float(np.mean(finite)), sd=float(np.std(finite, ddof=1)), n=len(finite), tolerance=manifest['smoke_gates']['null_z_mean_tolerance']))

    simulation = v.simulate(cfg, data, out)
    check('D3_three_modes_and_atomic_partials', all(simulation[mode] for mode in ('null', 'oracle', 'power')) and
          (out / 'simulation.partial.json').is_file() and (out / 'simulation.json').is_file())
    result['simulation'] = simulation
    result['verdict'] = v.classify(cfg, result, simulation)
    null_result = json.loads(json.dumps(result))
    null_result['evaluation']['chromosome']['bands'][v.PRIMARY_BAND]['ladder']['C1->phi'].update(observed=0., p_greater=1., status='ok')
    check('low_power_not_misclassified', v.classify(cfg, null_result, dict(power_adequate=False))['category'] == '판정 불능')
    significant_result = json.loads(json.dumps(null_result))
    significant_result['evaluation']['chromosome']['bands'][v.PRIMARY_BAND]['ladder']['C1->phi'].update(observed=1., p_greater=0.)
    check('significant_improvement_success', v.classify(cfg, significant_result)['category'] == '성공')
    for model in null_result['evaluation']['chromosome']['bands'][v.PRIMARY_BAND]['models'].values():
        model['weighted_ce']['observed'] = .5
    check('powered_all_flat_meaningful_failure', v.classify(cfg, null_result, dict(power_adequate=True))['category'] == '의미 있는 실패')
    v.verify_input_fingerprints(data)
    check('no_deprecated_classes_or_imports', static_purity(v))
    v.atomic_json(out / 'RESULT.json', dict(**result, guard=guard, smoke_checks=checks))
    v.atomic_json(out / 'SMOKE_SUITE_PASS', dict(synthetic=True, checks=len(checks),
                  result_sha256=v.sha_file(out / 'RESULT.json'), oracle_correlation=correlation,
                  null_z_mean=float(np.mean(finite)), code_sha256=v.sha_file(v.__file__)))
    v.atomic_json(out / 'L1_DONE', dict(kind='synthetic_smoke_suite', result_sha256=v.sha_file(out / 'RESULT.json')))
    print('SMOKE_SUITE_PASS', flush=True)
    print('L1_DONE smoke_suite', flush=True)
    return checks

def static_purity(v):
    tree = ast.parse(Path(v.__file__).read_text())
    forbidden = {'Residuals', 'BBJTruth', 'Evaluation', 'mix_terms', 'mix_em', 'mixed_prior',
                 'u17_metrics', 'conditional_diagnostics', 'load_ds', 'load_conditioning'}
    definitions = {n.name for n in ast.walk(tree) if isinstance(n, (ast.ClassDef, ast.FunctionDef))}
    imports = [n.module or '' for n in ast.walk(tree) if isinstance(n, ast.ImportFrom)]
    options = [n.value for n in ast.walk(tree) if isinstance(n, ast.Constant) and isinstance(n.value, str) and n.value.startswith('--')]
    return not (definitions & forbidden) and not any('l1_train_v9' in x or 'l1_train_v8' in x for x in imports) and not any(
        name in options for name in ('--resid', '--root', '--offset-dir', '--bbj-gene-z', '--kappa-grid'))
