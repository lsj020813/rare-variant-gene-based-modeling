
import os
import platform
import sys
import time
import traceback

_HERE = os.path.dirname(os.path.abspath(__file__))
_CODE = os.path.abspath(os.path.join(_HERE, "..", "code"))
if _CODE not in sys.path:
    sys.path.insert(0, _CODE)

import numpy as np

from hcf import align, phase, state
from hcf.align import AlignmentError, OrientationSourceError
from hcf.phase import PhaseSetCrossingError
from hcf.state import MISSING, OTHER, LabelLeakageError, StateDictionary

DESCRIPTIONS = {
    "test_01_cis_trans_sign": (
        "cis/trans 부호: [1,1]/[0,0] 과 [1,0]/[0,1] 은 DS/GT count 가 같고 Q 부호만 "
        "+1/-1 로 갈리며, 2C-GjGk 와 사본차 곱 두 산식이 일치한다"
    ),
    "test_02_copy_swap_invariance": (
        "copy-swap 불변: local window 전체에서 copy1/copy2 를 교환해도 상태 count Z·"
        "설계행렬·모든 pair 의 Q·선형 예측이 완전히 같다"
    ),
    "test_03_single_marker_subwindow": (
        "단일 마커 subwindow: 상태 count 가 G 만으로 결정되고 pair 가 0개여서 "
        "새 국소 위상 정보가 없음을 출력한다"
    ),
    "test_04_complete_ld": (
        "완전 LD: 두 마커가 항상 함께 움직이면 trans 사례가 0, 사본 r^2=1, 유효 rank 1 로 "
        "원인별 기여·phase-contrast 분리 제한을 탐지한다"
    ),
    "test_05_alt_major": (
        "ALT-major: INFO/AF>0.5 면 minor dosage = 2-DS 를 양쪽 사본·dosage 에 일관 적용하고, "
        "결측은 보존되며 Q 부호는 뒤집힌 마커 수에 따라 (-1)^flip 로만 변한다"
    ),
    "test_06_misalignment_fail_fast": (
        "sample/key misalignment: 표본 순서를 일부 섞거나 변이 키·행렬 모양이 어긋나면 "
        "조용히 재정렬하지 않고 AlignmentError 로 즉시 실패한다"
    ),
    "test_07_phase_set_crossing": (
        "phase-set crossing: 다른 PS·다른 염색체·다른 tile·20 kb 초과 gap 을 임의로 이은 "
        "subwindow 를 차단하고, PS 부재 시에만 가정 등급을 기록한다"
    ),
    "test_08_missing_not_reference": (
        "missingness: 결측 call 은 '.' 로 인코딩되어 all-reference state 로 세지지 않고 "
        "OTHER 로 가며, 결측이 섞인 Q·dosage 는 0 이 아니라 NA 다"
    ),
    "test_09_unseen_state_to_other": (
        "unseen state: B/C 에만 나타난 패턴은 A 사전의 OTHER 로 가고 사전은 확대되지 않으며 "
        "Z 행 합은 2 로 유지된다"
    ),
    "test_10_label_leakage_sentinel": (
        "label leakage sentinel: 표현형 열(y·TCHL·pheno*)이나 phenotype 인자가 상태·pair·"
        "방향 결정 함수에 들어오면 LabelLeakageError/OrientationSourceError 로 막는다"
    ),
    "test_11_constant_vector_na": (
        "constant vector: 상수 벡터의 correlation·표준화 effective rank 분모 0 은 NA(nan) 이며 "
        "r=1 완전 재현성으로 쓰지 않는다"
    ),
    "test_12_state_label_permutation": (
        "state label permutation: 상태 열 번호를 임의로 재배열해도 key 별 사람 count 와 "
        "선형 예측이 동일하다"
    ),
}

TEST_ORDER = [
    "test_01_cis_trans_sign",
    "test_02_copy_swap_invariance",
    "test_03_single_marker_subwindow",
    "test_04_complete_ld",
    "test_05_alt_major",
    "test_06_misalignment_fail_fast",
    "test_07_phase_set_crossing",
    "test_08_missing_not_reference",
    "test_09_unseen_state_to_other",
    "test_10_label_leakage_sentinel",
    "test_11_constant_vector_na",
    "test_12_state_label_permutation",
]

_DETAILS = {}

def _record(detail):
    _DETAILS[sys._getframe(1).f_code.co_name] = detail
    return None

def make_hap(rows):
    return np.array([[list(c1), list(c2)] for c1, c2 in rows], dtype=int)

def synthetic_window(n_people=60, n_markers=4, seed=20260922):
    rng = np.random.RandomState(seed)
    patterns = np.array(
        [
            [0, 0, 0, 0],
            [1, 1, 0, 0],
            [0, 1, 1, 0],
            [1, 0, 0, 1],
            [1, 1, 1, 1],
        ],
        dtype=int,
    )
    weights = np.array([0.40, 0.25, 0.15, 0.13, 0.07])
    idx = rng.choice(len(patterns), size=(n_people, 2), p=weights)
    hap = patterns[idx][:, :, :n_markers]
    return np.array(hap, dtype=int)

def marker_meta(n_markers, chrom="22", start=17000000, step=1000, ps=None, tile="22:170"):
    return [
        {
            "chrom": chrom,
            "pos": start + i * step,
            "ps": (ps[i] if isinstance(ps, (list, tuple)) else ps),
            "tile": tile,
        }
        for i in range(n_markers)
    ]

def assert_raises(exc_type, fn, *a, **kw):
    try:
        fn(*a, **kw)
    except exc_type:
        return True
    except Exception as other:
        raise AssertionError(
            "expected {0}, got {1}: {2}".format(exc_type.__name__, type(other).__name__, other)
        )
    raise AssertionError("expected {0}, nothing raised".format(exc_type.__name__))

def test_01_cis_trans_sign():
    hap = make_hap([([1, 1], [0, 0]), ([1, 0], [0, 1])])
    gj = phase.dosage(hap, 0)
    gk = phase.dosage(hap, 1)
    assert list(gj) == [1.0, 1.0], gj
    assert list(gk) == [1.0, 1.0], gk
    assert gj[0] == gj[1] and gk[0] == gk[1]

    q = phase.phase_contrast_q(hap, 0, 1, check_identity=True)
    assert q[0] == 1.0, q
    assert q[1] == -1.0, q
    assert q[0] == -q[1]

    c = phase.co_occurrence(hap, 0, 1)
    q_alt = 2.0 * c - gj * gk
    assert np.allclose(q, q_alt), (q, q_alt)

    cls = phase.phase_class(hap, 0, 1)
    assert list(cls) == ["cis", "trans"], list(cls)

    hom = make_hap([([1, 1], [0, 1])])
    assert phase.phase_contrast_q(hom, 0, 1)[0] == 0.0

    keys_cis = state.encode_subwindow(hap[:1])
    keys_trans = state.encode_subwindow(hap[1:])
    assert set(keys_cis[0]) == set(["11", "00"]), keys_cis
    assert set(keys_trans[0]) == set(["10", "01"]), keys_trans
    assert set(keys_cis[0]) != set(keys_trans[0])
    _record({
        "q_cis": float(q[0]),
        "q_trans": float(q[1]),
        "G_identical": True,
        "identity_2C_minus_GjGk": "verified",
    })

def test_02_copy_swap_invariance():
    hap = synthetic_window(n_people=80, n_markers=4, seed=20260922)
    swapped = phase.swap_copies(hap)
    assert not np.array_equal(hap, swapped), "합성 입력이 이미 대칭이면 시험 의미가 없다"

    keys = state.encode_subwindow(hap)
    keys_sw = state.encode_subwindow(swapped)
    sd = StateDictionary.fit(keys, min_distinct_people=5, max_states=32)
    Z, cols = sd.transform(keys)
    Z_sw, cols_sw = sd.transform(keys_sw)
    assert cols == cols_sw
    assert np.array_equal(Z, Z_sw), "state count 가 사본 교환에 의존한다"
    assert np.all(Z.sum(axis=1) == 2.0)

    X, xcols = sd.design_matrix(keys)
    X_sw, _ = sd.design_matrix(keys_sw)
    assert np.array_equal(X, X_sw)

    sd_sw = StateDictionary.fit(keys_sw, min_distinct_people=5, max_states=32)
    assert sd.states == sd_sw.states
    assert sd.reference_state == sd_sw.reference_state

    n_markers = hap.shape[2]
    for j in range(n_markers):
        for k in range(j + 1, n_markers):
            q = phase.phase_contrast_q(hap, j, k)
            q_sw = phase.phase_contrast_q(swapped, j, k)
            assert np.array_equal(q, q_sw), (j, k)
            sup = phase.pair_support(hap, j, k, min_double_het=1, min_cis=1, min_trans=1)
            sup_sw = phase.pair_support(swapped, j, k, min_double_het=1, min_cis=1, min_trans=1)
            assert sup == sup_sw, (j, k)

    coef = dict((c, 0.1 * (i + 1)) for i, c in enumerate(cols))
    pred = state.predict_linear(Z, cols, coef, intercept=0.5)
    pred_sw = state.predict_linear(Z_sw, cols_sw, coef, intercept=0.5)
    assert np.array_equal(pred, pred_sw)
    _record({
        "n_people": int(hap.shape[0]),
        "n_states": len(sd.states),
        "Z_identical": True,
        "Q_identical_all_pairs": True,
        "prediction_identical": True,
    })

def test_03_single_marker_subwindow():
    hap = make_hap(
        [([0], [0]), ([1], [0]), ([0], [1]), ([1], [1])] * 15
    )
    keys = state.encode_subwindow(hap)
    assert set(np.unique(keys)) == set(["0", "1"]), np.unique(keys)

    sd = StateDictionary.fit(keys, min_distinct_people=5, max_states=32)
    Z, cols = sd.transform(keys)
    g = phase.dosage(hap, 0)

    j1 = cols.index("1")
    j0 = cols.index("0")
    assert np.array_equal(Z[:, j1], g), (Z[:, j1], g)
    assert np.array_equal(Z[:, j0], 2.0 - g)

    rows = {}
    for i in range(Z.shape[0]):
        rows.setdefault(tuple(Z[i]), set()).add(float(g[i]))
    assert all(len(v) == 1 for v in rows.values())
    assert len(rows) == len(set(float(x) for x in g))

    sel = phase.select_pairs(["22:1000:A:G"], max_pairs=64)
    assert sel["n_possible_pairs"] == 0 and sel["n_selected"] == 0, sel
    assert_raises(phase.PhaseInputError, phase.phase_contrast_q, hap, 0, 0)
    _record({
        "states": sorted(sd.states),
        "state_count_determined_by_G": True,
        "n_pairs": 0,
        "new_local_phase_information": False,
    })

def test_04_complete_ld():
    rng = np.random.RandomState(20260922)
    bits = rng.randint(0, 2, size=(50, 2))
    hap = np.zeros((50, 2, 2), dtype=int)
    hap[:, :, 0] = bits
    hap[:, :, 1] = bits
    diag = phase.pair_diagnostics(hap, 0, 1, min_double_het=1, min_cis=1, min_trans=1)
    assert diag["n_double_het"] > 0, diag
    assert diag["n_trans"] == 0, diag
    assert abs(diag["haplotype_r2"] - 1.0) < 1e-12, diag
    assert diag["separable"] is False, diag
    assert "COMPLETE_LD_LIMITED_SEPARATION" in diag["flags"], diag
    assert "PHASE_CONTRAST_NOT_SEPARABLE" in diag["flags"], diag
    assert "TRANS_BELOW_THRESHOLD" in diag["flags"], diag
    assert diag["supported"] is False

    G = np.column_stack([phase.dosage(hap, 0), phase.dosage(hap, 1)])
    er = state.effective_rank(G, standardize=True)
    assert er["status"] == "OK", er
    assert abs(er["numerical_rank"] - 1.0) < 1e-9, er
    assert abs(er["effective_rank"] - 1.0) < 1e-9, er

    hap2 = synthetic_window(n_people=80, n_markers=4, seed=7)
    ok = phase.pair_diagnostics(hap2, 1, 3, min_double_het=1, min_cis=1, min_trans=1)
    assert ok["n_trans"] > 0 and ok["separable"] is True, ok
    _record({
        "complete_ld_r2": diag["haplotype_r2"],
        "n_cis": diag["n_cis"],
        "n_trans": diag["n_trans"],
        "flags": diag["flags"],
        "dosage_effective_rank": er["effective_rank"],
        "control_pair_separable": True,
    })

def test_05_alt_major():
    af = np.array([0.70, 0.20])
    flip = align.needs_flip(af)
    assert list(flip) == [True, False], flip
    assert list(align.minor_is_alt(af)) == [False, True]
    assert align.needs_flip(np.array([0.5]))[0] == False

    hap = make_hap([([1, 0], [1, 1]), ([0, 1], [1, 0]), ([1, 1], [0, 0])])
    oriented = align.orient_haplotypes(hap, af)
    assert np.array_equal(oriented[:, :, 0], 1 - hap[:, :, 0]), oriented
    assert np.array_equal(oriented[:, :, 1], hap[:, :, 1])

    ds_hard = align.dosage_from_haplotypes(hap)
    ds_oriented = align.dosage_from_haplotypes(oriented)
    assert np.allclose(ds_oriented[:, 0], 2.0 - ds_hard[:, 0])
    assert np.allclose(ds_oriented[:, 1], ds_hard[:, 1])
    assert np.allclose(align.minor_dosage(ds_hard, af), ds_oriented, equal_nan=True)

    assert np.allclose(align.minor_dosage(np.array([[1.4, 0.3]]), af), np.array([[0.6, 0.3]]))

    hap_miss = hap.copy()
    hap_miss[0, 1, 0] = MISSING
    om = align.orient_haplotypes(hap_miss, af)
    assert om[0, 1, 0] == MISSING, om[0]
    assert np.isnan(align.dosage_from_haplotypes(om)[0, 0])

    q_raw = phase.phase_contrast_q(hap, 0, 1)
    q_one = phase.phase_contrast_q(oriented, 0, 1)
    assert np.array_equal(q_one, -q_raw), (q_raw, q_one)
    both = align.orient_haplotypes(hap, np.array([0.70, 0.80]))
    q_both = phase.phase_contrast_q(both, 0, 1)
    assert np.array_equal(q_both, q_raw), (q_raw, q_both)

    k_raw = state.encode_subwindow(hap)
    k_or = state.encode_subwindow(oriented)
    assert k_raw[0, 0] == "10" and k_or[0, 0] == "00", (k_raw[0, 0], k_or[0, 0])

    assert_raises(OrientationSourceError, align.needs_flip, af, source="phenotype")
    _record({
        "flip_mask": [bool(x) for x in flip],
        "tie_af_0.5_minor": "ALT",
        "applied_to_both_copies": True,
        "missing_preserved": True,
        "q_sign_rule": "(-1)^(number of flipped markers)",
    })

def test_06_misalignment_fail_fast():
    ids = ["S{0:03d}".format(i) for i in range(10)]
    shuffled = list(ids)
    shuffled[3], shuffled[7] = shuffled[7], shuffled[3]
    assert align.assert_samples_aligned(ids, list(ids)) == 10
    assert set(shuffled) == set(ids)
    assert_raises(AlignmentError, align.assert_samples_aligned, ids, shuffled)
    assert_raises(AlignmentError, align.assert_samples_aligned, ids, ids[:-1])

    keys = [align.variant_key("chr22", 17000000 + i, "A", "G") for i in range(4)]
    assert keys[0] == "22:17000000:A:G", keys[0]
    bad_keys = list(keys)
    bad_keys[2] = align.variant_key("22", 17000002, "A", "T")
    assert align.assert_variant_keys_aligned(keys, list(keys)) == 4
    assert_raises(AlignmentError, align.assert_variant_keys_aligned, keys, bad_keys)
    assert_raises(AlignmentError, align.assert_variant_keys_aligned, keys, keys[::-1])

    hap = synthetic_window(n_people=10, n_markers=4, seed=1)
    assert align.check_inputs(ids, keys, hap=hap, af=np.full(4, 0.2))["n_markers"] == 4
    assert_raises(AlignmentError, align.check_inputs, ids[:9], keys, hap)
    assert_raises(AlignmentError, align.check_inputs, ids, keys[:3], hap)
    assert_raises(AlignmentError, align.check_inputs, ids, [keys[0]] * 4, hap)
    assert_raises(AlignmentError, align.check_inputs, ["S000"] * 10, keys, hap)
    assert_raises(AlignmentError, align.orient_haplotypes, hap, np.full(3, 0.2))
    _record({
        "order_permutation_blocked": True,
        "membership_equal_but_order_differs": "AlignmentError",
        "variant_key_mismatch_blocked": True,
        "shape_mismatch_blocked": True,
    })

def test_07_phase_set_crossing():
    meta_cross = marker_meta(3, ps=["17000000", "17000000", "18000000"])
    assert_raises(PhaseSetCrossingError, phase.check_phase_sets, meta_cross)
    ok_ps = phase.check_phase_sets(marker_meta(3, ps=["17000000"] * 3))
    assert ok_ps["provenance_grade"] == "PS_DECLARED", ok_ps
    assert_raises(
        PhaseSetCrossingError, phase.check_phase_sets, marker_meta(3, ps=["17000000", None, None])
    )
    m = marker_meta(2)
    m[1]["chrom"] = "21"
    assert_raises(PhaseSetCrossingError, phase.check_phase_sets, m)
    m2 = marker_meta(2)
    m2[1]["tile"] = "22:171"
    assert_raises(PhaseSetCrossingError, phase.check_phase_sets, m2)
    m3 = marker_meta(2)
    m3[1]["pos"] = m3[0]["pos"] + 25000
    assert_raises(PhaseSetCrossingError, phase.check_phase_sets, m3)
    crossed = phase.check_phase_sets(meta_cross, allow_cross=True)
    assert crossed["phase_set_source"] == "FORMAT/PS"

    ok = phase.check_phase_sets(marker_meta(16, ps=None))
    assert ok["provenance_grade"] == phase.ASSUMED_PHASE_SET_FLAG, ok
    assert ok["phase_set"] == "contig:22", ok
    assert ok["max_marker_gap_bp"] <= phase.DEFAULT_MAX_MARKER_GAP_BP
    unsorted_meta = marker_meta(3)
    unsorted_meta[0]["pos"], unsorted_meta[2]["pos"] = (
        unsorted_meta[2]["pos"],
        unsorted_meta[0]["pos"],
    )
    assert_raises(phase.PhaseInputError, phase.check_phase_sets, unsorted_meta)
    _record({
        "distinct_ps_blocked": True,
        "ps_and_none_mixed_blocked": True,
        "cross_chrom_blocked": True,
        "cross_tile_blocked": True,
        "gap_over_20kb_blocked": True,
        "ps_absent_grade": ok["provenance_grade"],
    })

def test_08_missing_not_reference():
    n_markers = 4
    allref = state.all_reference_key(n_markers)
    assert allref == "0000"

    rows = [([0] * 4, [0] * 4) for _ in range(30)]
    miss_copy = [0, MISSING, 0, 0]
    rows += [(miss_copy, [0] * 4) for _ in range(10)]
    rows += [(miss_copy, list(miss_copy)) for _ in range(10)]
    hap = make_hap(rows)
    keys = state.encode_subwindow(hap)
    assert keys[30, 0] == "0.00", keys[30, 0]
    assert state.has_missing(keys[30, 0]) and not state.has_missing(keys[0, 0])
    assert keys[30, 0] != allref

    sd = StateDictionary.fit(keys, min_distinct_people=5, max_states=32)
    assert sd.states == [allref], sd.states
    assert "0.00" not in sd.states
    assert sd.support[allref] == 40, sd.support
    assert sd.copy_counts[allref] == 70, sd.copy_counts
    assert sd.copy_counts["0.00"] == 30, sd.copy_counts
    assert sd.diagnostics["n_missing_bearing_copies"] == 30
    assert sd.diagnostics["n_people_with_any_missing_copy"] == 20
    assert sd.diagnostics["n_states_with_missing_markers"] == 1

    Z, cols = sd.transform(keys)
    jo = cols.index(OTHER)
    ja = cols.index(allref)
    assert Z[30, jo] == 1.0 and Z[30, ja] == 1.0, Z[30]
    assert Z[40, jo] == 2.0 and Z[40, ja] == 0.0, Z[40]
    assert Z[0, ja] == 2.0 and Z[0, jo] == 0.0
    assert np.all(Z.sum(axis=1) == 2.0)

    filled = hap.copy()
    filled[filled == MISSING] = 0
    sd_filled = StateDictionary.fit(state.encode_subwindow(filled), min_distinct_people=5)
    assert sd_filled.support[allref] == 50, sd_filled.support
    assert sd_filled.copy_counts[allref] == 100, sd_filled.copy_counts
    assert sd.support[allref] < sd_filled.support[allref], "결측이 reference 로 세어졌다"
    assert sd.copy_counts[allref] < sd_filled.copy_counts[allref]
    inv = state.state_inventory(sd, keys)
    assert inv["other_copy_share"] > 0.0

    hap2 = make_hap([([1, 1], [0, 0]), ([MISSING, 1], [0, 0])])
    q = phase.phase_contrast_q(hap2, 0, 1)
    assert q[0] == 1.0 and np.isnan(q[1]), q
    g = phase.dosage(hap2, 0)
    assert np.isnan(g[1]) and g[1] != 0
    sup = phase.pair_support(hap2, 0, 1, min_double_het=1, min_cis=1, min_trans=1)
    assert sup["n_missing"] == 1 and sup["n_double_het"] == 1, sup
    _record({
        "missing_encoding": state.MISSING_CHAR,
        "missing_key_is_not_all_reference": True,
        "missing_routed_to": OTHER,
        "allref_support_observed_vs_filled": [sd.support[allref], sd_filled.support[allref]],
        "Q_missing": "NA",
    })

def test_09_unseen_state_to_other():
    hap_a = synthetic_window(n_people=100, n_markers=4, seed=20260922)
    keys_a = state.encode_subwindow(hap_a)
    sd = StateDictionary.fit(keys_a, min_distinct_people=10, max_states=3)
    states_before = list(sd.states)
    assert len(states_before) == 3, states_before

    unseen1, unseen2 = "1010", "0101"
    assert unseen1 not in states_before and unseen2 not in states_before
    keys_b = np.array(
        [[states_before[0], unseen1], [unseen1, unseen2], [states_before[1], states_before[0]]],
        dtype=object,
    )
    Z, cols = sd.transform(keys_b)
    jo = cols.index(OTHER)
    assert Z[0, jo] == 1.0 and Z[1, jo] == 2.0 and Z[2, jo] == 0.0, Z
    assert np.all(Z.sum(axis=1) == 2.0)
    assert sd.states == states_before
    assert unseen1 not in sd.columns
    assert sd.map_key(unseen2) == OTHER

    dropped = sorted(set(str(k) for k in keys_a.ravel()) - set(states_before))
    assert dropped, "합성 입력에 탈락 state 가 없으면 시험 의미가 없다"
    assert all(sd.map_key(k) == OTHER for k in dropped)

    sd_no_other = StateDictionary.fit(
        keys_a, min_distinct_people=10, max_states=3, retain_other=False
    )
    assert_raises(state.StateEncodingError, sd_no_other.map_key, unseen1)
    _record({
        "n_states_A": len(states_before),
        "unseen_mapped_to": OTHER,
        "dictionary_expanded": False,
        "n_dropped_states_also_other": len(dropped),
    })

def test_10_label_leakage_sentinel():
    hap = synthetic_window(n_people=40, n_markers=4, seed=3)
    keys = state.encode_subwindow(hap)
    pheno_table = {"sample_id": ["S1"], "age": [50], "TCHL": [201.0]}
    covariate_table = {"sample_id": ["S1"], "age": [50], "sex_male": [1], "PC1": [0.01]}

    assert_raises(LabelLeakageError, state.encode_subwindow, hap, None, pheno_table)
    assert_raises(LabelLeakageError, StateDictionary.fit, keys, 5, 32, True, pheno_table)
    assert_raises(LabelLeakageError, state.assert_no_phenotype, {"y": [1, 2]})
    assert_raises(LabelLeakageError, state.assert_no_phenotype, ["snp1", "pheno_resid"])
    assert_raises(LabelLeakageError, state.assert_no_phenotype, {"y_res": [0.0]})
    assert state.assert_no_phenotype(covariate_table) is True
    sd = StateDictionary.fit(keys, min_distinct_people=5, table=covariate_table)

    assert_raises(LabelLeakageError, StateDictionary.fit, keys, 5, 32, True, None, np.zeros(40))
    assert_raises(LabelLeakageError, sd.transform, keys, None, np.zeros(40))
    assert_raises(
        LabelLeakageError, phase.pair_support, hap, 0, 1, 1, 1, 1, np.zeros(40)
    )
    assert_raises(
        LabelLeakageError, phase.select_pairs, ["22:1:A:G", "22:2:A:G"], 64, np.zeros(40)
    )
    assert_raises(LabelLeakageError, phase.select_pairs, ["22:1:A:G"], 64, None, pheno_table)
    assert_raises(LabelLeakageError, align.check_inputs, ["S1"], ["22:1:A:G"], None, None, None, pheno_table)
    for bad in ("phenotype", "Y", "TCHL", "B", "C"):
        assert_raises(OrientationSourceError, align.assert_orientation_source, bad)
    assert align.assert_orientation_source("INFO/AF") == "INFO/AF"

    Z1, cols1 = sd.transform(keys)
    Z2, cols2 = sd.transform(keys)
    assert cols1 == cols2 and np.array_equal(Z1, Z2)
    _record({
        "phenotype_columns_blocked": True,
        "phenotype_argument_blocked": True,
        "covariate_only_table_allowed": True,
        "orientation_from_phenotype_blocked": True,
    })

def test_11_constant_vector_na():
    const = np.full(20, 3.0)
    x = np.arange(20, dtype=float)
    r = state.safe_correlation(const, x)
    assert np.isnan(r), r
    assert r != 1.0
    assert np.isnan(state.safe_correlation(const, const))
    assert abs(state.safe_correlation(x, 2.0 * x + 1.0) - 1.0) < 1e-12

    M = np.column_stack([x, const])
    er = state.effective_rank(M, standardize=True)
    assert np.isnan(er["effective_rank"]), er
    assert er["status"] == "NA_CONSTANT_INPUT", er
    assert er["n_constant_columns"] == 1, er
    er_all = state.effective_rank(np.column_stack([const, const]), standardize=True)
    assert np.isnan(er_all["effective_rank"]) and er_all["status"] == "NA_CONSTANT_INPUT"
    assert state.effective_rank(np.zeros((0, 0)))["status"] == "NA_EMPTY"
    er_ok = state.effective_rank(np.column_stack([x, x[::-1] + 0.5 * x]), standardize=True)
    assert er_ok["status"] == "OK" and er_ok["effective_rank"] > 1.0, er_ok

    hap = np.zeros((20, 2, 2), dtype=int)
    hap[:, :, 0] = np.random.RandomState(0).randint(0, 2, size=(20, 2))
    diag = phase.pair_diagnostics(hap, 0, 1, min_double_het=1, min_cis=1, min_trans=1)
    assert np.isnan(diag["haplotype_r2"]), diag
    assert "COMPLETE_LD_LIMITED_SEPARATION" not in diag["flags"], diag
    _record({
        "constant_correlation": "NA",
        "standardized_effective_rank": "NA_CONSTANT_INPUT",
        "constant_not_reported_as_r1": True,
    })

def test_12_state_label_permutation():
    hap = synthetic_window(n_people=90, n_markers=4, seed=11)
    keys = state.encode_subwindow(hap)
    sd = StateDictionary.fit(keys, min_distinct_people=5, max_states=32)
    assert len(sd.states) >= 3, sd.states

    perm = list(range(len(sd.states)))[::-1]
    sd_perm = sd.reorder(perm)
    assert sd_perm.states != sd.states
    assert sorted(sd_perm.states) == sorted(sd.states)
    assert sd_perm.reference_state == sd.reference_state

    Z, cols = sd.transform(keys)
    Zp, colsp = sd_perm.transform(keys)
    assert cols != colsp
    by_key = sd.counts_by_key(keys)
    by_key_perm = sd_perm.counts_by_key(keys)
    assert set(by_key.keys()) == set(by_key_perm.keys())
    for k in by_key:
        assert np.array_equal(by_key[k], by_key_perm[k]), k
    order = [colsp.index(c) for c in cols]
    assert np.array_equal(Z, Zp[:, order])

    coef = dict((c, 0.37 * (i + 1)) for i, c in enumerate(cols))
    pred = state.predict_linear(Z, cols, coef, intercept=1.5)
    pred_perm = state.predict_linear(Zp, colsp, coef, intercept=1.5)
    assert np.allclose(pred, pred_perm), (pred[:3], pred_perm[:3])

    X, xc = sd.design_matrix(keys)
    Xp, xcp = sd_perm.design_matrix(keys)
    assert sorted(xc) == sorted(xcp) and sd.reference_state not in xc
    order2 = [xcp.index(c) for c in xc]
    assert np.array_equal(X, Xp[:, order2])
    _record({
        "n_states": len(sd.states),
        "counts_by_key_identical": True,
        "prediction_identical": True,
        "reference_state_stable": True,
    })

def _code_hashes():
    try:
        from hcf import workflow

        return workflow.code_sha256()
    except Exception:
        return []

def main(log_path=None, echo=True):
    if log_path is None:
        log_path = os.path.join(_HERE, "unit_tests.log")
    started = time.time()
    lines = []
    lines.append("# HCF-20260922-v1 unit_tests.log")
    lines.append("protocol: HCF-20260922-v1   track: skeleton_and_unit_tests")
    lines.append("scope: 합성 소형 입력만 (실자료·표현형·모형 적합 없음)")
    lines.append("utc_start: {0}".format(time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(started))))
    lines.append(
        "python: {0}   numpy: {1}   platform: {2}".format(
            sys.version.split()[0], np.__version__, platform.platform()
        )
    )
    lines.append("code_sha256:")
    for row in _code_hashes():
        lines.append("  {0}  {1}".format(row["sha256"], row["file"]))
    with open(os.path.abspath(__file__), "rb") as fh:
        import hashlib

        lines.append(
            "  {0}  {1}".format(
                hashlib.sha256(fh.read()).hexdigest(), os.path.basename(__file__)
            )
        )
    lines.append("")
    results = []
    passed = 0
    failed = 0
    module = sys.modules[__name__]
    for name in TEST_ORDER:
        fn = getattr(module, name)
        t0 = time.time()
        try:
            fn()
            detail = _DETAILS.get(name)
            status = "PASS"
            passed += 1
            err = ""
        except Exception as exc:
            status = "FAIL"
            failed += 1
            detail = None
            err = "{0}: {1}".format(type(exc).__name__, exc)
        dt = time.time() - t0
        results.append(
            {"test": name, "status": status, "seconds": round(dt, 4), "detail": detail, "error": err}
        )
        lines.append("[{0}] {1} ({2:.3f}s)".format(status, name, dt))
        lines.append("    검사 내용: {0}".format(DESCRIPTIONS[name]))
        if status == "PASS":
            lines.append("    관측: {0}".format(detail))
        else:
            lines.append("    오류: {0}".format(err))
            for tb_line in traceback.format_exc().rstrip().splitlines():
                lines.append("      " + tb_line)
        lines.append("")
    total = len(TEST_ORDER)
    lines.append("summary: total={0} passed={1} failed={2}".format(total, passed, failed))
    lines.append("verdict: {0}".format("ALL_PASS" if failed == 0 else "FAILURES_PRESENT"))
    lines.append("wall_seconds: {0:.3f}".format(time.time() - started))
    text = "\n".join(lines) + "\n"
    log_dir = os.path.dirname(os.path.abspath(log_path))
    if log_dir and not os.path.isdir(log_dir):
        os.makedirs(log_dir)
    with open(log_path, "w") as fh:
        fh.write(text)
    if echo:
        sys.stdout.write(text)
    return {
        "log_path": os.path.abspath(log_path),
        "total": total,
        "passed": passed,
        "failed": failed,
        "verdict": "ALL_PASS" if failed == 0 else "FAILURES_PRESENT",
        "results": results,
    }

if __name__ == "__main__":
    out = main()
    sys.exit(0 if out["failed"] == 0 else 1)
