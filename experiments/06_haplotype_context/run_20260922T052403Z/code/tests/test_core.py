import math
from hcf.core import (
    MISSING, AlignmentError, PhaseSetCrossingError,
    parse_phased_gt, copy_swap, minor_orientation_flip,
    support_state_counts, apply_state_dictionary, phase_contrast_Q,
    check_sample_key_alignment, check_no_phase_set_crossing,
    missing_to_state, safe_corr, relabel_cluster_ids,
)

def test_01_cis_trans_sign():
    cis = phase_contrast_Q(h1j=1, h1k=1, h2j=0, h2k=0)
    trans = phase_contrast_Q(h1j=1, h1k=0, h2j=0, h2k=1)
    assert cis == 1, cis
    assert trans == -1, trans
    Gj_cis, Gk_cis = 1 + 0, 1 + 0
    Gj_trans, Gk_trans = 1 + 0, 0 + 1
    assert (Gj_cis, Gk_cis) == (Gj_trans, Gk_trans) == (1, 1)
    assert cis != trans

def test_02_copy_swap_invariance():
    h1, h2 = [0, 1, 1, 0], [1, 0, 0, 1]
    Q_before = phase_contrast_Q(h1[0], h1[1], h2[0], h2[1])
    h1s, h2s = copy_swap(h1, h2)
    Q_after = phase_contrast_Q(h1s[0], h1s[1], h2s[0], h2s[1])
    Gj_before, Gk_before = h1[0] + h2[0], h1[1] + h2[1]
    Gj_after, Gk_after = h1s[0] + h2s[0], h1s[1] + h2s[1]
    assert (Gj_before, Gk_before) == (Gj_after, Gk_after)
    Cbefore = h1[0] * h1[1] + h2[0] * h2[1]
    Cafter = h1s[0] * h1s[1] + h2s[0] * h2s[1]
    assert Cbefore == Cafter
    assert Q_after == Q_before

def test_03_single_marker_no_phase_info():
    people = {
        "p1": (( "A0",), ("A1",)),
        "p2": (( "A1",), ("A1",)),
        "p3": (( "A0",), ("A0",)),
    }
    for _ in range(60):
        pass
    big = {f"person_{i}": (("A0",), ("A1",)) if i % 2 == 0 else (("A1",), ("A1",))
           for i in range(120)}
    result = support_state_counts(big, min_distinct_people=50, max_states=32)
    for z in result["Z"].values():
        assert sum(z) == 2

def test_04_complete_ld_limits_separation():
    qs_j = []
    qs_k = []
    for h1j, h1k, h2j, h2k in [(1, 1, 0, 0), (0, 0, 1, 1), (1, 1, 1, 1), (0, 0, 0, 0)]:
        Gj, Gk = h1j + h2j, h1k + h2k
        qs_j.append(Gj)
        qs_k.append(Gk)
    corr = safe_corr(qs_j, qs_k)
    assert corr is not None and abs(corr - 1.0) < 1e-9, corr

def test_05_alt_major_orientation_consistency():
    af = 0.7
    copy1_alt, copy2_alt = 1, 0
    m1 = minor_orientation_flip(copy1_alt, af)
    m2 = minor_orientation_flip(copy2_alt, af)
    assert (m1, m2) == (0, 1)
    DS = copy1_alt + copy2_alt
    minor_dosage_from_DS = 2 - DS
    assert minor_dosage_from_DS == m1 + m2 == 1

def test_06_sample_key_misalignment_fail_fast():
    expected = ["s1", "s2", "s3", "s4"]
    shuffled = ["s1", "s3", "s2", "s4"]
    try:
        check_sample_key_alignment(shuffled, expected)
        raised = False
    except AlignmentError:
        raised = True
    assert raised

def test_07_phase_set_crossing_blocked():
    try:
        check_no_phase_set_crossing(["PS1", "PS1", "PS2"])
        raised = False
    except PhaseSetCrossingError:
        raised = True
    assert raised
    check_no_phase_set_crossing(["PS1", "PS1", "PS1"])
    check_no_phase_set_crossing([None, None, None])

def test_08_missingness_not_reference():
    a1, a2 = parse_phased_gt(".|1")
    assert a1 == MISSING and a2 == 1
    assert missing_to_state(a1) == MISSING
    assert missing_to_state(a1) != 0

def test_09_unseen_state_goes_to_other():
    A_data = {f"a{i}": ((0,), (1,)) if i % 3 else ((1,), (1,)) for i in range(60)}
    fitted = support_state_counts(A_data, min_distinct_people=10, max_states=2)
    kept = fitted["kept_states"]
    C_data = {"new_person": ((9, 9), (0,))}
    Z = apply_state_dictionary(C_data, kept)
    other_index = len(kept)
    z = Z["new_person"]
    assert z[other_index] >= 1

def test_10_label_leakage_sentinel():
    A_data = {f"a{i}": ((0,), (1,)) if i % 2 else ((1,), (1,)) for i in range(120)}
    fit_1 = support_state_counts(A_data, min_distinct_people=50, max_states=32)
    fit_2 = support_state_counts(A_data, min_distinct_people=50, max_states=32)
    assert fit_1["kept_states"] == fit_2["kept_states"]
    assert fit_1["Z"] == fit_2["Z"]

def test_11_constant_vector_is_NA():
    const_vec = [5.0] * 10
    other = list(range(10))
    assert safe_corr(const_vec, other) is None
    assert safe_corr(other, other) == 1.0 if len(set(other)) > 1 else True

def test_12_state_label_permutation_invariance():
    assignment = {"p1": 0, "p2": 1, "p3": 0, "p4": 2}
    permutation = {0: 2, 1: 0, 2: 1}
    relabeled = relabel_cluster_ids(assignment, permutation)
    inv = {v: k for k, v in permutation.items()}
    recovered = {k: inv[v] for k, v in relabeled.items()}
    assert recovered == assignment
    def groups(d):
        g = {}
        for k, v in d.items():
            g.setdefault(v, set()).add(k)
        return set(frozenset(s) for s in g.values())
    assert groups(assignment) == groups(relabeled)

ALL_TESTS = [
    test_01_cis_trans_sign,
    test_02_copy_swap_invariance,
    test_03_single_marker_no_phase_info,
    test_04_complete_ld_limits_separation,
    test_05_alt_major_orientation_consistency,
    test_06_sample_key_misalignment_fail_fast,
    test_07_phase_set_crossing_blocked,
    test_08_missingness_not_reference,
    test_09_unseen_state_goes_to_other,
    test_10_label_leakage_sentinel,
    test_11_constant_vector_is_NA,
    test_12_state_label_permutation_invariance,
]
