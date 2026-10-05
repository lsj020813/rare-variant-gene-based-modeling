
import numpy as np

from hcf.state import MISSING, LabelLeakageError, assert_no_phenotype, safe_correlation

DEFAULT_MIN_DOUBLE_HET_PEOPLE = 100
DEFAULT_MIN_CIS_PEOPLE = 25
DEFAULT_MIN_TRANS_PEOPLE = 25
DEFAULT_MAX_PAIRS_PER_TILE = 64
DEFAULT_MAX_MARKER_GAP_BP = 20000

ASSUMED_PHASE_SET_FLAG = "VCF_ASSUMED_CONTIG_PHASE"

class PhaseSetCrossingError(ValueError):
    pass

class PhaseInputError(ValueError):
    pass

def _as_copy_pair(hap):
    arr = np.asarray(hap, dtype=int)
    if arr.ndim != 3 or arr.shape[1] != 2:
        raise PhaseInputError(
            "expected haplotype array of shape (n_people, 2, n_markers); got {0}".format(arr.shape)
        )
    return arr

def dosage(hap, j):
    arr = _as_copy_pair(hap)
    h1 = arr[:, 0, j].astype(float)
    h2 = arr[:, 1, j].astype(float)
    out = h1 + h2
    out[(arr[:, 0, j] == MISSING) | (arr[:, 1, j] == MISSING)] = np.nan
    return out

def co_occurrence(hap, j, k):
    arr = _as_copy_pair(hap)
    sel = arr[:, :, [j, k]]
    miss = np.any(sel == MISSING, axis=(1, 2))
    c = (arr[:, 0, j] * arr[:, 0, k] + arr[:, 1, j] * arr[:, 1, k]).astype(float)
    c[miss] = np.nan
    return c

def phase_contrast_q(hap, j, k, check_identity=True):
    arr = _as_copy_pair(hap)
    if j == k:
        raise PhaseInputError("phase contrast requires two distinct markers (j != k)")
    d1 = (arr[:, 0, j] - arr[:, 1, j]).astype(float)
    d2 = (arr[:, 0, k] - arr[:, 1, k]).astype(float)
    q = d1 * d2
    sel = arr[:, :, [j, k]]
    miss = np.any(sel == MISSING, axis=(1, 2))
    q[miss] = np.nan
    if check_identity:
        alt = 2.0 * co_occurrence(arr, j, k) - dosage(arr, j) * dosage(arr, k)
        both = ~miss
        if both.any() and not np.allclose(q[both], alt[both]):
            raise PhaseInputError("Q identity 2C - Gj*Gk != (h1j-h2j)(h1k-h2k) violated")
    return q

def phase_class(hap, j, k):
    q = phase_contrast_q(hap, j, k)
    out = np.empty(q.shape[0], dtype=object)
    for i, v in enumerate(q):
        if not np.isfinite(v):
            out[i] = "missing"
        elif v > 0:
            out[i] = "cis"
        elif v < 0:
            out[i] = "trans"
        else:
            out[i] = "homozygous_at_one_locus"
    return out

def is_double_het(hap, j, k):
    arr = _as_copy_pair(hap)
    gj = dosage(arr, j)
    gk = dosage(arr, k)
    return (gj == 1.0) & (gk == 1.0)

def pair_support(
    hap,
    j,
    k,
    min_double_het=DEFAULT_MIN_DOUBLE_HET_PEOPLE,
    min_cis=DEFAULT_MIN_CIS_PEOPLE,
    min_trans=DEFAULT_MIN_TRANS_PEOPLE,
    phenotype=None,
):
    if phenotype is not None:
        raise LabelLeakageError(
            "pair_support received a phenotype argument; pair selection is outcome-blind "
            "(resolved_config.phase_contrasts.pair_selection_uses_phenotype = false)"
        )
    q = phase_contrast_q(hap, j, k)
    finite = np.isfinite(q)
    n_cis = int(np.sum(finite & (q > 0)))
    n_trans = int(np.sum(finite & (q < 0)))
    n_dh = n_cis + n_trans
    reasons = []
    if n_dh < int(min_double_het):
        reasons.append("DOUBLE_HET_BELOW_THRESHOLD")
    if n_cis < int(min_cis):
        reasons.append("CIS_BELOW_THRESHOLD")
    if n_trans < int(min_trans):
        reasons.append("TRANS_BELOW_THRESHOLD")
    return {
        "marker_j": int(j),
        "marker_k": int(k),
        "n_people": int(q.shape[0]),
        "n_missing": int(np.sum(~finite)),
        "n_double_het": n_dh,
        "n_cis": n_cis,
        "n_trans": n_trans,
        "supported": len(reasons) == 0,
        "unsupported_reasons": reasons,
        "support_threshold_is_power_proof": False,
    }

def pair_diagnostics(hap, j, k, **kwargs):
    arr = _as_copy_pair(hap)
    sup = pair_support(arr, j, k, **kwargs)
    copies_j = np.concatenate([arr[:, 0, j], arr[:, 1, j]]).astype(float)
    copies_k = np.concatenate([arr[:, 0, k], arr[:, 1, k]]).astype(float)
    copies_j[copies_j == MISSING] = np.nan
    copies_k[copies_k == MISSING] = np.nan
    r = safe_correlation(copies_j, copies_k)
    r2 = float("nan") if not np.isfinite(r) else float(r * r)
    q = phase_contrast_q(arr, j, k)
    finite_nonzero = np.isfinite(q) & (q != 0)
    distinct_q = sorted(set(float(v) for v in q[finite_nonzero]))
    flags = list(sup["unsupported_reasons"])
    if np.isfinite(r2) and r2 >= 1.0 - 1e-12:
        flags.append("COMPLETE_LD_LIMITED_SEPARATION")
    if len(distinct_q) < 2:
        flags.append("PHASE_CONTRAST_NOT_SEPARABLE")
    out = dict(sup)
    out.update(
        {
            "haplotype_r2": r2,
            "haplotype_r": r,
            "distinct_q_values_among_double_het": distinct_q,
            "flags": sorted(set(flags)),
            "separable": len(distinct_q) >= 2,
        }
    )
    return out

def select_pairs(marker_keys, max_pairs=DEFAULT_MAX_PAIRS_PER_TILE, phenotype=None, table=None):
    if phenotype is not None:
        raise LabelLeakageError("select_pairs must not see the phenotype")
    if table is not None:
        assert_no_phenotype(table, context="phase.select_pairs(table=...)")
    from hcf.state import key_hash

    keys = [str(k) for k in marker_keys]
    pairs = []
    for a in range(len(keys)):
        for b in range(a + 1, len(keys)):
            pairs.append((a, b))
    pairs.sort(key=lambda ab: key_hash("{0}|{1}".format(keys[ab[0]], keys[ab[1]])))
    selected = pairs[: int(max_pairs)]
    inclusion_probability = (
        float(len(selected)) / float(len(pairs)) if pairs else float("nan")
    )
    return {
        "pairs": selected,
        "n_possible_pairs": len(pairs),
        "n_selected": len(selected),
        "inclusion_probability": inclusion_probability,
        "selection_rule": "sha256 hash order of joined marker keys; phenotype-blind",
    }

def check_phase_sets(
    marker_meta,
    allow_cross=False,
    max_marker_gap_bp=DEFAULT_MAX_MARKER_GAP_BP,
):
    if not marker_meta:
        raise PhaseInputError("marker_meta is empty")
    chroms = set(str(m["chrom"]) for m in marker_meta)
    if len(chroms) > 1:
        raise PhaseSetCrossingError(
            "subwindow spans multiple chromosomes {0}; phase cannot be joined".format(
                sorted(chroms)
            )
        )
    tiles = set(str(m.get("tile")) for m in marker_meta if m.get("tile") is not None)
    if len(tiles) > 1:
        raise PhaseSetCrossingError(
            "subwindow spans multiple tiles {0}; windows must not cross tile boundaries".format(
                sorted(tiles)
            )
        )
    ps_values = [m.get("ps") for m in marker_meta]
    present = set(str(v) for v in ps_values if v is not None)
    if len(present) > 1 and not allow_cross:
        raise PhaseSetCrossingError(
            "subwindow joins distinct phase sets {0}; resolved_config.phase.cross_phase_set "
            "is false".format(sorted(present))
        )
    if present and len(present) == 1 and any(v is None for v in ps_values):
        raise PhaseSetCrossingError(
            "subwindow mixes markers with and without PS; phase range is undetermined"
        )
    positions = [int(m["pos"]) for m in marker_meta]
    if any(positions[i] > positions[i + 1] for i in range(len(positions) - 1)):
        raise PhaseInputError("marker_meta must be sorted by position")
    gaps = [positions[i + 1] - positions[i] for i in range(len(positions) - 1)]
    max_gap = max(gaps) if gaps else 0
    if max_gap > int(max_marker_gap_bp):
        raise PhaseSetCrossingError(
            "marker gap {0} bp exceeds max_marker_gap_bp {1}; subwindow must not span it".format(
                max_gap, max_marker_gap_bp
            )
        )
    if present:
        return {
            "phase_set": sorted(present)[0],
            "phase_set_source": "FORMAT/PS",
            "provenance_grade": "PS_DECLARED",
            "chrom": sorted(chroms)[0],
            "n_markers": len(marker_meta),
            "max_marker_gap_bp": int(max_gap),
        }
    return {
        "phase_set": "contig:{0}".format(sorted(chroms)[0]),
        "phase_set_source": "ABSENT_PS_HEADER_PHASING_FULL",
        "provenance_grade": ASSUMED_PHASE_SET_FLAG,
        "chrom": sorted(chroms)[0],
        "n_markers": len(marker_meta),
        "max_marker_gap_bp": int(max_gap),
        "note": "PS 부재 -> 염색체 단위 위상 가정. switch-error rate 는 알려지지 않았다.",
    }

def swap_copies(hap):
    arr = _as_copy_pair(hap)
    return arr[:, ::-1, :].copy()

__all__ = [
    "DEFAULT_MIN_DOUBLE_HET_PEOPLE",
    "DEFAULT_MIN_CIS_PEOPLE",
    "DEFAULT_MIN_TRANS_PEOPLE",
    "DEFAULT_MAX_PAIRS_PER_TILE",
    "DEFAULT_MAX_MARKER_GAP_BP",
    "ASSUMED_PHASE_SET_FLAG",
    "PhaseSetCrossingError",
    "PhaseInputError",
    "dosage",
    "co_occurrence",
    "phase_contrast_q",
    "phase_class",
    "is_double_het",
    "pair_support",
    "pair_diagnostics",
    "select_pairs",
    "check_phase_sets",
    "swap_copies",
]
