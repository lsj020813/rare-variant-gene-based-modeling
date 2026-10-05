from __future__ import annotations
import hashlib
import math
from dataclasses import dataclass
from typing import Sequence, Optional

MISSING = -1

class AlignmentError(ValueError):
    pass

class PhaseSetCrossingError(ValueError):
    pass

def parse_phased_gt(gt: str) -> tuple[int, int]:
    if "/" in gt and "|" not in gt:
        raise AlignmentError(f"unphased GT separator in phased-required context: {gt!r}")
    sep = "|"
    if sep not in gt:
        raise AlignmentError(f"GT has neither '|' nor '/': {gt!r}")
    a, b = gt.split(sep, 1)
    c1 = MISSING if a == "." else int(a)
    c2 = MISSING if b == "." else int(b)
    return c1, c2

def copy_swap(h1: Sequence[int], h2: Sequence[int]) -> tuple[Sequence[int], Sequence[int]]:
    return list(h2), list(h1)

def minor_orientation_flip(allele: int, alt_af: float, *, tie_to_alt: bool = True) -> int:
    if allele == MISSING:
        return MISSING
    if alt_af > 0.5:
        return 1 - allele
    if alt_af == 0.5 and not tie_to_alt:
        return 1 - allele
    return allele

@dataclass(frozen=True)
class MarkerCall:
    copy1: int
    copy2: int

def support_state_counts(
    haplotype_bits_by_person: dict[str, tuple[int, ...]],
    *,
    min_distinct_people: int = 50,
    max_states: int = 32,
) -> dict:
    distinct_people_per_pattern: dict[tuple, set] = {}
    for person, (c1, c2) in haplotype_bits_by_person.items():
        for pat in (c1, c2):
            distinct_people_per_pattern.setdefault(pat, set()).add(person)

    def sort_key(pat):
        n = len(distinct_people_per_pattern[pat])
        h = hashlib.sha256(repr(pat).encode()).hexdigest()
        return (-n, h)

    supported = [p for p, s in distinct_people_per_pattern.items() if len(s) >= min_distinct_people]
    supported.sort(key=sort_key)
    kept = supported[:max_states]
    kept_set = set(kept)

    state_of = {pat: i for i, pat in enumerate(kept)}
    OTHER = len(kept)

    Z = {}
    for person, (c1, c2) in haplotype_bits_by_person.items():
        z = [0] * (len(kept) + 1)
        for pat in (c1, c2):
            idx = state_of.get(pat, OTHER)
            z[idx] += 1
        Z[person] = tuple(z)
        if sum(z) != 2:
            raise AssertionError(f"state-count vector for {person} does not sum to 2: {z}")

    return {"kept_states": kept, "OTHER_index": OTHER, "Z": Z,
            "n_distinct_people_per_kept_state": [len(distinct_people_per_pattern[p]) for p in kept]}

def apply_state_dictionary(
    haplotype_bits_by_person: dict[str, tuple[tuple, tuple]],
    kept_states: list[tuple],
) -> dict[str, tuple]:
    state_of = {pat: i for i, pat in enumerate(kept_states)}
    OTHER = len(kept_states)
    Z = {}
    for person, (c1, c2) in haplotype_bits_by_person.items():
        z = [0] * (len(kept_states) + 1)
        for pat in (c1, c2):
            z[state_of.get(pat, OTHER)] += 1
        Z[person] = tuple(z)
    return Z

def phase_contrast_Q(h1j: int, h1k: int, h2j: int, h2k: int) -> Optional[int]:
    if MISSING in (h1j, h1k, h2j, h2k):
        return None
    Gj = h1j + h2j
    Gk = h1k + h2k
    C = h1j * h1k + h2j * h2k
    Q = 2 * C - Gj * Gk
    assert Q == (h1j - h2j) * (h1k - h2k)
    return Q

def check_sample_key_alignment(vcf_sample_order: Sequence[str], expected_order: Sequence[str]) -> None:
    if list(vcf_sample_order) != list(expected_order):
        raise AlignmentError(
            f"sample order mismatch: len(vcf)={len(vcf_sample_order)} "
            f"len(expected)={len(expected_order)}, first mismatch at "
            f"{next((i for i, (a, b) in enumerate(zip(vcf_sample_order, expected_order)) if a != b), None)}"
        )

def check_no_phase_set_crossing(ps_by_fragment: Sequence[Optional[str]]) -> None:
    seen = {ps for ps in ps_by_fragment if ps is not None}
    if len(seen) > 1:
        raise PhaseSetCrossingError(f"multiple phase-sets in one window: {seen}")

def missing_to_state(allele: int) -> int:
    return allele

def safe_corr(x: Sequence[float], y: Sequence[float]) -> Optional[float]:
    n = len(x)
    if n == 0 or n != len(y):
        raise AlignmentError("length mismatch or empty vectors")
    mx = sum(x) / n
    my = sum(y) / n
    vx = sum((xi - mx) ** 2 for xi in x)
    vy = sum((yi - my) ** 2 for yi in y)
    if vx == 0 or vy == 0:
        return None
    cov = sum((xi - mx) * (yi - my) for xi, yi in zip(x, y))
    return cov / math.sqrt(vx * vy)

def relabel_cluster_ids(assignment: dict[str, int], permutation: dict[int, int]) -> dict[str, int]:
    return {k: permutation[v] for k, v in assignment.items()}
