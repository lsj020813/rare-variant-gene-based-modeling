
import hashlib
import math

import numpy as np

MISSING = -1

MISSING_CHAR = "."

OTHER = "OTHER"

PHENOTYPE_SENTINEL_NAMES = frozenset(
    [
        "y",
        "y_obs",
        "y_res",
        "y_syn",
        "tchl",
        "tchl_v3",
        "lip",
        "lip_v3",
        "pheno",
        "phenotype",
        "outcome",
        "label",
        "labels",
        "target",
        "trait",
        "response",
        "case_control",
    ]
)

PHENOTYPE_SENTINEL_PREFIXES = ("pheno", "tchl", "outcome_", "y_")

DEFAULT_MARKERS_PER_SUBWINDOW = 16
DEFAULT_MIN_DISTINCT_PEOPLE = 50
DEFAULT_MAX_SUPPORTED_STATES = 32

class LabelLeakageError(ValueError):
    pass

class StateEncodingError(ValueError):
    pass

def _candidate_names(obj):
    names = []
    cols = getattr(obj, "columns", None)
    if cols is not None:
        try:
            names.extend([str(c) for c in list(cols)])
        except TypeError:
            pass
    dtype = getattr(obj, "dtype", None)
    if dtype is not None and getattr(dtype, "names", None):
        names.extend([str(n) for n in dtype.names])
    if isinstance(obj, dict):
        names.extend([str(k) for k in obj.keys()])
    keys_attr = getattr(obj, "keys", None)
    if callable(keys_attr) and not isinstance(obj, dict):
        try:
            names.extend([str(k) for k in obj.keys()])
        except TypeError:
            pass
    if isinstance(obj, (list, tuple, set, frozenset)):
        names.extend([str(x) for x in obj if isinstance(x, str)])
    if isinstance(obj, str):
        names.append(obj)
    if getattr(obj, "name", None) is not None and isinstance(obj.name, str):
        names.append(obj.name)
    return names

def _is_phenotype_name(name):
    low = str(name).strip().lower()
    if low in PHENOTYPE_SENTINEL_NAMES:
        return True
    for pref in PHENOTYPE_SENTINEL_PREFIXES:
        if low.startswith(pref):
            return True
    return False

def assert_no_phenotype(obj, context="state"):
    hits = sorted(set(n for n in _candidate_names(obj) if _is_phenotype_name(n)))
    if hits:
        raise LabelLeakageError(
            "phenotype/label column(s) {0} reached an outcome-blind state function "
            "({1}); state construction must not see Y".format(hits, context)
        )
    return True

def key_hash(key):
    return hashlib.sha256(str(key).encode("utf-8")).hexdigest()

def encode_haplotype(bits):
    arr = np.asarray(bits, dtype=int).ravel()
    out = []
    for v in arr:
        if v == MISSING:
            out.append(MISSING_CHAR)
        elif v == 0:
            out.append("0")
        elif v == 1:
            out.append("1")
        else:
            raise StateEncodingError(
                "haplotype allele must be 0, 1 or MISSING(-1); got {0}".format(v)
            )
    return "".join(out)

def encode_subwindow(hap, sample_ids=None, table=None):
    if table is not None:
        assert_no_phenotype(table, context="state.encode_subwindow(table=...)")
    if sample_ids is not None:
        assert_no_phenotype(sample_ids, context="state.encode_subwindow(sample_ids=...)")
    arr = np.asarray(hap, dtype=int)
    if arr.ndim != 3 or arr.shape[1] != 2:
        raise StateEncodingError(
            "expected haplotype array of shape (n_people, 2, n_markers); got {0}".format(arr.shape)
        )
    n_people, _, n_markers = arr.shape
    if n_markers < 1:
        raise StateEncodingError("subwindow must contain at least 1 marker")
    keys = np.empty((n_people, 2), dtype=object)
    for i in range(n_people):
        for c in range(2):
            keys[i, c] = encode_haplotype(arr[i, c, :])
    return keys

def all_reference_key(n_markers):
    return "0" * int(n_markers)

def has_missing(key):
    return MISSING_CHAR in str(key)

class StateDictionary(object):

    def __init__(
        self,
        states,
        reference_state,
        n_markers,
        other_label=OTHER,
        support=None,
        copy_counts=None,
        diagnostics=None,
    ):
        self.states = list(states)
        if len(set(self.states)) != len(self.states):
            raise StateEncodingError("duplicate state keys in dictionary")
        if reference_state is not None and reference_state not in self.states + [other_label]:
            raise StateEncodingError("reference_state must be a retained state or OTHER")
        self.reference_state = reference_state
        self.n_markers = int(n_markers)
        self.other_label = other_label
        self.support = dict(support or {})
        self.copy_counts = dict(copy_counts or {})
        self.diagnostics = dict(diagnostics or {})

    @classmethod
    def fit(
        cls,
        hap_keys,
        min_distinct_people=DEFAULT_MIN_DISTINCT_PEOPLE,
        max_states=DEFAULT_MAX_SUPPORTED_STATES,
        retain_other=True,
        table=None,
        phenotype=None,
    ):
        if phenotype is not None:
            raise LabelLeakageError(
                "StateDictionary.fit received a phenotype argument; dictionary fitting is "
                "outcome-blind (brief §11.1, §15.1-10)"
            )
        if table is not None:
            assert_no_phenotype(table, context="StateDictionary.fit(table=...)")
        assert_no_phenotype(hap_keys, context="StateDictionary.fit(hap_keys=...)")

        keys = np.asarray(hap_keys, dtype=object)
        if keys.ndim != 2 or keys.shape[1] != 2:
            raise StateEncodingError(
                "hap_keys must have shape (n_people, 2); got {0}".format(keys.shape)
            )
        lengths = set(len(str(k)) for k in keys.ravel())
        if len(lengths) != 1:
            raise StateEncodingError(
                "state keys have inconsistent marker counts: {0}".format(sorted(lengths))
            )
        n_markers = lengths.pop()

        distinct = {}
        copies = {}
        n_missing_copies = 0
        n_people_with_missing = 0
        for i in range(keys.shape[0]):
            person_keys = [str(keys[i, 0]), str(keys[i, 1])]
            person_has_missing = False
            for k in person_keys:
                copies[k] = copies.get(k, 0) + 1
                if has_missing(k):
                    n_missing_copies += 1
                    person_has_missing = True
            if person_has_missing:
                n_people_with_missing += 1
            for k in set(person_keys):
                distinct[k] = distinct.get(k, 0) + 1

        eligible = dict((k, v) for k, v in distinct.items() if not has_missing(k))
        supported = dict((k, v) for k, v in eligible.items() if v >= int(min_distinct_people))

        ordered = sorted(supported.keys(), key=lambda k: (-supported[k], key_hash(k)))
        retained = ordered[: int(max_states)]
        reference = retained[0] if retained else (OTHER if retain_other else None)

        diagnostics = {
            "n_people_A": int(keys.shape[0]),
            "n_copies_A": int(keys.shape[0] * 2),
            "n_distinct_observed_states": int(len(distinct)),
            "n_states_with_missing_markers": int(len(distinct) - len(eligible)),
            "n_supported_states": int(len(supported)),
            "n_retained_states": int(len(retained)),
            "n_states_dropped_unsupported": int(len(eligible) - len(supported)),
            "n_states_dropped_over_cap": int(max(0, len(supported) - len(retained))),
            "n_missing_bearing_copies": int(n_missing_copies),
            "n_people_with_any_missing_copy": int(n_people_with_missing),
            "min_distinct_people_threshold": int(min_distinct_people),
            "max_supported_states": int(max_states),
            "all_reference_key_retained": all_reference_key(n_markers) in retained,
            "tie_break": "sha256(key) ascending",
            "support_threshold_is_power_proof": False,
        }
        return cls(
            states=retained,
            reference_state=reference,
            n_markers=n_markers,
            other_label=OTHER if retain_other else None,
            support=supported,
            copy_counts=copies,
            diagnostics=diagnostics,
        )

    @property
    def columns(self):
        cols = list(self.states)
        if self.other_label is not None:
            cols.append(self.other_label)
        return cols

    def map_key(self, key):
        k = str(key)
        if k in self.states:
            return k
        if self.other_label is None:
            raise StateEncodingError(
                "state {0!r} is not in the A dictionary and OTHER is disabled".format(k)
            )
        return self.other_label

    def transform(self, hap_keys, table=None, phenotype=None):
        if phenotype is not None:
            raise LabelLeakageError(
                "StateDictionary.transform received a phenotype argument; encoding is outcome-blind"
            )
        if table is not None:
            assert_no_phenotype(table, context="StateDictionary.transform(table=...)")
        keys = np.asarray(hap_keys, dtype=object)
        if keys.ndim != 2 or keys.shape[1] != 2:
            raise StateEncodingError("hap_keys must have shape (n_people, 2)")
        for k in keys.ravel():
            if len(str(k)) != self.n_markers:
                raise StateEncodingError(
                    "state key length {0} != dictionary n_markers {1}".format(
                        len(str(k)), self.n_markers
                    )
                )
        cols = self.columns
        index = dict((c, j) for j, c in enumerate(cols))
        Z = np.zeros((keys.shape[0], len(cols)), dtype=float)
        for i in range(keys.shape[0]):
            for c in range(2):
                Z[i, index[self.map_key(keys[i, c])]] += 1.0
        if not np.all(Z.sum(axis=1) == 2.0):
            raise StateEncodingError("state count rows must sum to 2")
        return Z, cols

    def design_matrix(self, hap_keys, table=None):
        Z, cols = self.transform(hap_keys, table=table)
        if self.reference_state is None:
            return Z, cols
        drop = cols.index(self.reference_state)
        keep = [j for j in range(len(cols)) if j != drop]
        return Z[:, keep], [cols[j] for j in keep]

    def reorder(self, new_order):
        new_states = [self.states[i] for i in new_order]
        if sorted(new_states) != sorted(self.states):
            raise StateEncodingError("reorder must be a permutation of the retained states")
        return StateDictionary(
            states=new_states,
            reference_state=self.reference_state,
            n_markers=self.n_markers,
            other_label=self.other_label,
            support=self.support,
            copy_counts=self.copy_counts,
            diagnostics=self.diagnostics,
        )

    def counts_by_key(self, hap_keys):
        Z, cols = self.transform(hap_keys)
        return dict((c, Z[:, j].copy()) for j, c in enumerate(cols))

def predict_linear(Z, columns, coef_by_key, intercept=0.0):
    Z = np.asarray(Z, dtype=float)
    beta = np.array([float(coef_by_key.get(c, 0.0)) for c in columns], dtype=float)
    return float(intercept) + Z.dot(beta)

def state_inventory(dictionary, hap_keys):
    Z, cols = dictionary.transform(hap_keys)
    copy_total = float(Z.sum())
    per_state = dict((c, float(Z[:, j].sum())) for j, c in enumerate(cols))
    allref = all_reference_key(dictionary.n_markers)
    allref_copies = per_state.get(allref, 0.0)
    nonref_total = copy_total - allref_copies
    nonref = dict((c, v) for c, v in per_state.items() if c != allref)
    top_nonref_share = (max(nonref.values()) / nonref_total) if (nonref and nonref_total > 0) else float("nan")
    return {
        "n_people": int(Z.shape[0]),
        "n_copies": int(copy_total),
        "n_columns": len(cols),
        "copies_per_state": per_state,
        "all_reference_key_present": allref in cols,
        "all_reference_copy_share": (allref_copies / copy_total) if copy_total > 0 else float("nan"),
        "other_copy_share": (per_state.get(dictionary.other_label, 0.0) / copy_total)
        if copy_total > 0
        else float("nan"),
        "top_share_excluding_all_reference": top_nonref_share,
        "diversity_claim": "NOT_ASSESSED_FROM_ALL_REFERENCE_MASS_ALONE",
    }

def safe_correlation(x, y):
    a = np.asarray(x, dtype=float).ravel()
    b = np.asarray(y, dtype=float).ravel()
    if a.size != b.size:
        raise ValueError("safe_correlation: length mismatch")
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 2:
        return float("nan")
    a, b = a[ok], b[ok]
    sa, sb = a.std(ddof=1), b.std(ddof=1)
    if sa == 0.0 or sb == 0.0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])

def effective_rank(matrix, standardize=True):
    M = np.asarray(matrix, dtype=float)
    if M.ndim != 2 or M.size == 0:
        return {
            "effective_rank": float("nan"),
            "numerical_rank": float("nan"),
            "status": "NA_EMPTY",
            "n_constant_columns": 0,
            "singular_values": [],
        }
    sd = M.std(axis=0, ddof=1) if M.shape[0] > 1 else np.zeros(M.shape[1])
    n_const = int(np.sum(sd == 0.0))
    X = M - M.mean(axis=0, keepdims=True)
    if standardize:
        if n_const > 0:
            return {
                "effective_rank": float("nan"),
                "numerical_rank": float("nan"),
                "status": "NA_CONSTANT_INPUT",
                "n_constant_columns": n_const,
                "singular_values": [],
            }
        X = X / sd
    sv = np.linalg.svd(X, compute_uv=False)
    total = float(sv.sum())
    if total <= 0.0:
        return {
            "effective_rank": float("nan"),
            "numerical_rank": 0.0,
            "status": "NA_CONSTANT_INPUT",
            "n_constant_columns": n_const,
            "singular_values": [float(s) for s in sv],
        }
    p = sv / total
    p = p[p > 0]
    ent = -float(np.sum(p * np.log(p)))
    return {
        "effective_rank": math.exp(ent),
        "numerical_rank": float(np.linalg.matrix_rank(X)),
        "status": "OK",
        "n_constant_columns": n_const,
        "singular_values": [float(s) for s in sv],
    }

__all__ = [
    "MISSING",
    "MISSING_CHAR",
    "OTHER",
    "PHENOTYPE_SENTINEL_NAMES",
    "DEFAULT_MARKERS_PER_SUBWINDOW",
    "DEFAULT_MIN_DISTINCT_PEOPLE",
    "DEFAULT_MAX_SUPPORTED_STATES",
    "LabelLeakageError",
    "StateEncodingError",
    "assert_no_phenotype",
    "key_hash",
    "encode_haplotype",
    "encode_subwindow",
    "all_reference_key",
    "has_missing",
    "StateDictionary",
    "predict_linear",
    "state_inventory",
    "safe_correlation",
    "effective_rank",
]
