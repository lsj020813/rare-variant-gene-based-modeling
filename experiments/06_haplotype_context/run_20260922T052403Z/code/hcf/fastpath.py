
import numpy as np

from hcf import state as hstate

_CODES = np.array([46, 48, 49], dtype=np.uint8)

def encode_keys(hap):
    arr = np.asarray(hap)
    if arr.ndim != 3 or arr.shape[1] != 2:
        raise ValueError("expected (n_people, 2, m); got %s" % (arr.shape,))
    m = arr.shape[2]
    idx = arr.astype(np.int16) + 1
    if idx.size and (idx.min() < 0 or idx.max() > 2):
        raise ValueError("haplotype allele must be 0, 1 or MISSING(-1)")
    b = _CODES[idx]
    s = np.ascontiguousarray(b).reshape(-1, m).view("S%d" % m).reshape(arr.shape[0], 2)
    return np.char.decode(s, "ascii")

def assert_equivalent_encoding(hap, keys, n_check=256, rng=None):
    n = hap.shape[0]
    if n == 0:
        return 0
    rng = rng or np.random.RandomState(20260922)
    k = min(n_check, n)
    sel = rng.choice(n, size=k, replace=False)
    ref = hstate.encode_subwindow(np.asarray(hap)[sel])
    got = keys[sel]
    for i in range(k):
        for c in range(2):
            if str(ref[i, c]) != str(got[i, c]):
                raise AssertionError(
                    "fastpath encoding differs from hcf.state.encode_subwindow at row %d copy %d"
                    % (int(sel[i]), c)
                )
    return k

def transform_counts(dictionary, keys):
    cols = list(dictionary.columns)
    index = dict((c, j) for j, c in enumerate(cols))
    other = index.get(dictionary.other_label) if dictionary.other_label is not None else None
    n = keys.shape[0]
    Z = np.zeros((n, len(cols)), dtype=np.float64)
    rows = np.arange(n)
    for cpy in range(2):
        col = keys[:, cpy]
        uniq, inv = np.unique(col, return_inverse=True)
        lut = np.empty(uniq.size, dtype=np.int64)
        for u in range(uniq.size):
            key = str(uniq[u])
            if key in index:
                lut[u] = index[key]
            elif other is None:
                raise hstate.StateEncodingError(
                    "state %r not in dictionary and OTHER disabled" % key
                )
            else:
                lut[u] = other
        np.add.at(Z, (rows, lut[inv]), 1.0)
    if not np.all(Z.sum(axis=1) == 2.0):
        raise hstate.StateEncodingError("state count rows must sum to 2")
    return Z, cols

def assert_equivalent_transform(dictionary, keys, Z, n_check=512, rng=None):
    n = keys.shape[0]
    if n == 0:
        return 0
    rng = rng or np.random.RandomState(20260923)
    k = min(n_check, n)
    sel = rng.choice(n, size=k, replace=False)
    ref, cols = dictionary.transform(np.asarray(keys, dtype=object)[sel])
    if not np.allclose(ref, Z[sel]):
        raise AssertionError("fastpath transform differs from StateDictionary.transform")
    return k
