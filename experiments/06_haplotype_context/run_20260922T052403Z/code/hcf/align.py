
import numpy as np

from hcf.state import MISSING, LabelLeakageError, assert_no_phenotype

ALLOWED_ORIENTATION_SOURCES = ("INFO/AF",)

FORBIDDEN_ORIENTATION_SOURCES = (
    "phenotype",
    "Y",
    "y",
    "TCHL",
    "B",
    "C",
    "development_split",
    "evaluation_split",
)

class AlignmentError(ValueError):
    pass

class OrientationSourceError(ValueError):
    pass

def variant_key(chrom, pos, ref, alt):
    chrom_s = str(chrom)
    if chrom_s.lower().startswith("chr"):
        chrom_s = chrom_s[3:]
    return "{0}:{1}:{2}:{3}".format(chrom_s, int(pos), str(ref).upper(), str(alt).upper())

def assert_orientation_source(source):
    if source in FORBIDDEN_ORIENTATION_SOURCES:
        raise OrientationSourceError(
            "allele orientation must come from upstream INFO/AF only, never from "
            "phenotype or B/C; got source={0!r}".format(source)
        )
    if source not in ALLOWED_ORIENTATION_SOURCES:
        raise OrientationSourceError(
            "unknown orientation source {0!r}; allowed={1}".format(
                source, list(ALLOWED_ORIENTATION_SOURCES)
            )
        )
    return source

def needs_flip(af, source="INFO/AF"):
    assert_orientation_source(source)
    af_arr = np.asarray(af, dtype=float)
    if np.any(~np.isfinite(af_arr)):
        raise AlignmentError("INFO/AF contains non-finite values; cannot decide orientation")
    if np.any(af_arr < 0.0) or np.any(af_arr > 1.0):
        raise AlignmentError("INFO/AF outside [0,1]")
    out = af_arr > 0.5
    return bool(out) if out.ndim == 0 else out

def minor_is_alt(af, source="INFO/AF"):
    return np.logical_not(needs_flip(af, source=source))

def minor_dosage(ds, af, source="INFO/AF"):
    ds_arr = np.asarray(ds, dtype=float)
    flip = np.asarray(needs_flip(af, source=source))
    if flip.ndim == 0:
        flip = np.full(ds_arr.shape[-1] if ds_arr.ndim else 1, bool(flip))
    if ds_arr.ndim and ds_arr.shape[-1] != flip.shape[-1]:
        raise AlignmentError(
            "DS marker axis ({0}) does not match AF length ({1})".format(
                ds_arr.shape[-1], flip.shape[-1]
            )
        )
    out = np.where(flip, 2.0 - ds_arr, ds_arr)
    return out

def orient_haplotypes(hap, af, source="INFO/AF"):
    hap_arr = np.asarray(hap)
    if hap_arr.ndim < 2 or hap_arr.shape[-2] != 2:
        raise AlignmentError(
            "haplotype array must have shape (..., 2, n_markers); got {0}".format(hap_arr.shape)
        )
    flip = np.atleast_1d(np.asarray(needs_flip(af, source=source)))
    if flip.shape[-1] != hap_arr.shape[-1]:
        raise AlignmentError(
            "haplotype marker axis ({0}) does not match AF length ({1})".format(
                hap_arr.shape[-1], flip.shape[-1]
            )
        )
    out = np.array(hap_arr, dtype=int, copy=True)
    miss = out == MISSING
    flipped = np.where(flip, 1 - out, out)
    flipped[miss] = MISSING
    return flipped

def dosage_from_haplotypes(hap):
    hap_arr = np.asarray(hap, dtype=int)
    if hap_arr.ndim < 2 or hap_arr.shape[-2] != 2:
        raise AlignmentError("haplotype array must have shape (..., 2, n_markers)")
    miss = np.any(hap_arr == MISSING, axis=-2)
    total = hap_arr.sum(axis=-2).astype(float)
    total[miss] = np.nan
    return total

def assert_samples_aligned(ids_left, ids_right, name_left="left", name_right="right"):
    left = [str(x) for x in ids_left]
    right = [str(x) for x in ids_right]
    if len(left) != len(right):
        raise AlignmentError(
            "sample count mismatch: {0}={1} vs {2}={3}".format(
                name_left, len(left), name_right, len(right)
            )
        )
    bad = [i for i in range(len(left)) if left[i] != right[i]]
    if bad:
        same_set = set(left) == set(right)
        raise AlignmentError(
            "sample order mismatch between {0} and {1}: {2} position(s) differ, "
            "first at index {3}; same_membership={4} (reordering is NOT performed "
            "silently)".format(name_left, name_right, len(bad), bad[0], same_set)
        )
    return len(left)

def assert_variant_keys_aligned(keys_left, keys_right, name_left="left", name_right="right"):
    left = [str(x) for x in keys_left]
    right = [str(x) for x in keys_right]
    if len(left) != len(right):
        raise AlignmentError(
            "variant count mismatch: {0}={1} vs {2}={3}".format(
                name_left, len(left), name_right, len(right)
            )
        )
    bad = [i for i in range(len(left)) if left[i] != right[i]]
    if bad:
        raise AlignmentError(
            "variant key order mismatch between {0} and {1}: {2} position(s) differ, "
            "first at index {3} ({4!r} vs {5!r})".format(
                name_left, name_right, len(bad), bad[0], left[bad[0]], right[bad[0]]
            )
        )
    return len(left)

def assert_matrix_aligned(matrix, sample_ids, variant_keys, name="matrix"):
    arr = np.asarray(matrix)
    if arr.shape[0] != len(sample_ids):
        raise AlignmentError(
            "{0} rows ({1}) != sample id count ({2})".format(name, arr.shape[0], len(sample_ids))
        )
    if arr.shape[-1] != len(variant_keys):
        raise AlignmentError(
            "{0} marker axis ({1}) != variant key count ({2})".format(
                name, arr.shape[-1], len(variant_keys)
            )
        )
    return arr.shape

def check_inputs(sample_ids, variant_keys, hap=None, dosage=None, af=None, table=None):
    if table is not None:
        assert_no_phenotype(table, context="align.check_inputs(table=...)")
    if len(set(str(k) for k in variant_keys)) != len(list(variant_keys)):
        raise AlignmentError("duplicate variant keys in subwindow")
    if len(set(str(s) for s in sample_ids)) != len(list(sample_ids)):
        raise AlignmentError("duplicate sample ids")
    if hap is not None:
        assert_matrix_aligned(hap, sample_ids, variant_keys, name="haplotypes")
    if dosage is not None:
        assert_matrix_aligned(dosage, sample_ids, variant_keys, name="dosage")
    if af is not None and len(np.atleast_1d(np.asarray(af))) != len(list(variant_keys)):
        raise AlignmentError("INFO/AF length != variant key count")
    return {
        "n_samples": len(list(sample_ids)),
        "n_markers": len(list(variant_keys)),
        "orientation_source": ALLOWED_ORIENTATION_SOURCES[0],
    }

__all__ = [
    "ALLOWED_ORIENTATION_SOURCES",
    "FORBIDDEN_ORIENTATION_SOURCES",
    "AlignmentError",
    "OrientationSourceError",
    "LabelLeakageError",
    "variant_key",
    "assert_orientation_source",
    "needs_flip",
    "minor_is_alt",
    "minor_dosage",
    "orient_haplotypes",
    "dosage_from_haplotypes",
    "assert_samples_aligned",
    "assert_variant_keys_aligned",
    "assert_matrix_aligned",
    "check_inputs",
]
