import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import os

SRC = _config_path('${PROJECT_ROOT}/work/fset/att/att_fit.py')
DST = _config_path('${PROJECT_ROOT}/work/fset/att/att_fit_cvmin.py')

src = open(SRC).read()
orig_len = len(src)

anchor = 'import att_common as C\n'
assert src.count(anchor) == 1
inject = anchor + """
# --- A5-CVMIN 변형: 산출 경로 분리 (원본 fset/att/ 는 읽기만) ---
import os as _os
OUT = C.ATT + '/cvmin'
PRIVOUT = OUT + '/private'
for _p in (OUT, PRIVOUT):
    _os.makedirs(_p, exist_ok=True)
"""
src = src.replace(anchor, inject)

old_lam = 'LAM = [1e-3, 1e-2, 1e-1, 1.0]'
assert src.count(old_lam) == 1
src = src.replace(old_lam, 'LAM = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0]')

old_sel = """            best = {}
            for hid in (0, HID):
                bm = max(res[(hid, l)][0] for l in LAM)
                bl = [l for l in LAM if res[(hid, l)][0] == bm][0]
                bse = res[(hid, bl)][1]
                pick = max([l for l in LAM if res[(hid, l)][0] >= bm - bse])
                best[hid] = (res[(hid, pick)][0], pick, bse)
            lin_m, lin_l, _ = best[0]
            mlp_m, mlp_l, mlp_se = best[HID]
            hid_use, lam_use = ((0, lin_l) if lin_m >= mlp_m - mlp_se
                                else (HID, mlp_l))
            sel.append(dict(fold=f, hid=int(hid_use), lam=float(lam_use),
                            inner_linear=lin_m, inner_mlp=mlp_m,
                            mlp_se=mlp_se,
                            rule='1-SE on inner 5-fold; 선형이 MLP 의 1-SE 안이면 선형'))"""
new_sel = """            best = {}
            for hid in (0, HID):
                bm = max(res[(hid, l)][0] for l in LAM)
                bl = [l for l in LAM if res[(hid, l)][0] == bm][0]
                bse = res[(hid, bl)][1]
                best[hid] = (bm, bl, bse)      # CV-min: 1-SE 완화 없이 최고점
            lin_m, lin_l, _ = best[0]
            mlp_m, mlp_l, mlp_se = best[HID]
            hid_use, lam_use = (HID, mlp_l)    # MLP 게이트 강제
            sel.append(dict(fold=f, hid=int(hid_use), lam=float(lam_use),
                            inner_linear=lin_m, inner_mlp=mlp_m,
                            inner_linear_lam=lin_l, mlp_se=mlp_se,
                            rule='CV-min on inner 5-fold; MLP 게이트 강제 (A5-CVMIN 변형)'))"""
assert src.count(old_sel) == 1, 'selection block not matched'
src = src.replace(old_sel, new_sel)

writes = [
    ("C.ATT + '/att_label_shuffle_null.csv'", "OUT + '/att_label_shuffle_null.csv'"),
    ("C.ATT + '/att_gene_results.csv'", "OUT + '/att_gene_results.csv'"),
    ("C.ATT + '/att_weights.csv'", "OUT + '/att_weights.csv'"),
    ("C.ATT + '/att_fit_meta.json'", "OUT + '/att_fit_meta.json'"),
    ("C.ATT + '/att_fit.done'", "OUT + '/att_fit.done'"),
    ("C.PRIV + '/person_fold.npy'", "PRIVOUT + '/person_fold.npy'"),
    ("C.PRIV + '/null_delta_int_%02d.npy'", "PRIVOUT + '/null_delta_int_%02d.npy'"),
    ("C.PRIV + '/null_delta_%02d.npy'", "PRIVOUT + '/null_delta_%02d.npy'"),
    ("C.PRIV + '/null_z_%02d.npy'", "PRIVOUT + '/null_z_%02d.npy'"),
    ("C.PRIV + '/att_raw.npz'", "PRIVOUT + '/att_raw.npz'"),
]
for a, b in writes:
    assert src.count(a) == 1, 'write path not matched: ' + a
    src = src.replace(a, b)

old_meta = "post_hoc=True, amendment='A5')"
assert src.count(old_meta) == 1
src = src.replace(old_meta,
                  "post_hoc=True, amendment='A5-CVMIN',\n"
                  "                       variant=dict(\n"
                  "                           lam_rule='CV-min (no 1-SE relaxation)',\n"
                  "                           gate='MLP forced (hid=%d)' % HID,\n"
                  "                           base_run='fset/att (1-SE, 2026-09-23 15:02)'))")

assert "np.load(C.ATT + '/att_labels.npz'" in src
assert "C.BURD + '/col%02d.npy'" in src
assert "C.ATT + '/att_" not in src.split('att_labels.npz')[1], 'C.ATT 쓰기 잔존'

open(DST, 'w').write(src)
print('wrote %s  (%d -> %d bytes)' % (DST, orig_len, len(src)))
print('변경 요약: LAM 6점(1e-5~1.0) / CV-min / MLP 강제 / 산출 fset/att/cvmin/')
