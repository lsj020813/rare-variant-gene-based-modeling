import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


SRC = _config_path('${PROJECT_ROOT}/work/fset/att/att_burden.py')
DST = _config_path('${PROJECT_ROOT}/work/fset/att/att_burden_std.py')
s = open(SRC).read()

s = s.replace('import att_common as C', 'import att_common_std as C', 1)

old = "    Z = np.load(C.ATT + '/att_labels.npz', allow_pickle=True)"
assert s.count(old) == 1
s = s.replace(old, "    Z = np.load(C.LABELS, allow_pickle=True)")

old = """    LAB = Z['labels']              # (21, 364430) int8
    ok = Z['dom_ok']"""
new = """    LAB = Z['labels']              # (NPART, 364430) int8
    ok = Z['dom_ok']
    assert LAB.shape[0] == C.NPART"""
assert s.count(old) == 1
s = s.replace(old, new)

old = """    doms = [g for g in range(C.NDOM)
            if dom_chrom[g].replace('chr', '') == str(ch)]"""
new = """    doms = [g for g in range(C.NDOM)
            if dom_chrom[g].replace('chr', '') == str(ch)
            and g % C.DOMSTEP == 0]"""
assert s.count(old) == 1
s = s.replace(old, new)

old = """        # 열 가중 행렬 W (M x NCOL): 0/1 지시자
        W = np.zeros((M, C.NCOL), dtype=np.float32)
        mm = maf[rid]
        cc = hcc[rid] > 0
        lo = mm < C.MAF_SEC
        lab0 = LAB[0][rid]
        W[:, C.C_TOT_ALL] = 1.0
        W[:, C.C_TOT_LAB] = (lab0 >= 0)
        W[:, C.C_CCRE] = cc
        W[:, C.C_TOT001_ALL] = lo
        W[:, C.C_TOT001_LAB] = lo & (lab0 >= 0)
        W[:, C.C_CCRE001] = lo & cc
        W[:, C.C_MAJ001] = lo & (lab0 == 0)
        for L in range(C.NSHUF + 1):
            W[:, C.C_MAJ0 + L] = (LAB[L][rid] == 0)
        W[~got, :] = 0.0
        out[ii, :C.C_K_ALL] = (D @ W[:, :C.C_K_ALL]).T
        Kw = np.zeros((M, 4), dtype=np.float32)
        Kw[:, 0] = 1.0
        Kw[:, 1] = (lab0 >= 0)
        Kw[:, 2] = (lab0 == 0)
        Kw[:, 3] = cc
        Kw[~got, :] = 0.0
        out[ii, C.C_K_ALL:] = ((D >= C.CARRIER_DS).astype(np.float32)
                               @ Kw).T"""
new = """        # 변이별 표준편차 (관측 dosage 기준). 표준화 가중 = 1/σ
        sd = D.std(0)
        usable = got & (sd > 1e-6)
        inv = np.zeros(M, dtype=np.float32)
        inv[usable] = 1.0 / sd[usable]
        mm = maf[rid]
        cc = hcc[rid] > 0
        W = np.zeros((M, C.NCOL), dtype=np.float32)
        W[:, C.C_RAW_ALL] = got.astype(np.float32)
        W[:, C.C_STD_ALL] = inv
        W[:, C.C_STD_COMMON] = inv * (mm >= C.MAF_COMMON)
        W[:, C.C_STD_LOW] = inv * ((mm >= C.MAF_LOW) & (mm < C.MAF_COMMON))
        W[:, C.C_STD_RARE] = inv * ((mm >= C.MAF_RARE) & (mm < C.MAF_LOW))
        W[:, C.C_STD_VRARE] = inv * (mm < C.MAF_RARE)
        W[:, C.C_RAW_COMMON] = got * (mm >= C.MAF_COMMON)
        W[:, C.C_STD_CCRE] = inv * cc
        W[:, C.C_STD_NOCCRE] = inv * (~cc)
        for L in range(C.NPART):
            lp = LAB[L][rid]
            W[:, C.C_P0 + 2 * L] = inv * (lp == 0)
            W[:, C.C_P0 + 2 * L + 1] = inv * (lp == 1)
        W[~got, :] = 0.0
        out[ii] = (D @ W).T"""
assert s.count(old) == 1, 'W 블록 미일치'
s = s.replace(old, new)

open(DST, 'w').write(s)
import ast
ast.parse(s)
print('wrote %s, 문법 OK' % DST)
