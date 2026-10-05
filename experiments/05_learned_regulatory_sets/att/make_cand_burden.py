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
DST = _config_path('${PROJECT_ROOT}/work/fset/att/att_burden_cand.py')
s = open(SRC).read()

old = 'import att_common as C'
assert s.count(old) == 1
s = s.replace(old, 'import att_common_cand as C')

old = "    Z = np.load(C.ATT + '/att_labels.npz', allow_pickle=True)"
new = "    Z = np.load(C.LABELS, allow_pickle=True)"
assert s.count(old) == 1
s = s.replace(old, new)

old = """    LAB = Z['labels']              # (21, 364430) int8
    ok = Z['dom_ok']"""
new = """    LAB = Z['labels']              # (NPART, 364430) int8
    ok = Z['dom_ok']
    assert LAB.shape[0] == C.NPART, 'NPART 불일치 %d' % LAB.shape[0]"""
assert s.count(old) == 1
s = s.replace(old, new)

old = """    doms = [g for g in range(C.NDOM)
            if dom_chrom[g].replace('chr', '') == str(ch)]"""
new = """    doms = [g for g in range(C.NDOM)
            if dom_chrom[g].replace('chr', '') == str(ch)
            and g % C.DOMSTEP == 0]"""
assert s.count(old) == 1
s = s.replace(old, new)

old = """        for L in range(C.NSHUF + 1):
            W[:, C.C_MAJ0 + L] = (LAB[L][rid] == 0)"""
new = """        for L in range(C.NPART):
            lp = LAB[L][rid]
            W[:, C.C_P0 + 2 * L] = (lp == 0)
            W[:, C.C_P0 + 2 * L + 1] = (lp == 1)"""
assert s.count(old) == 1
s = s.replace(old, new)

old = "        lab0 = LAB[0][rid]"
assert s.count(old) == 1

open(DST, 'w').write(s)
import ast
ast.parse(s)
print('wrote %s, 문법 OK' % DST)
