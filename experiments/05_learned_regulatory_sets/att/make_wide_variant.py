import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


SRC = _config_path('${PROJECT_ROOT}/work/fset/att/att_fit_lamfix.py')
DST = _config_path('${PROJECT_ROOT}/work/fset/att/att_fit_wide.py')
s = open(SRC).read()
n0 = len(s)

old = "OUT = C.ATT + '/lamfix'"
assert s.count(old) == 1
s = s.replace(old, "OUT = C.ATT + '/wide'")

old = "HID = 16"
assert s.count(old) == 1
s = s.replace(old, "HID = 16\n"
                   "# 게이트 출력 범위. (0,1) 이면 볼록 결합만 가능해 대조 방향이\n"
                   "# 배제된다. ceiling_w.py 실측상 [-0.5,1.5] 에서 천장이 포화.\n"
                   "WLO, WHI = -0.5, 1.5")

old = """def gate_w(net, Xmaj, Xmin):
    return torch.sigmoid(net(Xmaj).squeeze(-1) - net(Xmin).squeeze(-1))"""
new = """def gate_w(net, Xmaj, Xmin):
    s = torch.sigmoid(net(Xmaj).squeeze(-1) - net(Xmin).squeeze(-1))
    return WLO + (WHI - WLO) * s"""
assert s.count(old) == 1
s = s.replace(old, new)

old = """    out['ATT'] = zstats(S_att, y, Xc)"""
new = """    out['ATT'] = zstats(S_att, y, Xc)
    # 게이트 없는 고정 팔: 단독 모듈과 대조(major-minor)
    _one = np.ones(C.NDOM)
    _zero = np.zeros(C.NDOM)
    S_maj = arm_fixed_w(maj, mino, _one, y, fold)
    S_min = arm_fixed_w(maj, mino, _zero, y, fold)
    out['MAJONLY'] = zstats(S_maj, y, Xc)
    out['MINONLY'] = zstats(S_min, y, Xc)
    out['CONTRAST'] = zstats(S_maj - S_min, y, Xc)
    del S_maj, S_min"""
assert s.count(old) == 1
s = s.replace(old, new)

old = ("    arms = ['FLAT', 'FIXED', 'SIZE', 'ATT', 'INT', 'INT_product_only',\n"
       "            'ATT_GENEHOLD', 'ATT_MAF001', 'FLAT_MAF001', 'FIXED_MAF001']")
new = ("    arms = ['FLAT', 'FIXED', 'SIZE', 'ATT', 'INT', 'INT_product_only',\n"
       "            'ATT_GENEHOLD', 'ATT_MAF001', 'FLAT_MAF001', 'FIXED_MAF001',\n"
       "            'MAJONLY', 'MINONLY', 'CONTRAST']")
assert s.count(old) == 1, 'arms list not matched'
s = s.replace(old, new)

old = """        SIL, SPL = arm_int(tot, mj, mn, y, Xc, fold)"""
new = """        _o = np.ones(C.NDOM)
        _z = np.zeros(C.NDOM)
        _sm = arm_fixed_w(mj, mn, _o, y, fold)
        _sn = arm_fixed_w(mj, mn, _z, y, fold)
        _, zCL = zstats(_sm - _sn, y, Xc)
        qCL = bh(2.0 * norm.sf(np.abs(zCL)))
        dCL = np.abs(zCL) - np.abs(zflat)
        del _sm, _sn
        SIL, SPL = arm_int(tot, mj, mn, y, Xc, fold)"""
assert s.count(old) == 1
s = s.replace(old, new)

old = """                         n_sig_int_fdr05=int(np.nansum(qIL < 0.05)),"""
new = """                         n_sig_int_fdr05=int(np.nansum(qIL < 0.05)),
                         n_sig_contrast_fdr05=int(np.nansum(qCL < 0.05)),
                         mean_delta_contrast=float(np.nanmean(dCL)),
                         max_abs_z_contrast=float(np.nanmax(np.abs(zCL))),"""
assert s.count(old) == 1
s = s.replace(old, new)

s = s.replace("amendment='A5-LAMFIX'", "amendment='A5-WIDE'")
s = s.replace("gate='MLP forced (hid=%d)' % HID",
              "gate='MLP forced (hid=%d), w range (%.2f,%.2f)' % (HID, WLO, WHI)")

open(DST, 'w').write(s)
import ast
ast.parse(s)
print('wrote %s (%d -> %d bytes), 문법 OK' % (DST, n0, len(s)))
