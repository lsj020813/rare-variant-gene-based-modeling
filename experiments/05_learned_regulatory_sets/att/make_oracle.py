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
DST = _config_path('${PROJECT_ROOT}/work/fset/att/oracle_domain.py')
s = open(SRC).read()

s = s.replace('import att_common as C', 'import att_common_oracle as C', 1)

old = "    Z = np.load(C.ATT + '/att_labels.npz', allow_pickle=True)"
assert s.count(old) == 1
s = s.replace(old, "    Z = np.load(C.LABELS, allow_pickle=True)")

old = """    LAB = Z['labels']              # (21, 364430) int8
    ok = Z['dom_ok']"""
assert s.count(old) == 1
s = s.replace(old, """    ok = Z['dom_ok']
    PH = C.pheno()""")

old = """    doms = [g for g in range(C.NDOM)
            if dom_chrom[g].replace('chr', '') == str(ch)]"""
assert s.count(old) == 1
s = s.replace(old, """    doms = [g for g in range(C.NDOM)
            if dom_chrom[g].replace('chr', '') == str(ch)
            and g % C.DOMSTEP == 0]""")

old = """    from numpy.lib.format import open_memmap
    out = open_memmap(C.BURD + '/burden_chr%s.npy' % ch, mode='w+',
                      dtype=np.float32, shape=(len(doms), C.NCOL, C.NTOT))
    rep = []"""
assert s.count(old) == 1
s = s.replace(old, "    rep = []")

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
new = """        mm = maf[rid]
        rep.append(C.oracle_domain(D, got, mm, int(g), str(ch), PH))"""
assert s.count(old) == 1
s = s.replace(old, new)

old = """        del D
        rep.append(dict(domain_idx=int(g), chrom=str(ch), n_var=int(M),
                        n_matched=int(got.sum()), match_rate=mrate,
                        n_records_parsed=int(recs), n_intervals=len(iv),
                        dom_ok=bool(ok[g]),
                        n_ccre=int((cc & got).sum()),
                        n_maf_lt001=int((lo & got).sum()),
                        status='ok' if got.sum() >= 2 else 'too_few_matched'))
    out.flush()
    del out
    np.save(C.BURD + '/domidx_chr%s.npy' % ch,
            np.array(doms, dtype=np.int32))"""
new = """        rep[-1].update(dict(n_matched=int(got.sum()), match_rate=mrate,
                            dom_ok=bool(ok[g])))
        del D"""
assert s.count(old) == 1
s = s.replace(old, new)

old = """    with open(C.ATT + '/burden_chr%s.json' % ch, 'w') as f:
        json.dump(dict(chrom=str(ch), n_domains=len(doms),
                       wall_s=round(time.time() - t0, 1), per_domain=rep), f)"""
new = """    with open(C.ATT + '/oracle_chr%s.json' % ch, 'w') as f:
        json.dump(dict(chrom=str(ch), n_domains=len(doms),
                       wall_s=round(time.time() - t0, 1), per_domain=rep), f)"""
assert s.count(old) == 1
s = s.replace(old, new)

s = s.replace("with open(C.ATT + '/burden_chr%s.done' % ch, 'w') as f:",
              "with open(C.ATT + '/oracle_chr%s.done' % ch, 'w') as f:")

open(DST, 'w').write(s)
import ast
ast.parse(s)
print('wrote %s, 문법 OK' % DST)
