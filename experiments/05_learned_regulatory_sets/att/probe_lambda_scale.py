import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import sys
import time
import numpy as np
import torch

sys.path.insert(0, _config_path('${PROJECT_ROOT}/work/fset/att'))
import att_common as C
import att_fit_cvmin as F

torch.set_num_threads(4)

d = C.master()
Z = np.load(C.ATT + '/att_labels.npz', allow_pickle=True)
RAW, CM, gok = Z['raw'], Z['cmean'], Z['dom_ok']

fold_all = F.person_folds()
idx, y, Xc = F.load_pheno()
fold = fold_all[idx]

totlab = F.load_col(C.C_TOT_LAB, idx)
maj = F.load_col(C.C_MAJ0, idx)
mino = totlab - maj
F0, gm0 = F.feature_matrix(RAW[0], CM[0])
gmask = gm0 & gok
print('도메인 %d개, 개인 %d명' % (gmask.sum(), len(y)), flush=True)

Xmj = torch.tensor(F0[:, 0, :], dtype=torch.float32)
Xmn = torch.tensor(F0[:, 1, :], dtype=torch.float32)

f = 0
tr = np.where(fold != f)[0]
yt = F.resid_y(y, Xc, tr)
rng = np.random.default_rng(C.SEED + 1000 + f)
ip = rng.permutation(len(tr))
blocks = [tr[ip[k::F.NINNER]] for k in range(F.NINNER)]
o = F.block_stats(maj, mino, yt, blocks)
m = F.moments(o, list(range(F.NINNER)))
sa, sb = np.sqrt(m['va']), np.sqrt(m['vb'])
oa, ob = sa > 1e-9, sb > 1e-9
st = F.to_st(m, sa, sb, oa, ob)
base = gmask & oa & ob
msk = torch.tensor(base.astype(np.float32))
print('유효 도메인 %d' % int(base.sum()), flush=True)

print('\n%-10s %8s %8s %8s %8s %8s %8s' %
      ('lam', 'w_평균', 'w_표준편차', 'w_최소', 'w_최대', '|w-.5|>.05', 'evidence'))
for lam in (1.0, 1e-2, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 0.0):
    t0 = time.time()
    net = F.fit_gate(Xmj, Xmn, st, msk, lam, F.HID, C.SEED + f)
    with torch.no_grad():
        w = F.gate_w(net, Xmj, Xmn).numpy()
        pen = float(sum((p * p).sum() for p in net.parameters()))
    ev = F.eval_gate(net, Xmj, Xmn, st, msk)
    wb = w[base]
    print('%-10s %8.4f %8.4f %8.4f %8.4f %8d %10.3e   (‖p‖²=%.3f, %.0fs)'
          % (lam, wb.mean(), wb.std(), wb.min(), wb.max(),
             int((np.abs(wb - 0.5) > 0.05).sum()), ev, pen,
             time.time() - t0), flush=True)
