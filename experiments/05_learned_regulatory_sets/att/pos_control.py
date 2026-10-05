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
import csv
import json
import time
import numpy as np
import torch
from scipy.stats import norm

sys.path.insert(0, _config_path('${PROJECT_ROOT}/work/fset/att'))
import att_common as C
import att_fit_lamfix as F

torch.set_num_threads(4)

NREP = int(sys.argv[1]) if len(sys.argv) > 1 else 2
OUTD = C.ATT + '/' + (sys.argv[2] if len(sys.argv) > 2 else 'poscontrol')
import os
DELTAS = [float(x) for x in
          os.environ.get('POS_DELTAS', '0.0,0.005,0.01,0.02,0.05').split(',')]
NCAUSAL = int(os.environ.get('POS_NCAUSAL', '50'))
import os
os.makedirs(OUTD, exist_ok=True)

def zsc(A, fold):
    return F.arm_single(A, np.zeros(A.shape[1], dtype=np.float32), fold)

def main():
    t00 = time.time()
    d = C.master()
    Z = np.load(C.ATT + '/att_labels.npz', allow_pickle=True)
    RAW, CM, gok = Z['raw'], Z['cmean'], Z['dom_ok']

    fold_all = F.person_folds()
    idx = np.arange(C.NTOT)
    fold = fold_all[idx]
    N = len(idx)

    totlab = F.load_col(C.C_TOT_LAB, idx)
    maj = F.load_col(C.C_MAJ0, idx)
    mino = totlab - maj
    tot = F.load_col(C.C_TOT_ALL, idx)
    F0, gm0 = F.feature_matrix(RAW[0], CM[0])
    gmask = gm0 & gok
    Xc = np.zeros((N, 0), dtype=np.float32)
    print('[setup] N=%d 도메인=%d 유효=%d %.0fs'
          % (N, C.NDOM, gmask.sum(), time.time() - t00), flush=True)

    Zmaj = zsc(maj, fold)
    Zmin = zsc(mino, fold)

    frac_ccre_maj = RAW[0][:, 0, 5]
    frac_ccre_min = RAW[0][:, 1, 5]
    truth = {
        'P1_MAJOR': np.ones(C.NDOM, dtype=bool),
        'P2_CCRE': frac_ccre_maj >= frac_ccre_min,
        'P3_RANDOM': np.random.default_rng(C.SEED + 31337
                                           ).random(C.NDOM) < 0.5,
    }
    only = os.environ.get('POS_SCENARIOS')
    if only:
        truth = {k: v for k, v in truth.items() if k in only.split(',')}

    rows = []
    for scen, is_maj in truth.items():
        for rep in range(NREP):
            rng = np.random.default_rng(C.SEED + 7000 + rep)
            cand = np.where(gmask)[0]
            causal = rng.choice(cand, size=NCAUSAL, replace=False)
            for delta in DELTAS:
                ts = time.time()
                sig = np.zeros(N, dtype=np.float32)
                for g in causal:
                    sig += (Zmaj[g] if is_maj[g] else Zmin[g])
                sig = np.nan_to_num(sig)
                sig = sig / (sig.std() + 1e-12)
                ysyn = (delta * sig
                        + rng.standard_normal(N).astype(np.float32))
                ysyn = (ysyn - ysyn.mean()) / ysyn.std()

                S_att, w, sel = F.run_att(maj, mino, F0, gmask, ysyn, Xc,
                                          fold)
                _, z_att = F.zstats(S_att, ysyn, Xc)
                _, z_flat = F.zstats(F.arm_single(tot, ysyn, fold), ysyn, Xc)
                q_att = F.bh(2.0 * norm.sf(np.abs(z_att)))
                q_flat = F.bh(2.0 * norm.sf(np.abs(z_flat)))
                cm = np.zeros(C.NDOM, dtype=bool)
                cm[causal] = True
                wc = w[cm]
                want_maj = is_maj[cm]
                toward = np.where(want_maj, wc > 0.5, wc < 0.5)
                rows.append(dict(
                    scenario=scen, rep=rep, delta=delta,
                    n_causal=int(NCAUSAL),
                    w_mean=float(np.nanmean(w[gmask])),
                    w_sd=float(np.nanstd(w[gmask])),
                    w_causal_mean=float(np.nanmean(wc)),
                    frac_w_toward_truth=float(np.nanmean(toward)),
                    n_sig_att=int(np.nansum(q_att < 0.05)),
                    n_sig_flat=int(np.nansum(q_flat < 0.05)),
                    n_sig_att_causal=int(np.nansum(q_att[cm] < 0.05)),
                    n_sig_flat_causal=int(np.nansum(q_flat[cm] < 0.05)),
                    mean_delta_z=float(np.nanmean(np.abs(z_att)
                                                  - np.abs(z_flat))),
                    mean_delta_z_causal=float(np.nanmean(
                        np.abs(z_att[cm]) - np.abs(z_flat[cm]))),
                    lam0=float(sel[0]['lam']), hid0=int(sel[0]['hid']),
                    wall_s=round(time.time() - ts, 1)))
                print('[%s rep%d d=%.4f] w_causal=%.3f toward=%.2f '
                      'nsig ATT=%d FLAT=%d (인과 %d/%d) dz=%.4f lam=%g %.0fs'
                      % (scen, rep, delta, rows[-1]['w_causal_mean'],
                         rows[-1]['frac_w_toward_truth'],
                         rows[-1]['n_sig_att'], rows[-1]['n_sig_flat'],
                         rows[-1]['n_sig_att_causal'],
                         rows[-1]['n_sig_flat_causal'],
                         rows[-1]['mean_delta_z_causal'],
                         rows[-1]['lam0'], rows[-1]['wall_s']), flush=True)

    with open(OUTD + '/pos_control.csv', 'w') as fh:
        wtr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        wtr.writeheader()
        for r in rows:
            wtr.writerow(r)
    with open(OUTD + '/pos_control_meta.json', 'w') as fh:
        json.dump(dict(deltas=DELTAS, nrep=NREP, n_causal=NCAUSAL,
                       lam_grid=F.LAM, hid=F.HID, epochs=F.EPOCHS,
                       nfold=F.NFOLD, ninner=F.NINNER, seed=C.SEED,
                       n_individuals=int(N), covariates='none (synthetic)',
                       scenarios=list(truth.keys()),
                       wall_s=round(time.time() - t00, 1)), fh, indent=1)
    print('DONE %.1f s' % (time.time() - t00))

if __name__ == '__main__':
    main()
