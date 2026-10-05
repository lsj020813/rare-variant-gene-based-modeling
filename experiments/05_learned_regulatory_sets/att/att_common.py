import os as _cfg_os
import re as _cfg_re
import math as _cfg_math

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)

def _config_number(name, cast=int, positive=False):
    value = cast(_cfg_os.environ[name])
    if not _cfg_math.isfinite(value) or (positive and value <= 0):
        raise ValueError(name + " must be finite" + (" and positive" if positive else ""))
    return value

import os as _os
N_SAMPLES = _config_number("N_SAMPLES", int, True)
import os
import numpy as np

FSET = _config_path('${PROJECT_ROOT}/work/fset')
ATT = FSET + '/att'
PRIV = ATT + '/private'
BURD = PRIV + '/burden'
LOGS = ATT + '/logs'
MASTER = FSET + '/f3f6/cache/master.npz'
VCFD = _config_path('${PROJECT_ROOT}/work/ref/orig_index')
PHENO = _config_path('${PROJECT_ROOT}/work/ref/pheno_v3/tchl_v3.tsv')

SEED = 20260923
NSHUF = 20
NTOT = N_SAMPLES
NDOM = 1000
NCOL = 32
MAF_MAIN = 0.5
MAF_SEC = 0.01
COVARS = ['age', 'sex_male', 'CT', 'NC', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5']

C_TOT_ALL = 0
C_TOT_LAB = 1
C_CCRE = 2
C_TOT001_ALL = 3
C_TOT001_LAB = 4
C_CCRE001 = 5
C_MAJ001 = 6
C_MAJ0 = 7
C_K_ALL = 28
C_K_LAB = 29
C_K_MAJ = 30
C_K_CCRE = 31
CARRIER_DS = 0.5

def master():
    return np.load(MASTER, allow_pickle=True)

def ensure_dirs():
    for p in (ATT, PRIV, BURD, LOGS):
        os.makedirs(p, exist_ok=True)
