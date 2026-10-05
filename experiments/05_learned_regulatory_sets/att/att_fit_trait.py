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
import sys

sys.path.insert(0, _config_path('${PROJECT_ROOT}/work/fset/att'))
import att_common as C

TRAIT = os.environ['TRAIT']
C.PHENO = _config_path('${PROJECT_ROOT}/work/ref/pheno_v3/%s_v3.tsv') % TRAIT
assert os.path.exists(C.PHENO), C.PHENO

import att_fit_lamfix as F

F.OUT = C.ATT + '/trait_' + TRAIT
F.PRIVOUT = F.OUT + '/private'
for p in (F.OUT, F.PRIVOUT):
    os.makedirs(p, exist_ok=True)

if __name__ == '__main__':
    print('[trait] %s  pheno=%s  out=%s' % (TRAIT, C.PHENO, F.OUT), flush=True)
    F.main()
