import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


SRC = _config_path('${PROJECT_ROOT}/work/fset/att/')
OLD = 'import att_common as C'
NEW = 'import att_common_b200 as C'

for name in ('att_labels', 'att_burden', 'att_consolidate'):
    s = open(SRC + name + '.py').read()
    assert s.count(OLD) == 1, name
    open(SRC + name + '_b200.py', 'w').write(s.replace(OLD, NEW))
    print('wrote %s_b200.py' % name)

s = open(SRC + 'att_fit_cvmin.py').read()
assert s.count(OLD) == 1
s = s.replace(OLD, NEW)

old_out = """OUT = C.ATT + '/cvmin'
PRIVOUT = OUT + '/private'"""
new_out = """OUT = C.ATT
PRIVOUT = C.PRIV"""
assert s.count(old_out) == 1
s = s.replace(old_out, new_out)

old_arg = "    nshuf = int(sys.argv[1]) if len(sys.argv) > 1 else C.NSHUF"
new_arg = """    nshuf = int(sys.argv[1]) if len(sys.argv) > 1 else C.NSHUF
    L0 = int(sys.argv[2]) if len(sys.argv) > 3 else 1
    L1 = int(sys.argv[3]) if len(sys.argv) > 3 else nshuf
    is_part = len(sys.argv) > 3
    print('[part] L %d..%d (is_part=%s)' % (L0, L1, is_part), flush=True)"""
assert s.count(old_arg) == 1
s = s.replace(old_arg, new_arg)

old_loop = "    for L in range(1, nshuf + 1):"
assert s.count(old_loop) == 1
s = s.replace(old_loop, "    for L in range(L0, L1 + 1):")

old_csv = "    with open(OUT + '/att_label_shuffle_null.csv', 'w') as fh:"
new_csv = ("    _nullcsv = (OUT + '/null_part_%03d_%03d.csv' % (L0, L1)) \\\n"
           "        if is_part else (OUT + '/att_label_shuffle_null.csv')\n"
           "    with open(_nullcsv, 'w') as fh:")
assert s.count(old_csv) == 1
s = s.replace(old_csv, new_csv)

for path in ("OUT + '/att_gene_results.csv'", "OUT + '/att_weights.csv'",
             "OUT + '/att_fit_meta.json'", "OUT + '/att_fit.done'"):
    old = "    with open(%s, 'w') as fh:" % path
    new = ("    if L0 == 1:\n"
           "      with open(%s, 'w') as fh:" % path)
    assert s.count(old) == 1, path
    s = s.replace(old, new)

lines = s.split('\n')
out = []
i = 0
while i < len(lines):
    out.append(lines[i])
    if lines[i].startswith('    if L0 == 1:'):
        out.append(lines[i + 1])
        i += 2
        while i < len(lines) and (lines[i].startswith('        ')
                                  or lines[i].strip() == ''):
            if lines[i].strip() == '' and i + 1 < len(lines) \
                    and not lines[i + 1].startswith('        '):
                break
            out.append('  ' + lines[i] if lines[i].strip() else lines[i])
            i += 1
        continue
    i += 1
s = '\n'.join(out)

s = s.replace("amendment='A5-CVMIN'", "amendment='A5-CVMIN-B200'")
s = s.replace("base_run='fset/att (1-SE, 2026-09-23 15:02)'",
              "base_run='fset/att/cvmin (CV-min, 20 shuffles)'")

open(SRC + 'att_fit_b200.py', 'w').write(s)
print('wrote att_fit_b200.py')

import ast
ast.parse(open(SRC + 'att_fit_b200.py').read())
for name in ('att_labels', 'att_burden', 'att_consolidate'):
    ast.parse(open(SRC + name + '_b200.py').read())
print('문법 검사 OK')
