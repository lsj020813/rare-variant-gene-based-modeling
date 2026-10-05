import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import gzip
import io
import json
import zipfile
import numpy as np
import att_common as C

ZIP = _config_path('${PROJECT_ROOT}/work/ref/bbj_fm/hum0197.v5.finemap.TC.v1.zip')
MEMBER = 'BBJ.TC.Kanai2021.FINEMAP.tsv.gz'
PIP_MIN = 0.1
WIN = 500000

def main():
    hits = {}
    with zipfile.ZipFile(ZIP) as z:
        with z.open(MEMBER) as fh:
            gz = gzip.GzipFile(fileobj=io.BytesIO(fh.read()))
            head = gz.readline().decode().rstrip('\n').split('\t')
            I = dict((c, i) for i, c in enumerate(head))
            ci = I.get('chromosome', 0)
            pi = I.get('position', 1)
            qi = I.get('pip', 11)
            for line in gz:
                p = line.decode().rstrip('\n').split('\t')
                try:
                    if float(p[qi]) < PIP_MIN:
                        continue
                    hits.setdefault(p[ci].replace('chr', ''), []).append(
                        int(p[pi]))
                except (ValueError, IndexError):
                    continue
    for k in hits:
        hits[k] = np.sort(np.asarray(hits[k], dtype=np.int64))
    d = C.master()
    di = np.asarray(d['dom_idx'])
    p19 = np.asarray(d['pos19'])
    udom = np.asarray(d['udom'])
    dch = np.asarray(d['dom_chrom'])
    rows = []
    for g in range(C.NDOM):
        s = p19[di == g]
        lo, hi = int(s.min()) - WIN, int(s.max()) + WIN
        ch = dch[g].replace('chr', '')
        arr = hits.get(ch)
        n = 0 if arr is None else int(
            np.searchsorted(arr, hi, 'right') - np.searchsorted(arr, lo))
        rows.append((g, str(udom[g]), ch, int(np.median(s)), n, int(n > 0)))
    with open(C.ATT + '/att_bbj_tc_domains.csv', 'w') as f:
        f.write('domain_idx,gene,chrom,center_hg19,'
                'n_bbj_tc_pip01_in_window,bbj_tc_domain\n')
        for r in rows:
            f.write(','.join(str(x) for x in r) + '\n')
    meta = dict(source=ZIP, member=MEMBER, pip_min=PIP_MIN, window_bp=WIN,
                rule=('도메인 U-b 변이 span ±500 kb 안에 BBJ TC FINEMAP '
                      'PIP>=0.1 변이 1개 이상'),
                n_bbj_variants_pip01=int(sum(len(v) for v in hits.values())),
                n_domains_flagged=int(sum(r[5] for r in rows)))
    with open(C.ATT + '/att_bbj_meta.json', 'w') as f:
        json.dump(meta, f, indent=1, ensure_ascii=False)
    print(json.dumps(meta, indent=1, ensure_ascii=False))

if __name__ == '__main__':
    main()
