#!/usr/bin/env python3
import sys, os, subprocess, json, time
thr = os.environ.get('LD_THREADS', '4')
for k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[k] = thr
import numpy as np
vcf, keep, out = sys.argv[1:4]
step2 = sys.argv[4] if len(sys.argv) > 4 else None
allow = None
if step2:
    with open(step2) as f:
        hdr = f.readline().rstrip('\n').split('\t'); ic, ip, ia1, ia2 = hdr.index('CHR'), hdr.index('POS'), hdr.index('Allele1'), hdr.index('Allele2')
        allow = set()
        for line in f:
            t = line.rstrip('\n').split('\t'); allow.add(f"{t[ic].removeprefix('chr')}:{t[ip]}:{t[ia1]}:{t[ia2]}")
BCF = 'bcftools'
t0 = time.time()
samples = subprocess.check_output([BCF, 'query', '-S', keep, '-l', vcf]).decode().split()
N = len(samples); assert N > 1000, N
nvar = int(subprocess.check_output([BCF, 'index', '-n', vcf]).decode().strip())
X = np.empty((nvar, N), dtype=np.float32)
keys = []; n_slow = 0; n_missing = 0; i = 0; nread = 0
p = subprocess.Popen([BCF, 'query', '-S', keep, '-f', '%CHROM:%POS:%REF:%ALT[\t%DS]\n', vcf], stdout=subprocess.PIPE, bufsize=1 << 22)
for line in p.stdout:
    key, rest = line.rstrip(b'\n').split(b'\t', 1)
    arr = np.fromstring(rest, dtype=np.float32, sep='\t')
    if arr.shape[0] != N:
        n_slow += 1
        toks = rest.split(b'\t'); assert len(toks) == N, (len(toks), N)
        arr = np.array([float(t) if t not in (b'.', b'') else np.nan for t in toks], dtype=np.float32)
        m = np.isnan(arr); n_missing += int(m.sum()); arr[m] = np.nanmean(arr)
    nread += 1
    if allow is not None and key.decode() not in allow:
        continue
    assert i < nvar, 'index -n 보다 레코드가 많음'
    X[i] = arr; keys.append(key.decode()); i += 1
rc = p.wait(); assert rc == 0, f'bcftools query rc={rc}'
assert nread == nvar, f'레코드 수 {nread} != index {nvar}'
if allow is not None:
    assert i == len(allow), f'step2 변이 {len(allow)} 중 VCF 에서 {i} 만 매칭'
    X = X[:i]
t_read = time.time() - t0
mu = X.mean(axis=1, dtype=np.float64).astype(np.float32)
X -= mu[:, None]
ss = np.einsum('ij,ij->i', X, X, dtype=np.float64)
sd = np.sqrt(ss / N)
keepv = sd > 0
n_mono = int((~keepv).sum())
if n_mono:
    X = X[keepv]; keys = [k for k, f in zip(keys, keepv) if f]; sd = sd[keepv]
X /= sd[:, None].astype(np.float32)
R = (X @ X.T) / np.float32(N)
del X
d = np.diag(R).astype(np.float64); max_diag_dev = float(np.abs(d - 1).max())
R = (R + R.T) * np.float32(0.5)
np.fill_diagonal(R, 1.0)
np.clip(R, -1.0, 1.0, out=R)
assert np.isfinite(R).all(), 'R 에 비유한값'
R.astype(np.float32, copy=False).tofile(out + '.ld.bin')
with open(out + '.ld.vars', 'w') as f:
    f.write('\n'.join(keys) + '\n')
info = dict(n_samples=N, n_var_in=nvar, n_var_step2=(len(allow) if allow is not None else None), n_var_ld=len(keys), n_mono_dropped=n_mono, n_slow_lines=n_slow, n_missing_ds=n_missing,
            max_diag_dev_prefix=max_diag_dev, read_s=round(t_read, 1), total_s=round(time.time() - t0, 1))
json.dump(info, open(out + '.ld.json', 'w'))
print('LD_OK', json.dumps(info))
