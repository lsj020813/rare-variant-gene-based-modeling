
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import numpy as np, os, time, resource, sys
R=_config_path("${PROJECT_ROOT}/work/ref15")
def rss(): return int(open("/proc/self/status").read().split("VmRSS:")[1].split()[0])/1048576
chrs=[str(i) for i in range(1,23)]
metas=[]
for c in chrs:
    z=np.load(f"{R}/annot/ds/chr{c}.ds.npz", allow_pickle=True); metas.append((c,len(z["keys"]),int(z["indptr"][-1]))); del z
NV,NNZ=sum(m[1] for m in metas),sum(m[2] for m in metas)
print(f"NV {NV:,} NNZ {NNZ:,} expected data+idx {NNZ*8/1e9:.1f} GB | rss after metas {rss():.1f}G", flush=True)
data=np.empty(NNZ,np.float32); indices=np.empty(NNZ,np.int32); indptr=np.empty(NV+1,np.int64); indptr[0]=0
print(f"after prealloc rss {rss():.1f}G (untouched pages)", flush=True)
off=noff=0; t0=time.time()
for c,nv,nnz in metas:
    z=np.load(f"{R}/annot/ds/chr{c}.ds.npz", allow_pickle=True)
    d=z["data"]; r1=rss(); data[noff:noff+nnz]=d; del d
    ix=z["indices"]; r2=rss(); indices[noff:noff+nnz]=ix; del ix
    indptr[off+1:off+nv+1]=z["indptr"][1:]+noff; off+=nv; noff+=nnz; del z
    print(f"chr{c:>2} nnz {nnz/1e9:.2f}G  rss: after data-temp {r1:.1f}G  after idx-temp {r2:.1f}G  end {rss():.1f}G  filled {noff*8/1e9:.1f}G  t={time.time()-t0:.0f}s", flush=True)
    if rss() > 150: print("STOP >150G"); sys.exit(3)
print("LOAD_PROBE_DONE", f"final rss {rss():.1f}G")
