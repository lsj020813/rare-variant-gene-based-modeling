#!/usr/bin/env python3
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import shutil, sys, hashlib, py_compile
P = _config_path("${PROJECT_ROOT}/work/run_band15/model_v10_out/l1_train_v10.py")
src = open(P).read()
if "E6C_PATCH" in src: print("already patched"); sys.exit(0)
assert "E6B_PATCH" in src
shutil.copy(P, P + ".pre_e6c_bak")
def rep(old, new):
    global src
    assert src.count(old) == 1, f"count {src.count(old)}: {old[:70]}"
    src = src.replace(old, new)
rep("        self.A = self.matrix(self.train)\n",
    "        # E6C_PATCH: basis-only float64 matrix; trait intercepts implicit via row->trait position (identical loss/gradient/Hessian)\n"
    "        self.A = self.basis_for(self.train)\n"
    "        self.tidx = np.searchsorted(np.asarray(self.trait_ids), data.ti[self.train])\n"
    "        require(np.array_equal(np.asarray(self.trait_ids)[self.tidx], data.ti[self.train]), 'trait index mapping')\n")
rep("        eta = self.A @ a\n        loss = np.sum(self.w * (np.logaddexp(0., eta) - self.y[self.train] * eta))\n        grad = self.A.T @ (self.w * (expit(eta) - self.y[self.train]))\n",
    "        eta = self.A @ a[:self.ncoef] + a[self.ncoef + self.tidx]  # E6C_PATCH\n"
    "        loss = np.sum(self.w * (np.logaddexp(0., eta) - self.y[self.train] * eta))\n"
    "        resid = self.w * (expit(eta) - self.y[self.train])\n"
    "        grad = np.concatenate([self.A.T @ resid, np.bincount(self.tidx, weights=resid, minlength=len(self.trait_ids))])\n")
rep("        basis = self.A[:, :self.ncoef]\n", "        basis = self.A  # E6C_PATCH: A is basis-only\n")
rep("        p0 = expit(self.A @ self.initial)\n        hessian = self.A.T @ ((self.w * p0 * (1 - p0))[:, None] * self.A)\n        hessian[:self.ncoef, :self.ncoef] += 2 * lam * self.P\n",
    "        p0 = expit(self.A @ self.initial[:self.ncoef] + self.initial[self.ncoef + self.tidx])  # E6C_PATCH chunked Hessian\n"
    "        wpp = self.w * p0 * (1 - p0)\n"
    "        q, n_t = self.ncoef, len(self.trait_ids)\n"
    "        hessian = np.zeros((q + n_t, q + n_t))\n"
    "        step = max(1, self.args.design_chunk * 16)\n"
    "        for lo in range(0, len(wpp), step):\n"
    "            Ac = self.A[lo:lo + step]; WA = wpp[lo:lo + step, None] * Ac\n"
    "            hessian[:q, :q] += Ac.T @ WA\n"
    "            onehot = (self.tidx[lo:lo + step, None] == np.arange(n_t)[None, :]).astype(float)\n"
    "            hessian[q:, :q] += onehot.T @ WA\n"
    "        hessian[:q, q:] = hessian[q:, :q].T\n"
    "        hessian[q:, q:] = np.diag(np.bincount(self.tidx, weights=wpp, minlength=n_t))\n"
    "        hessian[:self.ncoef, :self.ncoef] += 2 * lam * self.P\n")
src = src.replace("data = Data(args)\n", "data = Data(args); _trim_memory()  # E6C_PATCH\n")
rep("class Trainer:\n",
    "def _trim_memory():\n    \"\"\"E6C_PATCH: return parser garbage to the OS after Data load (no analysis effect).\"\"\"\n"
    "    import gc, ctypes\n    gc.collect()\n    try:\n        ctypes.CDLL('libc.so.6').malloc_trim(0)\n    except Exception:\n        pass\n\n\nclass Trainer:\n")
open(P, "w").write(src); py_compile.compile(P, doraise=True)
print("PATCHED3", hashlib.sha256(open(P,'rb').read()).hexdigest()[:16], src.count("E6C_PATCH"))
