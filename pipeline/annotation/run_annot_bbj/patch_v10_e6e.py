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
if "E6E_PATCH" in src: print("already patched"); sys.exit(0)
assert "E6D_PATCH" in src
shutil.copy(P, P + ".pre_e6e_bak")
def rep(old, new):
    global src
    assert src.count(old) == 1, f"count {src.count(old)}: {old[:70]}"
    src = src.replace(old, new)
rep("        self.A = self.basis_for(self.train)\n",
    "        # E6E_PATCH: keep the float32 basis rows (exact values of design.B); products run in float64 per row chunk\n"
    "        if self.group_names == ['shared']:\n"
    "            self.A = self.design.B[self.train if self.control else self.data.vi[self.train]]\n"
    "        else:\n"
    "            self.A = self.basis_for(self.train).astype(np.float32)\n"
    "        self._step = max(1, args.design_chunk * 16)\n")
rep("    def basis_for(self, indices):\n",
    "    def _Av(self, coef):\n"
    "        \"\"\"E6E_PATCH: A @ coef with float64 arithmetic, chunked over rows.\"\"\"\n"
    "        out = np.empty(self.A.shape[0])\n"
    "        for lo in range(0, self.A.shape[0], self._step):\n"
    "            out[lo:lo + self._step] = self.A[lo:lo + self._step].astype(np.float64) @ coef\n"
    "        return out\n\n"
    "    def _Atv(self, v):\n"
    "        \"\"\"E6E_PATCH: A.T @ v with float64 arithmetic, chunked over rows.\"\"\"\n"
    "        out = np.zeros(self.A.shape[1])\n"
    "        for lo in range(0, self.A.shape[0], self._step):\n"
    "            out += self.A[lo:lo + self._step].astype(np.float64).T @ v[lo:lo + self._step]\n"
    "        return out\n\n"
    "    def basis_for(self, indices):\n")
rep("        eta = self.A @ a[:self.ncoef] + a[self.ncoef + self.tidx]  # E6C_PATCH\n",
    "        eta = self._Av(a[:self.ncoef]) + a[self.ncoef + self.tidx]  # E6C/E6E_PATCH\n")
rep("        grad = np.concatenate([self.A.T @ resid, np.bincount(self.tidx, weights=resid, minlength=len(self.trait_ids))])\n",
    "        grad = np.concatenate([self._Atv(resid), np.bincount(self.tidx, weights=resid, minlength=len(self.trait_ids))])\n")
rep("        basis = self.A  # E6C_PATCH: A is basis-only\n        rms = np.sqrt(np.mean(np.square(basis @ pg)))\n",
    "        rms = np.sqrt(np.mean(np.square(self._Av(pg))))  # E6C/E6E_PATCH\n")
rep("        p0 = expit(self.A @ self.initial[:self.ncoef] + self.initial[self.ncoef + self.tidx])  # E6C_PATCH chunked Hessian\n",
    "        p0 = expit(self._Av(self.initial[:self.ncoef]) + self.initial[self.ncoef + self.tidx])  # E6C/E6E_PATCH chunked Hessian\n")
rep("        step = max(1, self.args.design_chunk * 16)\n        for lo in range(0, len(wpp), step):\n            Ac = self.A[lo:lo + step]; WA = wpp[lo:lo + step, None] * Ac\n",
    "        step = self._step\n        for lo in range(0, len(wpp), step):\n            Ac = self.A[lo:lo + step].astype(np.float64); WA = wpp[lo:lo + step, None] * Ac\n")
open(P, "w").write(src); py_compile.compile(P, doraise=True)
print("PATCHED5", hashlib.sha256(open(P,'rb').read()).hexdigest()[:16], src.count("E6E_PATCH"))
