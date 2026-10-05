
s = open("main.py").read()

old_load = """gv = collections.defaultdict(list)"""
new_load = """CACHE = f"{R}/l0/l0_tensors.pt"
gv = collections.defaultdict(list)"""
assert s.count(old_load) == 1
s = s.replace(old_load, new_load, 1)

old_t = """X = torch.tensor(np.stack(X_rows), device=DEV)                        # [V,11] fp32
Dh = torch.tensor(np.stack(ds_rows), device=DEV)                      # [V,N] fp16 ~6.5GB
GID = torch.tensor(gene_id, device=DEV)"""
new_t = """X = torch.tensor(np.stack(X_rows), device=DEV)                        # [V,11] fp32
Dh = torch.tensor(np.stack(ds_rows), device=DEV)                      # [V,N] fp16 ~6.5GB
GID = torch.tensor(gene_id, device=DEV)
torch.save({"X": X.cpu(), "Dh": Dh.cpu(), "GID": GID.cpu()}, CACHE)
print("tensor cache saved", flush=True)"""
assert s.count(old_t) == 1
s = s.replace(old_t, new_t)

old_ds = """DS = {}
for ch, keys in byc.items():"""
new_ds = """import os
USE_CACHE = os.path.exists(CACHE)
DS = {}
for ch, keys in (byc.items() if not USE_CACHE else []):"""
assert s.count(old_ds) == 1
s = s.replace(old_ds, new_ds)
old_g = """assert len(DS) == len(need), f"GATE FAIL dosage {len(DS)}/{len(need)}"
print(f"dosage: {len(DS):,} variants x configured samples fp16", flush=True)"""
new_g = """if not USE_CACHE:
    assert len(DS) == len(need), f"GATE FAIL dosage {len(DS)}/{len(need)}"
    print(f"dosage: {len(DS):,} variants x configured samples fp16", flush=True)
else:
    print("loading tensor cache", flush=True)"""
assert s.count(old_g) == 1
s = s.replace(old_g, new_g)

old_build = """X_rows, ds_rows, gene_id = [], [], []
for gi, g in enumerate(genes_all):
    for k, sv in gv[g]:
        X_rows.append(sv); ds_rows.append(DS[k]); gene_id.append(gi)"""
new_build = """if USE_CACHE:
    _c = torch.load(CACHE)
    X, Dh, GID = _c["X"].to(DEV), _c["Dh"].to(DEV), _c["GID"].to(DEV)
X_rows, ds_rows, gene_id = [], [], []
for gi, g in enumerate(genes_all):
    if USE_CACHE: break
    for k, sv in gv[g]:
        X_rows.append(sv); ds_rows.append(DS[k]); gene_id.append(gi)"""
assert s.count(old_build) == 1
s = s.replace(old_build, new_build)
old_t2 = """X = torch.tensor(np.stack(X_rows), device=DEV)"""
s = s.replace(old_t2, """if not USE_CACHE:
    X = torch.tensor(np.stack(X_rows), device=DEV)""", 1)
s = s.replace("""    X = torch.tensor(np.stack(X_rows), device=DEV)                        # [V,11] fp32
Dh = torch.tensor(np.stack(ds_rows), device=DEV)                      # [V,N] fp16 ~6.5GB
GID = torch.tensor(gene_id, device=DEV)
torch.save({"X": X.cpu(), "Dh": Dh.cpu(), "GID": GID.cpu()}, CACHE)
print("tensor cache saved", flush=True)""",
"""    X = torch.tensor(np.stack(X_rows), device=DEV)                        # [V,11] fp32
    Dh = torch.tensor(np.stack(ds_rows), device=DEV)                      # [V,N] fp16 ~6.5GB
    GID = torch.tensor(gene_id, device=DEV)
    torch.save({"X": X.cpu(), "Dh": Dh.cpu(), "GID": GID.cpu()}, CACHE)
    print("tensor cache saved", flush=True)""")

old_loop = """        for t in TRAITS:
            cols, resid = tr[t]["cols"], tr[t]["resid"]
            m = phi_fn()                                   # small graph, keep
            S0 = burden_nograd(m.detach(), cols)[train_gi] # [Gtr, n] fp32, no graph"""
new_loop = """        for t in TRAITS:
            cols_full, resid_full = tr[t]["cols"], tr[t]["resid"]
            # column minibatch: 30k samples per step, reshuffled per epoch [OOM fix;
            # all configured participants contribute across epochs]
            bidx = torch.randperm(len(cols_full), device=DEV)[:30000]
            cols, resid = cols_full[bidx], resid_full[bidx]
            m = phi_fn()                                   # small graph, keep
            S0 = burden_nograd(m.detach(), cols)[train_gi] # [Gtr, n] fp32, no graph"""
assert s.count(old_loop) == 1
s = s.replace(old_loop, new_loop)

open("main.py","w").write(s)
import ast; ast.parse(s)
print("patched: cache + column minibatch, syntax OK")
