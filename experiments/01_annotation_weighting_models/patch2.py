
s = open("main.py").read()

old_burden = """def burden(m, blk=4096):
    # block-wise: Dh stays fp16; only a blk x NS fp32 slice materializes (~1.4GB at 4096)
    S = torch.zeros(NG, NS, device=DEV)
    V = Dh.shape[0]
    for a in range(0, V, blk):
        b = min(a + blk, V)
        S.index_add_(0, GID[a:b], m[a:b, None] * Dh[a:b].float())
    return S"""
new_burden = """BLK = 2048
def burden_nograd(m, cols=None):
    # forward only — no autograd storage. cols: optional sample subset.
    with torch.no_grad():
        ns = NS if cols is None else len(cols)
        S = torch.zeros(NG, ns, device=DEV)
        for a in range(0, Dh.shape[0], BLK):
            b = min(a + BLK, Dh.shape[0])
            Db = Dh[a:b] if cols is None else Dh[a:b][:, cols]
            S.index_add_(0, GID[a:b], m[a:b, None] * Db.float())
        return S

def grad_m_from_gS(gS, cols=None):
    # dL/dm_v = sum_i gS[gene(v), i] * D[v, i]  — block-wise, no big storage
    with torch.no_grad():
        gm = torch.zeros(Dh.shape[0], device=DEV)
        for a in range(0, Dh.shape[0], BLK):
            b = min(a + BLK, Dh.shape[0])
            Db = Dh[a:b] if cols is None else Dh[a:b][:, cols]
            gm[a:b] = (gS[GID[a:b]] * Db.float()).sum(1)
        return gm"""
assert s.count(old_burden) == 1
s = s.replace(old_burden, new_burden)

old_train = """def train_arm(phi_fn, params, train_gi, lam, tag, fold):
    opt = torch.optim.AdamW(params, lr=5e-3)
    for ep in range(EPOCHS):
        opt.zero_grad()
        for t in TRAITS:
            m = phi_fn()
            S = zrow(burden(m)[train_gi][:, tr[t]["cols"]])
            resid = tr[t]["resid"]
            G = S @ S.T + lam * torch.eye(len(train_gi), device=DEV)
            beta = torch.linalg.solve(G, S @ resid)
            ((resid - beta @ S) ** 2).mean().backward()
            del m, S, G, beta
            torch.cuda.empty_cache()
        opt.step()
        if ep % 100 == 0: print(f"[f{fold} {tag} lam{lam}] ep {ep}", flush=True)"""
new_train = """def train_arm(phi_fn, params, train_gi, lam, tag, fold):
    # manual chain rule: big products run under no_grad; S is re-leafed; autograd
    # covers only (i) the tiny phi graph X->m and (ii) the S->loss head.
    opt = torch.optim.AdamW(params, lr=5e-3)
    for ep in range(EPOCHS):
        opt.zero_grad()
        for t in TRAITS:
            cols, resid = tr[t]["cols"], tr[t]["resid"]
            m = phi_fn()                                   # small graph, keep
            S0 = burden_nograd(m.detach(), cols)[train_gi] # [Gtr, n] fp32, no graph
            S_leaf = S0.detach().requires_grad_(True)
            S = zrow(S_leaf)
            G = S @ S.T + lam * torch.eye(len(train_gi), device=DEV)
            beta = torch.linalg.solve(G, S @ resid)
            loss = ((resid - beta @ S) ** 2).mean()
            loss.backward()                                # -> S_leaf.grad [Gtr, n]
            gS_full = torch.zeros(NG, len(cols), device=DEV)
            gS_full[train_gi] = S_leaf.grad
            gm = grad_m_from_gS(gS_full, cols)             # [V]
            m.backward(gradient=gm)                        # -> phi params
            del m, S0, S_leaf, S, G, beta, loss, gS_full, gm
            torch.cuda.empty_cache()
        opt.step()
        if ep % 100 == 0: print(f"[f{fold} {tag} lam{lam}] ep {ep}", flush=True)"""
assert s.count(old_train) == 1
s = s.replace(old_train, new_train)

s = s.replace("S_all = burden(m)", "S_all = burden_nograd(m)")
open("main.py","w").write(s)
import ast; ast.parse(s)
print("rewritten — manual chain rule, syntax OK")
