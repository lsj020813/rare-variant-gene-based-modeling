
import re
s = open("main.py").read()
old = """def burden(m):
    S = torch.zeros(NG, NS, device=DEV)
    S.index_add_(0, GID, (m[:, None] * Dh.float()))
    return S"""
new = """def burden(m, blk=4096):
    # block-wise: Dh stays fp16; only a blk x NS fp32 slice materializes (~1.4GB at 4096)
    S = torch.zeros(NG, NS, device=DEV)
    V = Dh.shape[0]
    for a in range(0, V, blk):
        b = min(a + blk, V)
        S.index_add_(0, GID[a:b], m[a:b, None] * Dh[a:b].float())
    return S"""
assert s.count(old) == 1
open("main.py","w").write(s.replace(old, new))
import ast; ast.parse(open("main.py").read())
print("patched")
