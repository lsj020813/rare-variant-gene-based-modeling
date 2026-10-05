import torch, torch.nn as nn, time
torch.manual_seed(0)
N_ANNO = 15
PHI_H = 20
RHO_H = 10
BATCH = 1024
MAX_VAR = 200

class DeepSet(nn.Module):

    def __init__(s):
        super().__init__()
        s.phi = nn.Sequential(nn.Linear(N_ANNO, PHI_H), nn.LeakyReLU(), nn.Linear(PHI_H, PHI_H), nn.LeakyReLU())
        s.rho = nn.Sequential(nn.Linear(PHI_H, RHO_H), nn.LeakyReLU(), nn.Linear(RHO_H, RHO_H), nn.LeakyReLU(), nn.Linear(RHO_H, 1))

    def forward(s, x):
        h = s.phi(x)
        h, _ = h.max(dim=1)
        return s.rho(h)

def bench(dev, iters=200):
    m = DeepSet().to(dev)
    opt = torch.optim.AdamW(m.parameters())
    x = torch.randn(BATCH, MAX_VAR, N_ANNO, device=dev)
    y = torch.randn(BATCH, 1, device=dev)
    lossf = nn.MSELoss()
    for _ in range(10):
        opt.zero_grad()
        l = lossf(m(x), y)
        l.backward()
        opt.step()
    if dev == 'cuda':
        torch.cuda.synchronize()
    t = time.time()
    for _ in range(iters):
        opt.zero_grad()
        l = lossf(m(x), y)
        l.backward()
        opt.step()
    if dev == 'cuda':
        torch.cuda.synchronize()
    return (time.time() - t) / iters * 1000
torch.set_num_threads(8)
cpu8 = bench('cpu')
torch.set_num_threads(1)
cpu1 = bench('cpu')
gpu = bench('cuda')
print(f'model params: {sum((p.numel() for p in DeepSet().parameters())):,}')
print(f'batch {BATCH} x {MAX_VAR} variants x {N_ANNO} anno')
print(f'CPU  1-thread : {cpu1:7.2f} ms/step')
print(f'CPU  8-thread : {cpu8:7.2f} ms/step')
print(f'GPU  A2       : {gpu:7.2f} ms/step')
print(f'speedup GPU vs 1-thread CPU: {cpu1 / gpu:.1f}x')
print(f'speedup GPU vs 8-thread CPU: {cpu8 / gpu:.1f}x')
print(f'GPU mem allocated: {torch.cuda.max_memory_allocated() / 1000000.0:.0f} MB')
