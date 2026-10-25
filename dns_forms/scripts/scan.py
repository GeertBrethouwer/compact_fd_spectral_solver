import numpy as np, re, sys, os
def ld(c):
    rows=[]
    for l in open(c+'/data/log.dat'):
        if l.startswith('#'): continue
        rows.append([float(re.sub(r'(\d)([+-]\d{3})$', r'\1E\2', x)) for x in l.split()])
    return np.array(rows)
a=sys.argv[1]
print(f'Scalars in turbulent channel, tanh alpha_s = {a}, 32x65x32, Re_tau ~ 200, momentum rotational/D1, t = 0..20')
print('s1: Pr = 0.71, s = +1/-1 at the walls, linear initially: wall gradient at t = 20 (bottom, top), and the max over t of its')
print('    relative deviation from skew/D1m (all variants start from the same field; s1 is still developing)')
print('s2: inviscid (Pr = 1e10), s = 0 at the walls, noise initially: <s2^2>(t)/<s2^2>(0), max over t and value at t = 20')
print(f'{"scal_form":10s} {"d/dy":4s} | {"dS1/dy bottom":>14s} {"top":>9s} {"max dev":>9s} | {"s2 var max":>10s} {"at t=20":>9s}')
ref=ld(f'scal_a{a}_skew_d1m')
for f,sd in (('advective','d1'),('advective','d1m'),('divergence','d1'),('divergence','d1m'),('skewd1','d1'),('skew','d1m')):
    c=f'scal_a{a}_{f}_{sd}'
    if not os.path.exists(c+'/data/log.dat'): continue
    d=ld(c); t=d[:,1]; m=(t>5)
    gb,gt,v2=d[:,10],d[:,11],d[:,15]
    ok=np.isfinite(v2)&(v2<1e30)
    tb='' if ok.all() and t.max()>19.9 else f'  blow-up/stopped at t={t[~ok][0] if (~ok).any() else t.max():.1f}'
    n=min(len(d),len(ref)); dev=max(abs(d[:n,10]/ref[:n,10]-1).max(),abs(d[:n,11]/ref[:n,11]-1).max())
    print(f'{("skew" if f.startswith("skew") else f):10s} {sd:4s} | {gb[-1]:14.5f} {gt[-1]:9.5f} {dev:9.1e} | {(v2[ok]/v2[0]).max():10.3g} {v2[ok][-1]/v2[0]:9.3g}{tb}')
