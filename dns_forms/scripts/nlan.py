import numpy as np, re
def ld(c):
    rows=[]
    for l in open(c+'/data/log.dat'):
        if l.startswith('#'): continue
        v=[float(re.sub(r'(\d)([+-]\d{3})$', r'\1E\2', x)) for x in l.split()]
        rows.append(v)
    d=np.array(rows); return d[:,1],d[:,9],d[:,8],d[:,5]
ref=ld('nlref_a2_385')
def err(c, r=ref):
    t,em,ek,tw=ld(c)
    bad=~np.isfinite(ek) | (ek>1e3*r[2].max())
    if bad.any() or t.max()<49.9:
        tb=t[bad][0] if bad.any() else t.max()
        return f'{"blow-up at t="+format(tb,".1f"):>29s}'
    i=lambda y: np.interp(r[0], t, y); m=r[0]<=50
    e=[abs(i(em)-r[1])[m].max()/r[1][m].max(), abs(i(ek)-r[2])[m].max()/r[2][m].max(), abs(i(tw)-r[3])[m].max()/r[3][m].max()]
    return ' '.join(f'{x:9.1e}' for x in e)
print('Finite-amplitude TS wave (Re = 5000, alpha = 1.12, amplitude 0.05, nx = 16, 0 <= t <= 50)')
print('max over t of the relative deviation from the reference (Ny = 385, tanh alpha_s = 2, rotational, D1)')
print('columns per grid: E_TS (energy of the TS mode), Ekin (fluctuation energy), tau_w (wall shear)')
print('reference uncertainty: Ny=257 a2:', err('nlref_a2_257'), '| Ny=257 a4:', err('nlref_a4_257'))
print(f'{"form":11s} {"d/dy":4s} | ' + ' | '.join(f'{f"a{a} Ny={n}":>29s}' for a in (4,5) for n in (33,65)))
for nl in ('rotational','advective','divergence','skew'):
  for dy in ('d1','d1m'):
    print(f'{nl:11s} {dy:4s} | '+' | '.join(err(f'nl_a{a}_{n}_{nl}_{dy}') for a in (4,5) for n in (33,65)))
t,em,ek,tw=ref; m=t<=50; print('reference: E_TS(0) = %.3e, E_TS(50) = %.3e, Ekin(50) = %.3e'%(em[0],em[m][-1],ek[m][-1]))
