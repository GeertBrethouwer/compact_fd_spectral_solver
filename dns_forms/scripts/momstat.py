import numpy as np, re, os
def ld(c):
    rows=[]
    for l in open(c+'/data/log.dat'):
        if l.startswith('#'): continue
        rows.append([float(re.sub(r'(\d)([+-]\d{3})$', r'\1E\2', x)) for x in l.split()])
    return np.array(rows)
def stats(c):
    # profiles from stats.dat (written at the end), averaged from tstats_start
    f=c+'/data/stats.dat'
    if not os.path.exists(f): return None
    hdr=[l for l in open(f) if l.startswith('#')]
    d=np.loadtxt(f); return hdr,d
print('Turbulent channel, 32x65x32, Re_b = 2800, restarted from developed turbulence (rotational/D1)')
print('alpha_s = 4: t = 0..40, alpha_s = 5: t = 0..30 (from the alpha_s = 4 field interpolated in y); averages over t >= 10')
print('max div: over t > 0 (the alpha_s = 5 start field, interpolated in y, has div = 2e-2 at t = 0 only)')
print('U_c+, u_rms+, -uv+: from stats.dat (u_tau from the averaged Re_tau), centreline and peak values')
print(f'{"form":11s} {"d/dy":4s} | {"a":>2s} {"status":>18s} {"<Re_tau>":>9s} {"<Ekin>":>9s} {"max div":>9s} | {"U_c+":>6s} {"urms+":>6s} {"vrms+":>6s} {"wrms+":>6s} {"-uv+":>6s}')
for a in (4,5):
  for nl in ('rotational','advective','divergence','skew'):
    for dy in ('d1','d1m'):
        c=f'mom_a{a}_{nl}_{dy}'; d=ld(c); t=d[:,1]; ek=d[:,8]
        bad=~np.isfinite(ek)|(ek>1.)
        st='blow-up at t=%.1f'%t[bad][0] if bad.any() else 'stable to t=%.0f'%t.max()
        m=(t>=10)&~bad
        rt=d[m,6].mean() if m.any() and not bad.any() else np.nan; e=ek[m].mean() if m.any() and not bad.any() else np.nan
        sp=''
        if not bad.any():
            ut=rt/2800.; p=np.loadtxt(c+'/data/stats.dat'); ny=len(p); h=p[:(ny+1)//2]
            sp=f' | {p[ny//2,1]/ut:6.2f} {np.sqrt(h[:,4].max())/ut:6.2f} {np.sqrt(h[:,5].max())/ut:6.2f} {np.sqrt(h[:,6].max())/ut:6.2f} {(-h[:,7]).max()/ut**2:6.3f}'
        print(f'{nl:11s} {dy:4s} | {a:2d} {st:>18s} {rt:9.1f} {e:9.2e} {d[(~bad)&(t>0),7].max():9.1e}'+sp)
