"""
Corbino-Castillo mimetic operators (as implemented in MOLE, https://github.com/csrc-sdsu/mole) on uniform and
hyperbolic-tangent grids, against hichan's compact second derivative D2.

Staggered layout with m cells: faces x_0..x_m, cell centres plus the two wall values (m+2 values).
  G: centres+walls -> faces,  D: faces -> centres (first and last rows zero),  Laplacian L = D G.
Non-uniform grids by the same diagonal metric as MOLE's gradNonUniform/divNonUniform.

MOLE is GPL-3.0, so no MOLE code is included here: the boundary matrices A and the interior stencils are read
at run time from MOLE's gradNonPeriodic.m and divNonPeriodic.m (downloaded into the cache directory).

Usage: python compare_mimetic.py OPSDIR CACHEDIR   (OPSDIR: ops_<grid>_<ny>_o<order>.bin from dump_ops.f90)
"""
import sys, os, re, urllib.request, numpy as np
from fractions import Fraction
ops, cache = sys.argv[1], sys.argv[2]
URL = 'https://raw.githubusercontent.com/csrc-sdsu/mole/main/src/octave/{}.m'
GRIDS = {'uniform': 0, 'tanh2': 2, 'tanh3': 3, 'tanh4': 4}
NYS = [33, 65, 129, 257]

def mole_tables(name):
    """{k: (A, interior stencil)} from MOLE's gradNonPeriodic.m / divNonPeriodic.m"""
    fn = os.path.join(cache, name + '.m')
    if not os.path.exists(fn):
        os.makedirs(cache, exist_ok=True); urllib.request.urlretrieve(URL.format(name), fn)
    src = open(fn).read()
    out = {}
    for k in (4, 6, 8):
        block = re.search(rf'case {k}\b(.*?)(?=case \d|\n\s*end\s*\n)', src, re.S).group(1)
        num = lambda t: float(Fraction(t))
        A = re.search(r'A = \[(.*?)\];', block, re.S).group(1).replace('...', '')
        A = np.array([[num(t) for t in row.split()] for row in A.split(';') if row.strip()])
        sten = re.findall(r'=\s*\[([^\]A]*)\];', block)[-1]
        out[k] = (A, np.array([num(t) for t in sten.split()]))
    return out

GT, DT = mole_tables('gradNonPeriodic'), mole_tables('divNonPeriodic')

def mimetic(k, m):
    """MOLE's G ((m+1) x (m+2)) and D ((m+2) x (m+1)) for unit spacing (index patterns as in MOLE)"""
    A, s = GT[k]; nb = A.shape[0]; G = np.zeros((m+1, m+2))
    G[:nb, :k+1] = A; G[m+1-nb:, m+1-k:] = -np.rot90(A, 2)
    for i in range(nb, m+1-nb):                       # MOLE: i = k/2+1 .. m-k/2+1 (1-based)
        G[i, i-k//2+1:i+k//2+1] = s
    A, s = DT[k]; nb = A.shape[0]; D = np.zeros((m+2, m+1))
    D[1:nb+1, :k+1] = A; D[m+1-nb:m+1, m-k:] = -np.rot90(A, 2)
    for i in range(nb+1, m+1-nb):                     # MOLE: i = k/2+1 .. m-k/2+2 ... (1-based rows nb+2..)
        D[i, i-k//2:i+k//2] = s
    return G, D

def ymap(xi, a):
    return 2*xi-1 if a == 0 else np.tanh(a*(xi-.5))/np.tanh(.5*a)

f = lambda y: np.sin(2.3*y+.4); df = lambda y: 2.3*np.cos(2.3*y+.4); d2f = lambda y: -2.3**2*f(y)

def gauss(G, D):
    """MOLE's weights (weightsP: G^T p = (-1,0..0,1); weightsQ: interior rows of D^T q = (-1,0..0,1), q = 1 at the
    walls) and the boundary operator B = Q D + G^T P of the extended Gauss divergence theorem.
    Returns: B's deviation from the plain boundary term (-1 at (0,0), +1 at (m+1,m)) relative to 1, the depth
    (rows from the wall) to which B is non-zero, and whether all weights are positive."""
    nq, npf = D.shape
    b = np.zeros(nq); b[0], b[-1] = -1, 1
    p = np.linalg.lstsq(G.T, b, rcond=None)[0]
    bq = np.zeros(npf); bq[0], bq[-1] = -1, 1
    q = np.r_[1, np.linalg.lstsq(D[1:-1].T, bq, rcond=None)[0], 1]
    B = q[:, None]*D + G.T*p[None, :]
    B0 = np.zeros_like(B); B0[0, 0], B0[-1, -1] = -1, 1
    dev = abs(B-B0); rows = np.where(dev.max(1) > 1e-10*abs(B).max())[0]
    depth = 0 if len(rows) == 0 else int(min(rows.max(), nq-1-rows.min()) if False else max(min(r, nq-1-r) for r in rows))+1
    return dev.max(), depth, bool((p > 0).all() and (q > 0).all())

def load_d2(grid, ny, order):
    a = np.fromfile(f'{ops}/ops_{grid}_{ny}_o{order}.bin', dtype='<f8', offset=4)
    return a[:ny], a[2*ny+2*ny*ny:].reshape(ny, ny).T

rows = []
for grid, a in GRIDS.items():
    for ny in NYS:
        m = ny-1
        xf = ymap(np.arange(m+1)/m, a); xc = np.r_[xf[0], ymap((np.arange(m)+.5)/m, a), xf[-1]]
        for k in (4, 6, 8):
            G, D = mimetic(k, m)
            G = G/(G @ xc)[:, None]                          # gradNonUniform: J = (G x)^-1
            jd = D @ xf; jd[[0, -1]] = 1; D = D/jd[:, None]  # divNonUniform (wall rows stay zero)
            L = D @ G
            eg = abs(G @ f(xc) - df(xf)).max(); ed = abs((D @ df(xf) - d2f(xc))[1:-1]).max()
            el = abs((L @ f(xc) - d2f(xc))[1:-1]).max()
            lam = np.linalg.eigvals(L[1:-1, 1:-1]); dy = np.diff(xf).min()
            r = dict(op=f'mimetic k={k}', grid=grid, ny=ny, eg=eg, ed=ed, el=el, neg=bool((lam.real < 0).all()),
                     rho=abs(lam).max()*dy**2, im=abs(lam.imag).max()/abs(lam).max())
            r['bdev'], r['bdepth'], r['wpos'] = gauss(G, D)
            rows.append(r)
        for order in (4, 6):
            y, D2 = load_d2(grid, ny, order)
            el = abs((D2 @ f(y) - d2f(y))[1:-1]).max()
            lam = np.linalg.eigvals(D2[1:-1, 1:-1]); dy = np.diff(y).min()
            rows.append(dict(op=f'hichan D2-{order}', grid=grid, ny=ny, eg=np.nan, ed=np.nan, el=el,
                             neg=bool((lam.real < 0).all()), rho=abs(lam).max()*dy**2,
                             im=abs(lam.imag).max()/abs(lam).max()))

for grid in GRIDS:
    print(f'\n=== grid {grid}   (Ny = faces = m+1 for mimetic; collocated points for hichan)')
    print(f'{"operator":14s} {"G error Ny=33..257":^40s} ord | {"Laplacian error Ny=33..257":^40s} ord |'
          f' {"B dev":>7s} {"depth":>5s} {"w>0":>4s} | {"Re<0":>4s} {"rho*dy^2":>8s} {"|Im|/|lam|":>10s}')
    for op in [f'mimetic k={k}' for k in (4, 6, 8)] + ['hichan D2-4', 'hichan D2-6']:
        rr = [r for r in rows if r['op'] == op and r['grid'] == grid]
        eg = [r['eg'] for r in rr]; el = [r['el'] for r in rr]; r129 = [r for r in rr if r['ny'] == 129][0]
        g = f'{r129["bdev"]:7.2f} {r129["bdepth"]:5d} {"yes" if r129["wpos"] else "no":>4s}' if 'bdev' in r129 else f'{"-":>7s} {"-":>5s} {"-":>4s}'
        print(f'{op:14s} ' + ' '.join(f'{e:9.2e}' for e in eg) + f' {np.log2(eg[-2]/eg[-1]):4.1f} | '
              + ' '.join(f'{e:9.2e}' for e in el) + f' {np.log2(el[-2]/el[-1]):4.1f} | {g} |'
              f' {"yes" if r129["neg"] else "NO":>4s} {r129["rho"]:8.2f} {r129["im"]:10.1e}')
print('\nG error: max |G f - f\'| at the faces. Laplacian error: max |L f - f\'\'| over the interior (mimetic L = D G with the'
      '\nwall values of f; hichan: D2). B dev, depth: largest deviation of MOLE\'s boundary operator B = Q D + G^T P from'
      '\nthe plain boundary term, and how many rows from each wall B reaches; w>0: MOLE\'s weights P, Q positive (Ny = 129).'
      '\nRe<0, rho*dy^2, |Im|/|lam|: spectrum of the Dirichlet Laplacian at Ny = 129: all eigenvalues in the left half'
      '\nplane, spectral radius times dy_min^2 (explicit time-step limit), and the largest imaginary part relative to it'
      '\n(0 for a self-adjoint operator).')
