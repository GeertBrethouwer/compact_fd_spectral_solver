"""
Compare first-derivative operators on uniform and hyperbolic-tangent grids (alpha_s = 2, 3, 4):

  D1     hichan's compact 4th-order scheme with coefficients computed on the stretched grid
  D1m    hichan's compact 4th-order scheme on the uniform xi grid, D1m = diag(1/y_xi) D_xi
  CSBP   compact 4th-order (Pade) interior with a summation-by-parts closure, norm = LHS matrix P:
         P f' = Q f, P symmetric positive definite, Q + Q^T = diag(-1,0,...,0,1); mapped like D1m
  SBP24  explicit diagonal-norm SBP operator (4th-order interior, 2nd-order boundary), mapped like D1m

Usage: python compare_sbp.py OPSDIR   (OPSDIR holds the ops_<grid>_<ny>.bin files from dump_ops.f90)
"""
import sys, numpy as np, scipy.optimize as so, scipy.linalg as sl
ops = sys.argv[1]
GRIDS = ['uniform', 'tanh2', 'tanh3', 'tanh4']
NYS = [33, 65, 129, 257, 513]

def load(grid, ny):
    a = np.fromfile(f'{ops}/ops_{grid}_{ny}.bin', dtype='<f8', offset=4)
    y, wq = a[:ny], a[ny:2*ny]
    d1 = a[2*ny:2*ny+ny*ny].reshape(ny, ny).T          # Fortran column order
    d1m = a[2*ny+ny*ny:].reshape(ny, ny).T
    return y, wq, d1, d1m

# ---------------------------------------------------------------- CSBP closure (unit spacing)
R, BP, BQ = 4, 1, 2                                     # boundary rows, half-bandwidths of P and Q
def csbp_block(x):
    n = R+4; P = np.zeros((n, n)); Q = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i >= R or j >= R:
                P[i, j] = {0: 2/3, 1: 1/6}.get(abs(i-j), 0.)
                Q[i, j] = {1: .5, -1: -.5}.get(j-i, 0.)
    k = 0
    for i in range(R):
        for j in range(i, R):
            if j-i <= BP: P[i, j] = P[j, i] = x[k]; k += 1
    for i in range(R):
        for j in range(i+1, R):
            if j-i <= BQ: Q[i, j] = x[k]; Q[j, i] = -x[k]; k += 1
    Q[0, 0] = -.5
    return P, Q, k
def csbp_res(x, degs):
    P, Q, _ = csbp_block(x); xi = np.arange(P.shape[0], dtype=float)
    return np.concatenate([(P @ (m*xi**max(m-1, 0)) - Q @ xi**m)[:R] for m in degs])
NU = csbp_block(np.zeros(99))[2]
A = np.array([csbp_res(e, (0, 1)) - csbp_res(0*e, (0, 1)) for e in np.eye(NU)]).T
x0 = np.linalg.lstsq(A, -csbp_res(np.zeros(NU), (0, 1)), rcond=None)[0]
N = sl.null_space(A)
def csbp_obj(t):                                        # exact to degree 1; smallest degree-2..4 error, P SPD
    x = x0 + N @ t; P = csbp_block(x)[0]
    pen = max(0., 0.2 - np.linalg.eigvalsh(P[:R+2, :R+2]).min())
    return np.sum(csbp_res(x, (2, 3, 4))**2 * np.repeat([1, .1, .01], R)) + 1e3*pen**2
t = min((so.minimize(csbp_obj, np.random.default_rng(s).normal(size=N.shape[1]), method='Nelder-Mead',
                     options={'maxiter': 40000, 'xatol': 1e-13, 'fatol': 1e-16}) for s in range(10)),
        key=lambda o: o.fun).x
PB, QB, _ = csbp_block(x0 + N @ t)

def csbp_xi(n):
    """P, Q on n points of unit spacing (both boundaries)"""
    P = np.diag(np.full(n, 2/3)) + np.diag(np.full(n-1, 1/6), 1) + np.diag(np.full(n-1, 1/6), -1)
    Q = np.diag(np.full(n-1, .5), 1) - np.diag(np.full(n-1, .5), -1)
    P[:R+1, :R+1] = PB[:R+1, :R+1]; Q[:R+1, :R+1] = QB[:R+1, :R+1]
    J = np.eye(n)[::-1]                                  # mirror for the upper boundary
    P[-R-1:, -R-1:] = (J @ P @ J)[-R-1:, -R-1:]; Q[-R-1:, -R-1:] = -(J @ Q @ J)[-R-1:, -R-1:]
    return P, Q

# ---------------------------------------------------------------- explicit diagonal-norm SBP 2-4
def sbp24_xi(n):
    H = np.ones(n); H[:4] = H[-4:][::-1] = [17/48, 59/48, 43/48, 49/48]
    D = np.zeros((n, n))
    for i in range(n):
        if 2 <= i <= n-3: D[i, i-2:i+3] = [1/12, -2/3, 0, 2/3, -1/12]
    B = [[-24/17, 59/34, -4/17, -3/34, 0, 0], [-1/2, 0, 1/2, 0, 0, 0],
         [4/43, -59/86, 0, 59/86, -4/43, 0], [3/98, 0, -59/98, 0, 32/49, -4/49]]
    for i in range(4):
        D[i, :6] = B[i]; D[n-1-i, n-6:] = -np.array(B[i])[::-1]
    return np.diag(H), D

def mapped(Dxi, h, y):
    """D_y = diag(1/y_xi) D_xi / h with y_xi from the same operator"""
    yx = Dxi @ y / h
    assert (yx > 0).all()
    return Dxi / h / yx[:, None], yx

E = lambda n: np.diag(np.r_[-1., np.zeros(n-2), 1.])
def skew_defect(D, H):
    return np.linalg.norm(H @ D + D.T @ H - E(len(D)), 2) / np.linalg.norm(H @ D, 2)
def best_diag_defect(D):
    """smallest ||diag(w) D + D^T diag(w) - E||_2 / ||diag(w) D||_2 over diagonal norms (LS fit)"""
    n = len(D); M = np.array([(np.outer(np.eye(n)[i], D[i]) + np.outer(D[i], np.eye(n)[i])).ravel() for i in range(n)]).T
    w = np.linalg.lstsq(M, E(n).ravel(), rcond=None)[0]
    return skew_defect(D, np.diag(w)), w.min() > 0

def stability(D, y):
    """u_t + u_y = 0, inflow u(-1)=0 (rows/cols 1 removed) and u_t - u_y = 0, inflow u(1)=0;
    max Re(lambda) and spectral radius, both times dy_min"""
    dy = np.diff(y).min()
    l1 = np.linalg.eigvals(-D[1:, 1:]); l2 = np.linalg.eigvals(D[:-1, :-1])
    return max(l1.real.max(), l2.real.max())*dy, max(abs(l1).max(), abs(l2).max())*dy

def growth(D, y, T=4.):
    """max over t in (0,T] of ||exp(t L)||_2 for u_t + u_y = 0 with inflow at y=-1 (transient growth)"""
    L = -D[1:, 1:]; return max(np.linalg.norm(sl.expm(t*L), 2) for t in np.linspace(.05, T, 40))

def skew_advection(D, y):
    """u_t = -1/2 [v D u + D (v u)], v = (1-y^2)(1/2+y) (zero at the walls, changes sign):
    energy-conserving in the continuous case; max Re(lambda)*dy_min and max_t ||exp(tL)|| for t <= 20"""
    v = (1-y**2)*(.5+y); L = -.5*(v[:, None]*D + D*v[None, :])
    dy = np.diff(y).min()
    return np.linalg.eigvals(L).real.max()*dy, max(np.linalg.norm(sl.expm(t*L), 2) for t in np.linspace(.5, 20, 40))

f = lambda y: np.sin(2.3*y+.4); df = lambda y: 2.3*np.cos(2.3*y+.4)
res = {}
for grid in GRIDS:
    for ny in NYS:
        y, wq, d1, d1m = load(grid, ny)
        h = 2/(ny-1); xi = np.linspace(-1, 1, ny)
        P, Q = csbp_xi(ny); Dc, yxc = mapped(np.linalg.solve(P, Q), h, y)
        Hs, Ds = sbp24_xi(ny); Dsb, yxs = mapped(Ds, h, y)
        ops_ = {'D1': d1, 'D1m': d1m, 'CSBP': Dc, 'SBP24': Dsb}
        own = {'D1': np.diag(wq), 'D1m': np.diag(wq), 'CSBP': h*np.diag(np.sqrt(yxc)) @ P @ np.diag(np.sqrt(yxc)),
               'SBP24': h*Hs*yxs[None, :]}
        for k, D in ops_.items():
            err = abs(D @ f(y) - df(y))
            r = dict(emax=err.max(), eint=err[abs(y) < .5].max())
            if ny <= 257:
                r['own'] = skew_defect(D, own[k]); r['diag'], r['wpos'] = best_diag_defect(D) if ny <= 129 else (np.nan, True)
                r['reλ'], r['ρ'] = stability(D, y)
            if ny == 129: r['G'] = growth(D, y); r['sk_reλ'], r['sk_G'] = skew_advection(D, y)
            res[grid, ny, k] = r

np.set_printoptions(precision=5, suppress=True, linewidth=160)
print('CSBP boundary block P (unit spacing, rows 1..5):\n', PB[:R+1, :R+2]); print('Q:\n', QB[:R+1, :R+2])
print('min eig P block:', np.linalg.eigvalsh(PB[:R+2, :R+2]).min())
for grid in GRIDS:
    print(f'\n=== grid {grid}')
    print(f'{"op":6s} {"max error Ny=33..513":^55s} {"order":>5s} | {"|y|<.5 err @257":>14s} {"order":>5s} |'
          f' {"skew own":>9s} {"best diag":>9s} {"maxRe*dy":>9s} {"rho*dy":>7s} {"G":>9s} | {"skw Re*dy":>9s} {"skw G":>9s}')
    for k in ['D1', 'D1m', 'CSBP', 'SBP24']:
        em = [res[grid, n, k]['emax'] for n in NYS]; ei = [res[grid, n, k]['eint'] for n in NYS]
        r = res[grid, 129, k]
        print(f'{k:6s} ' + ' '.join(f'{e:9.2e}' for e in em) + f' {np.log2(em[-2]/em[-1]):5.2f} | {ei[3]:14.2e}'
              f' {np.log2(ei[2]/ei[3]):5.2f} | {r["own"]:9.1e} {r["diag"]:9.1e}{"" if r["wpos"] else "*"} {r["reλ"]:9.1e} {r["ρ"]:7.3f} {r["G"]:9.2g} | {r["sk_reλ"]:9.1e} {r["sk_G"]:9.3g}')
print('\nskew own: ||H D + D^T H - E||/||H D|| in the operator\'s own norm (D1, D1m: quadrature wq; CSBP: h Y^1/2 P Y^1/2; '
      'SBP24: h H diag(y_xi)), Ny = 129.\nbest diag: same, minimised over all diagonal norms (* = some weight <= 0).'
      '\nmaxRe*dy, rho*dy: max real part and spectral radius of the advection operator with inflow BC, times dy_min (Ny = 129).'
      '\nG: max_t ||exp(t L)|| for t <= 4 (transient growth; 1 for a contraction), Ny = 129.'
      '\nskw Re*dy, skw G: the same for skew-symmetric advection -1/2[v D u + D(v u)] with v = (1-y^2)(1/2+y), t <= 20, Ny = 129'
      '\n(energy-conserving in the continuous case: max Re = 0 and G = 1 in the right norm).')
