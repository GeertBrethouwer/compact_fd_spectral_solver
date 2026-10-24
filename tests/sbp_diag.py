"""
Diagonal-norm summation-by-parts first-derivative operators D = H^{-1} Q (Strand 1994; Mattsson & Nordstrom 2004)
on a uniform grid of unit spacing: central interior stencil of order 2s, r boundary rows exact for polynomials
of degree <= s, H = diag(h_1..h_r, 1, 1, ...), Q + Q^T = diag(-1, 0, ..., 0, 1).

The accuracy conditions Q x^m = H (x^m)' are linear in (h, Q), so a closure follows from one linear solve.
Free parameters left after that are chosen to minimise the boundary truncation error (degree s+1) and the
spectral radius, with all h_i > 0.  For s = 2 and 3 this reproduces Strand's norms (17/48, 59/48, 43/48, 49/48
and 13649/43200, 12013/8640, ...).
"""
import numpy as np, scipy.optimize as so, scipy.linalg as sl

CENTRAL = {2: [1/2], 4: [2/3, -1/12], 6: [3/4, -3/20, 1/60], 8: [4/5, -1/5, 4/105, -1/280]}

def free_entries(r, s):
    """upper-triangle Q entries (i<j) that are unknown: i < r and j < r+s"""
    return [(i, j) for i in range(r) for j in range(i+1, r+s)]

def assemble(x, r, s, n):
    """h (length n) and Q (n x n) for the left boundary; x = (h_1..h_r, free Q entries)"""
    c = CENTRAL[2*s]
    h = np.ones(n); h[:r] = x[:r]
    Q = np.zeros((n, n))
    for i in range(n):
        for k, ck in enumerate(c, 1):
            if i+k < n and not (i < r and i+k < r+s): Q[i, i+k] = ck; Q[i+k, i] = -ck
    for (i, j), v in zip(free_entries(r, s), x[r:]):
        Q[i, j] = v; Q[j, i] = -v
    Q[0, 0] = -1/2
    return h, Q

def residual(x, r, s):
    """boundary rows exact to degree s, rows r..r+s-1 (coupled to the closure) exact to degree 2s"""
    n = r+3*s+2; h, Q = assemble(x, r, s, n); L = n-1; xs = np.arange(n)/L   # scaled for conditioning
    res = lambda m: Q @ xs**m - h*m*xs**max(m-1, 0)/L
    return np.r_[np.concatenate([res(m)[:r] for m in range(s+1)]),
                 np.concatenate([res(m)[r:r+s] for m in range(2*s+1)])]

def closure(r, s, seeds=6):
    nu = r + len(free_entries(r, s))
    b0 = residual(np.zeros(nu), r, s)
    A = np.array([residual(e, r, s) - b0 for e in np.eye(nu)]).T
    x0 = np.linalg.lstsq(A, -b0, rcond=None)[0]
    if np.linalg.norm(A @ x0 + b0) > 1e-12*np.linalg.norm(b0):
        return None
    N = sl.null_space(A)
    if N.shape[1] == 0:
        return x0
    xi = np.arange(r+3*s+2, dtype=float)
    def obj(t):
        x = x0 + N @ t; h, Q = assemble(x, r, s, 2*(r+s)+4); D = Q / h[:, None]
        hs, Qs = assemble(x, r, s, r+3*s+2)
        err = (Qs @ xi**(s+1) - hs*(s+1)*xi**s)[:r]
        rho = max(abs(np.linalg.eigvals(D[1:, 1:])))
        return np.sum(err**2) + 0.05*rho**2 + 1e3*np.sum(np.minimum(h[:r]-0.05, 0)**2)
    best = min((so.minimize(obj, np.random.default_rng(k).normal(size=N.shape[1])*.1, method='Nelder-Mead',
                            options={'maxiter': 6000, 'xatol': 1e-12, 'fatol': 1e-15}) for k in range(seeds)),
               key=lambda o: o.fun)
    x = x0 + N @ best.x
    return x if (x[:r] > 0).all() else None

def sbp_xi(n, s, x, r):
    """diagonal norm H and D = H^{-1} Q on n points of unit spacing, both boundaries"""
    h, Q = assemble(x, r, s, n)
    J = np.eye(n)[::-1]
    hR, QR = J @ h, -(J @ Q @ J)
    m = r+s
    h[-r:] = hR[-r:]
    Q[-m:, -m:] = QR[-m:, -m:]
    return np.diag(h), Q / h[:, None]

ORDERS = {2: 4, 3: 6, 4: 8}            # s -> interior order
ROWS = {2: (4,), 3: (6,), 4: (8, 9, 10)}

def get(s):
    for r in ROWS[s]:
        x = closure(r, s)
        if x is not None:
            return x, r
    raise RuntimeError(f'no diagonal-norm closure for s = {s}')

if __name__ == '__main__':
    np.set_printoptions(precision=6, suppress=True, linewidth=160)
    for s in (2, 3, 4):
        x, r = get(s)
        n = 41; H, D = sbp_xi(n, s, x, r)
        E = np.diag(np.r_[-1., np.zeros(n-2), 1.]); xs = np.arange(n)/(n-1)
        acc = max(abs(D @ xs**m - m*xs**max(m-1, 0)/(n-1)).max() for m in range(s+1))
        print(f'SBP {s}-{2*s} (r = {r}): h = {np.diag(H)[:r]}\n  SBP defect {np.linalg.norm(H @ D + D.T @ H - E):.1e},'
              f' exactness to degree {s} on [0,1]: {acc:.1e}')
