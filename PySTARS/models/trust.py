import numpy as np
from scipy.linalg import eigh

def seceqn(lambda_, eigval, alpha, delta):
    lambda_ = np.atleast_1d(lambda_)
    value = []
    for lam in lambda_:
        denom = eigval + lam
        with np.errstate(divide='ignore', invalid='ignore'):
            term = np.where(denom != 0, (alpha / denom) ** 2, np.inf)
        s_norm = np.sqrt(np.sum(term))
        if np.isnan(s_norm): s_norm = 0
        ## This is to ignore the divisible by 0 warning 
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            value.append(1.0 / delta - 1.0 / s_norm)

        # value.append(1.0 / delta - 1.0 / s_norm)
    return np.array(value)

def rfzero(Fun, x, itbnd, eigval, alpha, delta, tol=1e-12):
    dx = abs(x) / 2 if x != 0 else 0.5
    a = x
    fa = Fun(np.array([a]), eigval, alpha, delta)[0]
    b = x + dx
    fb = Fun(np.array([b]), eigval, alpha, delta)[0]

    count = 2
    while np.sign(fa) == np.sign(fb):
        dx *= 2
        b = x + dx
        fb = Fun(np.array([b]), eigval, alpha, delta)[0]
        count += 1
        if count > itbnd:
            break

    c = b
    fc = fb

    for _ in range(itbnd):
        if fb == 0 or abs(b - c) < tol:
            break
        if np.sign(fb) == np.sign(fc):
            c, fc = a, fa

        m = 0.5 * (c - b)
        toler = 2.0 * tol * max(abs(b), 1.0)
        if abs(m) <= toler:
            break

        d = m
        a = b
        fa = fb
        b += d if abs(d) > toler else (toler if b <= c else -toler)
        fb = Fun(np.array([b]), eigval, alpha, delta)[0]
        count += 1

    return b, c, count

def trust(g, H, delta):
    tol = 1e-12
    tol2 = 1e-8
    n = len(g)
    H = np.array(H)
    eigval, V = eigh(H)
    mineig = np.min(eigval)
    jmin = np.argmin(eigval)
    alpha = -V.T @ g
    sig = np.sign(alpha[jmin]) if alpha[jmin] != 0 else 1

    if mineig > 0:
        coeff = alpha / eigval
        s = V @ coeff
        nrms = np.linalg.norm(s)
        if nrms <= 1.2 * delta:
            val = g.T @ s + 0.5 * s.T @ H @ s
            return s, val, 1, 0, 0
        laminit = 0
    else:
        laminit = -mineig

    if seceqn(np.array([laminit]), eigval, alpha, delta)[0] > 0:
        b, _, count = rfzero(seceqn, laminit, 50, eigval, alpha, delta, tol)
        lambda_ = b
    else:
        lambda_ = -mineig
        count = 0

    w = eigval + lambda_
    coeff = np.where(w != 0, alpha / w, 0)
    s = V @ coeff
    nrms = np.linalg.norm(s)

    if nrms < 0.8 * delta:
        beta = np.sqrt(delta**2 - nrms**2)
        s += beta * sig * V[:, jmin]

    val = g.T @ s + 0.5 * s.T @ H @ s
    posdef = 1 if mineig > 0 else 0
    return s, val, posdef, count, lambda_
