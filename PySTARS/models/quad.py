import numpy as np
from numpy.linalg import lstsq, qr

def phi2eval(X):
    X = np.atleast_2d(X)
    m, n = X.shape
    num_terms = int(n * (n + 1) / 2)
    Phi = np.zeros((m, num_terms))
    idx = 0
    for i in range(n):
        Phi[:, idx] = 0.5 * X[:, i] ** 2
        idx += 1
        for j in range(i + 1, n):
            Phi[:, idx] = X[:, i] * X[:, j] / np.sqrt(2)
            idx += 1
    return Phi

def fitfroquad(poised_set, f_poised):
    f_poised = f_poised.flatten()
    m, n = poised_set.shape

    if m <= n + 1:
        M = np.hstack((np.ones((m, 1)), poised_set))  # m x (n+1)
        coeffs, *_ = lstsq(M, f_poised, rcond=None)
        G = coeffs[1:]
        H = np.zeros((n, n))
        return G, H

    M = np.hstack((np.ones((m, 1)), poised_set))  # m x (n+1)
    Q, _ = qr(M, mode='reduced')                 # QR of M
    Z = np.eye(m) - Q @ Q.T                      # Projection onto null space

    Phi = phi2eval(poised_set)  # m x q
    N = Phi.T                   # q x m
    L = np.dot(N,Z)             # q x m           
    rhs = np.dot(Z,f_poised)    # m x 1   
    Beta, *_ = lstsq(L.T, rhs, rcond=None)  # q x 1

    residual = f_poised - Phi @ Beta
    Alpha, *_ = lstsq(M, residual, rcond=None)
    G = Alpha[1:]

    H = np.zeros((n, n))
    idx = 0
    for i in range(n):
        H[i, i] = Beta[idx]
        idx += 1
        for j in range(i + 1, n):
            H[i, j] = Beta[idx] / np.sqrt(2)
            H[j, i] = H[i, j]
            idx += 1

    return G, H
