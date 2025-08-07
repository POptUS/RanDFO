import numpy as np

def polynomial_to_matrix(L, n):
    b = L[0]
    c = L[1:n+1]
    L = L[n+1:]

    Q = np.zeros((n, n))
    for k in range(n):
        Q[k, k:] = L[:n - k]
        L = L[n - k:]

    D = np.diag(Q)
    Q = Q - np.diag(D)
    Q = Q + Q.T + np.diag(D)
    return Q, c, b
