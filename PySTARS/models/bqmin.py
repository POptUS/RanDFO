import numpy as np

def bqmin(A, B, L, U, tol=1e-13, maxit=5000):
    n = A.shape[0]
    X = np.zeros(n)  # Start from zero
    G = A @ X + B
    Projg = X - np.clip(X - G, L, U)
    it = 0

    while it < maxit and np.linalg.norm(Projg) > tol:
        it += 1
        t = 1
        pap = Projg @ A @ Projg
        if pap > 0:
            t = min(1, (Projg @ G) / pap)

        X = X - t * Projg
        G = A @ X + B
        Projg = X - np.clip(X - G, L, U)

    f = 0.5 * X @ A @ X + B @ X
    return X, f