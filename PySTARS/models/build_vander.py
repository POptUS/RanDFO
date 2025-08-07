import numpy as np

def build_vander(Y):
    m, n = Y.shape
    X = np.copy(Y)
    vander_terms = []

    for i in range(n):
        for j in range(i, n):
            term = 0.5 * Y[:, i] * Y[:, j] if i == j else Y[:, i] * Y[:, j]
            vander_terms.append(term.reshape(-1, 1))

    quadratic_part = np.hstack(vander_terms)
    return np.hstack((np.ones((m, 1)), X, quadratic_part))
