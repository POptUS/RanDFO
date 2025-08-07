import numpy as np

def generate_random_Q(n, p):
    """Generating a random Gaussian matrix Q ∈ ℝⁿˣᵖ"""
    return np.random.normal(0, 1 / p, size=(n, p))

def haar_orthog_matrix(dimension):
    A = np.random.randn(dimension, dimension)
    Q, R = np.linalg.qr(A)
    Q = Q @ np.diag(np.sign(np.diag(R)))
    return Q

def generate_random_hashing_matrix(n, p, seed=None):
    rng = np.random.default_rng(seed)
    return rng.choice([-1, 1], size=(p, n)) / np.sqrt(p)

def generate_subspace_matrix(jlm_type, n, p, hash_param=None):
    if jlm_type == 2:
        assert n == p, "Identity matrix can only be used when n == p"
        return np.eye(p)
    elif jlm_type == 3:
        Q_full = haar_orthog_matrix(n)
        return np.sqrt(n / p) * Q_full[:, :p]
    elif jlm_type == 4:
        Q_full = haar_orthog_matrix(n)
        return Q_full[:, :p]
    else:
        return generate_random_hashing_matrix(n, p, seed=hash_param).T

