import numpy as np
from models.build_vander import build_vander
from models.polynomial_to_matrix import polynomial_to_matrix
from models.trust import trust

def maximize_lagrange_polynomial(L, x_0, Delta):
    max_val = 0
    Q, c, _ = polynomial_to_matrix(L, x_0.shape[0])

    for j in [0, 1]:
        linear_term = Q @ x_0 + c.reshape(-1, 1)
        Q = (-1)**j * Q
        linear_term = (-1)**j * linear_term

        if not np.any(Q):
            test_step = -linear_term / np.linalg.norm(linear_term)
            test_step = test_step * Delta
        else:
            test_step, *_ = trust(linear_term.flatten(), Q, Delta)

        test_x = x_0 + test_step.reshape(-1,1)
        test_val = np.abs(L @ build_vander(test_x.T).T)[0]

        if test_val > max_val:
            max_val = test_val
            argmax = test_x

    return argmax, max_val
