import numpy as np
from models.build_vander import build_vander
from models.maximize_lagrange_polynomial import maximize_lagrange_polynomial
from models.shift_and_scale import shift_and_scale
from models.undo_shift_and_scale import undo_shift_and_scale

def algorithm_6_4(Y, Delta, f):
    x_0 = Y[0, :].reshape(-1, 1)
    Y, scalefactor = shift_and_scale(Y)

    p_ini, n = Y.shape
    p = int((n + 1) * (n + 2) / 2)

    X = build_vander(Y)
    U = np.eye(p)

    for i in range(p):
        ab_values = np.abs(X @ U[i, :].reshape(-1, 1)).flatten()
        ab_values[:i] = 0
        value = np.max(ab_values)
        j_i = np.argmax(ab_values)

        if j_i < i and i < p_ini:
            j_i = i
            print("Point selection error in algorithm_6_4")

        if value > 0 and i < p_ini:
            Y[[i, j_i]] = Y[[j_i, i]]
            X[[i, j_i]] = X[[j_i, i]]
            f[[i, j_i]] = f[[j_i, i]]
        else:
            argmax, _ = maximize_lagrange_polynomial(U[i, :], np.zeros_like(x_0), Delta / scalefactor)
            new_point = argmax.flatten().reshape(1, -1)
            Y = np.vstack([Y, new_point])
            f = np.vstack([f, [[np.inf]]])
            X = build_vander(Y)

        for j in range(i + 1, p):
            U[j, :] -= ((U[j, :] @ X[i, :]) / (U[i, :] @ X[i, :])) * U[i, :]

    Y = undo_shift_and_scale(Y, x_0.flatten(), scalefactor)
    return Y, f
