import numpy as np

def shift_and_scale(Y):
    Y_hat = Y - Y[0, :]
    if Y_hat.shape[0] > 1:
        maxnorm = max(np.linalg.norm(row) for row in Y_hat)
        Y_hat = Y_hat / maxnorm
    else:
        maxnorm = 1
    return Y_hat, maxnorm
