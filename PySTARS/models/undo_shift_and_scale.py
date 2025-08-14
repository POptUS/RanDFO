import numpy as np

def undo_shift_and_scale(Y_hat, x, maxnorm):
    return Y_hat * maxnorm + x
