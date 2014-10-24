import numpy as np

def rec_to_mat(rec):
    if rec.dtype.names is not None:
        rec = rec.astype([(k, np.float64) for k in rec.dtype.names])
        return np.mat(rec.view((np.float64, len(rec.dtype.names))))
    else:
        return np.matrix(rec.astype(np.float64))
