import numpy as np
import test_matrices
from TPD import tpd, itpd, TPD, iTPD
from qiskit.quantum_info.operators.channel import SuperOp, PTM
import sys

np.set_printoptions(suppress=True, linewidth=sys.maxsize, threshold=sys.maxsize)


def can_to_ptm_old(matrix):
    # matDim = matrix.shape[0]
    qBitDim = int(np.log(matrix.shape[0]) / np.log(4))
    W = np.apply_along_axis(TPD, 1, matrix, qBitDim)
    PTM = np.apply_along_axis(iTPD, 0, W, qBitDim)
    return PTM


def can_to_ptm(matrix):
    # matDim = matrix.shape[0]
    qBitDim = int(np.log(matrix.shape[0]) / np.log(4))
    W = np.apply_along_axis(tpd, 1, matrix, qBitDim)
    PTM = np.apply_along_axis(itpd, 0, W, qBitDim)
    return PTM


if __name__ == "__main__":
    qDim = 2
    data = test_matrices.rand_dense_mat(4 ** qDim)
    print(np.array_str(data, precision=1))
    mat = can_to_ptm(data)

    SOp = SuperOp(data)
    mat2 = PTM(SOp).data
    dif = mat2 - mat
    print("dif\n", np.array_str(dif, precision=1))
