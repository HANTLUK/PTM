import numpy as np
import test_matrices
from TPD import tpd, itpd, tpd_new, itpd_new
# from IterTPD import tpd, itpd
from qiskit.quantum_info.operators.channel import SuperOp, PTM
import sys

np.set_printoptions(suppress=True, linewidth=sys.maxsize)


def can_to_ptm_new(matrix):
    num_of_qubits: int = int(np.log(matrix.shape[0]) / np.log(4))
    W = np.apply_along_axis(tpd_new, 1, matrix, num_of_qubits)
    PTM = np.apply_along_axis(itpd_new, 0, W, num_of_qubits)
    return PTM


def can_to_ptm(matrix):
    num_of_qubits: int = int(np.log(matrix.shape[0]) / np.log(4))
    W = np.apply_along_axis(tpd, 1, matrix, num_of_qubits)
    PTM = np.apply_along_axis(itpd, 0, W, num_of_qubits)
    return PTM


if __name__ == "__main__":
    qDim = 2
    data = test_matrices.rand_dense_mat(4 ** qDim)
    print(np.array_str(data, precision=1))
    mat = can_to_ptm_new(data)

    SOp = SuperOp(data)
    mat2 = PTM(SOp).data
    dif = mat2 - mat
    print("dif\n", np.array_str(dif, precision=1))
