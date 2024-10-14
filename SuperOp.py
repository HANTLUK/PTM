import numpy as np
import numba as nb
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import test_matrices
from TPD import tpd, itpd
from qiskit.quantum_info.operators.channel import SuperOp, PTM
import sys

np.set_printoptions(suppress=True, linewidth=sys.maxsize)


def can_to_ptm(matrix):
    num_of_qubits: int = (matrix.shape[0].bit_length() - 1) >> 1
    if num_of_qubits < 6:
        intermediate = np.apply_along_axis(tpd, 1, matrix, num_of_qubits)
        return np.apply_along_axis(itpd, 0, intermediate, num_of_qubits)
    else:
        intermediate = apply_along_row_parallel(tpd, matrix, num_of_qubits)
        return apply_along_row_parallel(itpd, intermediate.T, num_of_qubits).T


def apply_along_row_parallel(func, data, num_of_qubits):
    with ProcessPoolExecutor() as executor:
        results = list(executor.map(func_with_args, data, [func]*len(data), [num_of_qubits]*len(data)))
    return np.array(results)


def func_with_args(row, func, num_of_qubits):
    return func(row, num_of_qubits)


if __name__ == "__main__":
    NUM_OF_QUBITS: int = 2
    DATA = test_matrices.rand_dense_mat(4 ** NUM_OF_QUBITS)
    print(np.array_str(DATA, precision=1))
    mat = can_to_ptm(DATA)

    SOp = SuperOp(DATA)
    mat2 = PTM(SOp).data
    dif = mat2 - mat
    print("dif\n", np.array_str(dif, precision=1))
