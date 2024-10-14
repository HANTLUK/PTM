import torch
import numpy as np
import numba as nb
import test_matrices
import sys

np.set_printoptions(suppress=True, linewidth=sys.maxsize)


# short-circuiting replacement for np.any()
@nb.jit(nopython=True)
def sc_any(array):
    for x in array.flat:
        if x:
            return True
    return False


def tpd(vector: np.ndarray, num_of_qubits: int) -> np.ndarray:
    if num_of_qubits == 0:
        return np.array(vector[0])
    else:
        matrix_dim: int = 1 << num_of_qubits
        matrix: np.ndarray = vector.reshape(matrix_dim, matrix_dim)
        return (0.5 ** num_of_qubits) * recursive_tpd(matrix, num_of_qubits, matrix_dim)


def recursive_tpd(matrix: np.ndarray, num_of_qubits: int, dimension: int) -> np.ndarray:
    if num_of_qubits == 0:
        return np.array(matrix[0][0])
    else:
        pauli_weights: np.ndarray = np.zeros(1 << (2 * num_of_qubits), dtype=np.complex64)
        num_of_qubits -= 1
        dimension >>= 1
        top_left: np.ndarray = matrix[0:dimension, 0:dimension]
        top_right: np.ndarray = matrix[0:dimension, dimension:]
        bottom_left: np.ndarray = matrix[dimension:, 0:dimension]
        bottom_right: np.ndarray = matrix[dimension:, dimension:]
        cmw_1: np.ndarray = top_left + bottom_right
        cmw_x: np.ndarray = bottom_left + top_right
        cmw_y: np.ndarray = -1.0j * (bottom_left - top_right)
        cmw_z: np.ndarray = top_left - bottom_right
        cmws: list[np.ndarray] = [cmw_1, cmw_x, cmw_y, cmw_z]

        index: int = 0
        for cmw in cmws:
            if sc_any(cmw):
                pauli_weights[index:index + dimension ** 2] = recursive_tpd(cmw, num_of_qubits, dimension)
            index += dimension ** 2

        return pauli_weights


def itpd(vector: np.ndarray, num_of_qubits: int) -> np.ndarray:
    if num_of_qubits == 0:
        return np.array(vector[0])
    else:
        matrix_dim: int = 1 << num_of_qubits
        matrix: np.ndarray = vector.reshape(matrix_dim, matrix_dim)
        return recursive_itpd(matrix, num_of_qubits, matrix_dim)


def recursive_itpd(matrix: np.ndarray, num_of_qubits: int, dimension: int) -> np.ndarray:
    if num_of_qubits == 0:
        return np.array(matrix[0][0])
    else:
        can_weights: np.ndarray = np.zeros(1 << (2 * num_of_qubits), dtype=np.complex64)
        num_of_qubits -= 1
        dimension >>= 1
        top_left: np.ndarray = matrix[0:dimension, 0:dimension]
        top_right: np.ndarray = matrix[0:dimension, dimension:]
        bottom_left: np.ndarray = matrix[dimension:, 0:dimension]
        bottom_right: np.ndarray = matrix[dimension:, dimension:]
        cmw_1: np.ndarray = top_left + bottom_right
        cmw_2: np.ndarray = bottom_left + top_right
        cmw_3: np.ndarray = 1.0j * (bottom_left - top_right)
        cmw_4: np.ndarray = top_left - bottom_right
        cmws: list[np.ndarray] = [cmw_1, cmw_2, cmw_3, cmw_4]

        index: int = 0
        for cmw in cmws:
            if sc_any(cmw):
                can_weights[index:index + dimension ** 2] = recursive_itpd(cmw, num_of_qubits, dimension)
            index += dimension ** 2

        return can_weights


if __name__ == "__main__":
    NUM_OF_QUBITS = 2
    mat = test_matrices.rand_diag_mat(2 ** NUM_OF_QUBITS)
    res = itpd(mat, NUM_OF_QUBITS)
    print(res, res.shape)
