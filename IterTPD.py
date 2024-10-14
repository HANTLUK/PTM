import numpy as np
import numba as nb
import test_matrices
import sys

np.set_printoptions(suppress=True, linewidth=sys.maxsize)


@nb.jit(nopython=True)
def sc_any(array):
    for x in array.flat:
        if x:
            return True
    return False


def iter_tpd(vector: np.ndarray, num_of_qubits: int) -> np.ndarray:
    if num_of_qubits == 0:
        return np.array(vector[0])
    else:
        matrix_dim: int = 1 << num_of_qubits
        vector_dim: int = 1 << (2 * num_of_qubits)
        matrix: np.ndarray = vector.reshape(matrix_dim, matrix_dim)
        pauli_weights: np.ndarray = np.zeros(vector_dim, dtype=np.complex64)
        recursive_tpd(matrix, num_of_qubits, matrix_dim, pauli_weights)
        return (0.5 ** num_of_qubits) * pauli_weights


def recursive_tpd(matrix: np.ndarray, num_of_qubits: int, dimension: int,
                  pauli_weights: np.ndarray, index_outer: int = 0):
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

    if num_of_qubits == 0:
        index: int = 0
        for cmw in cmws:
            if cmw.any():
                pauli_weights[index_outer + index] = cmw[0][0]
            index += 1
    else:
        index: int = 0
        for cmw in cmws:
            if sc_any(cmw):
                recursive_tpd(cmw, num_of_qubits, dimension, pauli_weights, index_outer + index)
            index += dimension ** 2

    return


def itpd(vector: np.ndarray, num_of_qubits: int):
    if num_of_qubits == 0:
        return np.array(vector[0])
    else:
        matrix_dim: int = 1 << num_of_qubits
        vector_dim: int = 1 << (2 * num_of_qubits)
        matrix: np.ndarray = vector.reshape(matrix_dim, matrix_dim)
        can_weights: np.ndarray = np.zeros(vector_dim, dtype=np.complex64)
        recursive_itpd(matrix, num_of_qubits, matrix_dim, can_weights)
        return can_weights

def recursive_itpd(matrix: np.ndarray, num_of_qubits: int, dimension: int,
                   can_weights: np.ndarray, index_outer: int = 0):

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

    if num_of_qubits == 0:
        index: int = 0
        for cmw in cmws:
            if cmw.any():
                can_weights[index_outer + index] = cmw[0][0]
            index += 1
    else:
        index: int = 0
        for cmw in cmws:
            if sc_any(cmw):
                recursive_itpd(cmw, num_of_qubits, dimension, can_weights, index_outer + index)
            index += dimension ** 2

    return


if __name__ == "__main__":
    NUM_OF_QUBITS = 2
    mat = test_matrices.rand_diag_mat(2 ** NUM_OF_QUBITS)
    print(mat, mat.shape)
    res = iter_itpd(mat, NUM_OF_QUBITS)
    print(res)
