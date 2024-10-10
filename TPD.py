import numpy as np
import math

import test_matrices as TM

import sys

np.set_printoptions(suppress=True, linewidth=sys.maxsize, threshold=sys.maxsize)


def tpd(vector: np.ndarray, num_of_qubits: int):
    if num_of_qubits == 0:
        return np.array(vector[0])
    else:
        matrix_dim: int = 1 << num_of_qubits
        matrix: np.ndarray = vector.reshape(matrix_dim, matrix_dim)
        return recursive_tpd(matrix, num_of_qubits)


def recursive_tpd(matrix: np.ndarray, num_of_qubits: int):
    if num_of_qubits == 0:
        return np.array(matrix[0][0])
    else:
        pauli_weights: np.ndarray = np.zeros(1 << (2 * num_of_qubits), dtype=np.complex64)
        num_of_qubits -= 1
        halved_dim: int = 1 << num_of_qubits
        top_left: np.ndarray = matrix[0:halved_dim, 0:halved_dim]
        top_right: np.ndarray = matrix[0:halved_dim, halved_dim:]
        bottom_left: np.ndarray = matrix[halved_dim:, 0:halved_dim]
        bottom_right: np.ndarray = matrix[halved_dim:, halved_dim:]
        cmw_1: np.ndarray = 0.5 * (top_left + bottom_right)
        cmw_x: np.ndarray = 0.5 * (bottom_left + top_right)
        cmw_y: np.ndarray = -0.5j * (bottom_left - top_right)
        cmw_z: np.ndarray = 0.5 * (top_left - bottom_right)
        cmws: list[np.ndarray] = [cmw_1, cmw_x, cmw_y, cmw_z]

        index: int = 0
        for cmw in cmws:
            if cmw.any():
                pauli_weights[index:index + halved_dim ** 2] = recursive_tpd(cmw, num_of_qubits)
            index += halved_dim ** 2

        return pauli_weights


def itpd(vector: np.ndarray, num_of_qubits: int):
    if num_of_qubits == 0:
        return np.array(vector[0])
    else:
        matrix_dim: int = 1 << num_of_qubits
        matrix: np.ndarray = vector.reshape(matrix_dim, matrix_dim)
        return recursive_itpd(matrix, num_of_qubits)


def recursive_itpd(matrix: np.ndarray, num_of_qubits: int):
    if num_of_qubits == 0:
        return np.array(matrix[0][0])
    else:
        can_weights: np.ndarray = np.zeros(1 << (2 * num_of_qubits), dtype=np.complex64)
        num_of_qubits -= 1
        halved_dim: int = 1 << num_of_qubits
        top_left: np.ndarray = matrix[0:halved_dim, 0:halved_dim]
        top_right: np.ndarray = matrix[0:halved_dim, halved_dim:]
        bottom_left: np.ndarray = matrix[halved_dim:, 0:halved_dim]
        bottom_right: np.ndarray = matrix[halved_dim:, halved_dim:]
        cmw_1: np.ndarray = top_left + bottom_right
        cmw_2: np.ndarray = bottom_left + top_right
        cmw_3: np.ndarray = 1.0j * (bottom_left - top_right)
        cmw_4: np.ndarray = top_left - bottom_right
        cmws: list[np.ndarray] = [cmw_1, cmw_2, cmw_3, cmw_4]

        index: int = 0
        for cmw in cmws:
            if cmw.any():
                can_weights[index:index + halved_dim ** 2] = recursive_itpd(cmw, num_of_qubits)
            index += halved_dim ** 2

        return can_weights


def TPD(matrix, qBitDim, flag=True):
    """
        Computes the Pauli decomposition of a square matrix.

        Iteratively splits tensor factors off and decomposes those smaller
        matrices. This is done using submatrices of the original matrix.
        The Pauli strings are generated in each step.

        Args:
            matrix: Matrix to be decomposed
    """
    # Array in row column convention
    if len(matrix.shape) == 1:
        matDim = 1 << qBitDim
        matrix = matrix.reshape(matDim, matDim)

    decomposition = []
    # Output for dimension 1

    if qBitDim == 0:
        decomposition = [matrix[0, 0]]
        del matrix

    # Calculates the tensor product coefficients via the sliced submatrices.
    # If one of these components is zero that coefficient is ignored.

    if qBitDim > 0:
        qBitDim -= 1
        halfDim = 1 << qBitDim

        coeff1 = 0.5 * (matrix[0:halfDim, 0:halfDim]
                        + matrix[halfDim:, halfDim:])
        coeffX = 0.5 * (matrix[halfDim:, 0:halfDim]
                        + matrix[0:halfDim, halfDim:])
        coeffY = -1.j * 0.5 * (matrix[halfDim:, 0:halfDim]
                               - matrix[0:halfDim, halfDim:])
        coeffZ = 0.5 * (matrix[0:halfDim, 0:halfDim]
                        - matrix[halfDim:, halfDim:])

        coefficients = [coeff1, coeffX, coeffY, coeffZ]
        del matrix

        # Recursion for the Submatrices

        for coeff in coefficients:
            mat = coeff
            if mat.any() != 0:
                subDec = TPD(mat, qBitDim, flag)
            else:
                subDec = np.zeros(halfDim ** 2)
            decomposition.extend(subDec)
    if flag:
        return np.array(decomposition)
    return decomposition


def iTPD(matrix, qBitDim, flag=None):
    """
        Computes the inverse Pauli decomposition of a square matrix.

        Iteratively splits tensor factors off and decomposes those smaller
        matrices. This is done using submatrices of the original matrix.
        The Pauli strings are generated in each step.

        Args:
            matrix: Matrix to be decomposed
    """

    # Array in row column convention
    if len(matrix.shape) == 1:
        matDim = 1 << qBitDim
        matrix = matrix.reshape(matDim, matDim)

    if flag is None:
        flag = True

    decomposition = []
    # Output for dimension 1

    if qBitDim == 0:
        decomposition = [matrix[0, 0]]
        del matrix

    # Calculates the tensor product coefficients via the sliced submatrices.
    # If one of these components is zero that coefficient is ignored.

    if qBitDim > 0:
        qBitDim -= 1
        halfDim = 1 << qBitDim

        coeff1 = (matrix[0:halfDim, 0:halfDim]
                  + matrix[halfDim:, halfDim:])
        coeff2 = (matrix[halfDim:, 0:halfDim]
                  + matrix[0:halfDim, halfDim:])
        coeff3 = (1.j * matrix[halfDim:, 0:halfDim]
                  - 1.j * matrix[0:halfDim, halfDim:])
        coeff4 = (matrix[0:halfDim, 0:halfDim]
                  - matrix[halfDim:, halfDim:])

        coefficients = [coeff1, coeff2, coeff3, coeff4]
        del matrix

        # Recursion for the Submatrices

        for coeff in coefficients:
            mat = coeff
            if mat.any() != 0:
                subDec = iTPD(mat, qBitDim, flag)
            else:
                subDec = np.zeros(halfDim ** 2)
            decomposition.extend(subDec)
    if flag:
        return np.array(decomposition)
    return decomposition


if __name__ == "__main__":
    qDim = 2
    mat = TM.rand_diag_mat(2 ** qDim)
    res = iTPD(mat, qDim)
    print(res, res.shape)
