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


# def iterative_tpd(vector: np.ndarray, num_of_qubits: int):
#     if num_of_qubits == 0:
#         return np.array(vector[0])
#     matrix_dim: int = 1 << num_of_qubits
#     cmws = [vector.reshape(matrix_dim, matrix_dim)]
#     for _ in range(num_of_qubits):
#         matrix_dim /= 2
#         num_of_qubits -= 1


def tpd(vector: np.ndarray, num_of_qubits: int) -> np.ndarray:
    if num_of_qubits == 0:
        return np.array(vector[0])
    else:
        matrix_dim: int = 1 << num_of_qubits
        matrix: np.ndarray = vector.reshape(matrix_dim, matrix_dim)
        return (0.5 ** num_of_qubits) * recursive_tpd(matrix, num_of_qubits, matrix_dim)


def tpd_new(vector: np.ndarray, num_of_qubits: int) -> np.ndarray:
    if num_of_qubits == 0:
        return np.array(vector[0])
    else:
        matrix_dim: int = 1 << num_of_qubits
        matrix: torch.Tensor = torch.from_numpy(vector.reshape(matrix_dim, matrix_dim))
        return (0.5 ** num_of_qubits) * recursive_tpd_new(matrix, num_of_qubits, matrix_dim)


def recursive_tpd_new(matrix: torch.Tensor, num_of_qubits: int, dimension: int) -> np.ndarray:
    if num_of_qubits == 0:
        return np.array(matrix[0][0])
    else:
        pauli_weights: np.ndarray = np.zeros(1 << (2 * num_of_qubits), dtype=np.complex64)
        num_of_qubits -= 1
        dimension >>= 1
        top_left: torch.Tensor = matrix[0:dimension, 0:dimension]
        top_right: torch.Tensor = matrix[0:dimension, dimension:]
        bottom_left: torch.Tensor = matrix[dimension:, 0:dimension]
        bottom_right: torch.Tensor = matrix[dimension:, dimension:]
        cmw_1: torch.Tensor = top_left + bottom_right
        cmw_x: torch.Tensor = bottom_left + top_right
        cmw_y: torch.Tensor = -1.0j * (bottom_left - top_right)
        cmw_z: torch.Tensor = top_left - bottom_right
        cmws: list[torch.Tensor] = [cmw_1, cmw_x, cmw_y, cmw_z]

        index: int = 0
        for cmw in cmws:
            if torch.any(cmw):
                pauli_weights[index:index + dimension ** 2] = recursive_tpd_new(cmw, num_of_qubits, dimension)
            index += dimension ** 2

        return pauli_weights


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


def itpd_new(vector: np.ndarray, num_of_qubits: int) -> np.ndarray:
    if num_of_qubits == 0:
        return np.array(vector[0])
    else:
        matrix_dim: int = 1 << num_of_qubits
        matrix: torch.Tensor = torch.from_numpy(vector.reshape(matrix_dim, matrix_dim))
        return recursive_itpd_new(matrix, num_of_qubits, matrix_dim)


def recursive_itpd_new(matrix: torch.Tensor, num_of_qubits: int, dimension: int) -> np.ndarray:
    if num_of_qubits == 0:
        return np.array(matrix[0][0])
    else:
        pauli_weights: np.ndarray = np.zeros(1 << (2 * num_of_qubits), dtype=np.complex64)
        num_of_qubits -= 1
        dimension >>= 1
        top_left: torch.Tensor = matrix[0:dimension, 0:dimension]
        top_right: torch.Tensor = matrix[0:dimension, dimension:]
        bottom_left: torch.Tensor = matrix[dimension:, 0:dimension]
        bottom_right: torch.Tensor = matrix[dimension:, dimension:]
        cmw_1: torch.Tensor = top_left + bottom_right
        cmw_2: torch.Tensor = bottom_left + top_right
        cmw_3: torch.Tensor = 1.0j * (bottom_left - top_right)
        cmw_4: torch.Tensor = top_left - bottom_right
        cmws: list[torch.Tensor] = [cmw_1, cmw_2, cmw_3, cmw_4]

        index: int = 0
        for cmw in cmws:
            if torch.any(cmw):
                pauli_weights[index:index + dimension ** 2] = recursive_itpd_new(cmw, num_of_qubits, dimension)
            index += dimension ** 2

        return pauli_weights


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
    num_of_qubits = 2
    mat = test_matrices.rand_diag_mat(2 ** num_of_qubits)
    res = iTPD(mat, num_of_qubits)
    print(res, res.shape)
