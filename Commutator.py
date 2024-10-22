import numpy as np
import test_matrices
import sys

import logging
logger = logging.getLogger(__name__)

np.set_printoptions(suppress=True, linewidth=sys.maxsize, threshold=sys.maxsize)

comX: np.ndarray = np.array([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, -1j], [0, 0, 1j, 0]])
comY: np.ndarray = np.array([[0, 0, 0, 0], [0, 0, 0, 1j], [0, 0, 0, 0], [0, -1j, 0, 0]])
comZ: np.ndarray = np.array([[0, 0, 0, 0], [0, 0, -1j, 0], [0, 1j, 0, 0], [0, 0, 0, 0]])

SINGLE_COMMUTATOR: dict = {"X": comX, "Y": comY, "Z": comZ}

antiI: np.ndarray = np.eye(4)
antiX: np.ndarray = np.array([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
antiY: np.ndarray = np.array([[0, 0, 1, 0], [0, 0, 0, 0], [1, 0, 0, 0], [0, 0, 0, 0]])
antiZ: np.ndarray = np.array([[0, 0, 0, 1], [0, 0, 0, 0], [0, 0, 0, 0], [1, 0, 0, 0]])

SINGLE_ANTICOMMUTATOR: dict = {"I": antiI, "X": antiX, "Y": antiY, "Z": antiZ}

pauliI: list = [1., np.add, [0, 3]]
pauliX: list = [1., np.add, [1, 2]]
pauliY: list = [1.j, np.subtract, [1, 2]]
pauliZ: list = [1., np.subtract, [0, 3]]
paulis: dict = {"I": pauliI, "X": pauliX, "Y": pauliY, "Z": pauliZ}

PAULI_LIST: list = ["I", "X", "Y", "Z"]


def commutator_to_ptm(matrix: np.ndarray):
    matrix_dim: int = matrix.shape[0]
    num_of_qubits: int = matrix_dim.bit_length() - 1
    PTM_dim: int = 1 << (2 * num_of_qubits)
    PTM = recursive_commutator(matrix, num_of_qubits)
    return PTM

def recursive_commutator(matrix: np.ndarray, num_of_qubits: int):
    PTM_dim: int = 1 << (2*num_of_qubits)
    PTM: np.ndarray = np.zeros((PTM_dim, PTM_dim), dtype=np.complex64)
    if num_of_qubits == 0:
        return matrix

    num_of_qubits -= 1
    dimension: int = 1 << num_of_qubits

    top_left: np.ndarray = matrix[0:dimension, 0:dimension]
    top_right: np.ndarray = matrix[0:dimension, dimension:]
    bottom_left: np.ndarray = matrix[dimension:, 0:dimension]
    bottom_right: np.ndarray = matrix[dimension:, dimension:]
    cmw_1: np.ndarray = top_left + bottom_right
    cmw_x: np.ndarray = bottom_left + top_right
    cmw_y: np.ndarray = -1.0j * (bottom_left - top_right)
    cmw_z: np.ndarray = top_left - bottom_right
    cmws: list[np.ndarray] = [cmw_1, cmw_x, cmw_y, cmw_z]

    for pauliInd,cmw in zip(PAULI_LIST,cmws):
        if cmw.any():
            PTM_rec_com: np.ndarray = recursive_commutator(cmw, num_of_qubits)
            PTM_single_anti: np.ndarray = SINGLE_ANTICOMMUTATOR[pauliInd]
            PTM += np.reshape(np.tensordot(PTM_single_anti, PTM_rec_com, axes=0), (PTM_dim, PTM_dim))

            if pauliInd != "I":
                PTM_rec_anti: np.ndarray = recursive_anticommutator(cmw, num_of_qubits)
                PTM_single_com: np.ndarray = SINGLE_COMMUTATOR[pauliInd]
                PTM += np.reshape(np.tensordot(PTM_single_com, PTM_rec_anti, axes=0), (PTM_dim, PTM_dim))
    return PTM


def anticommutator_to_ptm(matrix: np.ndarray):
    matrix_dim: int = matrix.shape[0]
    num_of_qubits: int = matrix_dim.bit_length() - 1
    PTM_dim: int = 1 << (2 * num_of_qubits)
    PTM: np.ndarray = recursive_anticommutator(matrix, num_of_qubits)
    return PTM

def recursive_anticommutator(matrix: np.ndarray, num_of_qubits: int):
    PTM_dim: int = 1 << (2*num_of_qubits)
    PTM: np.ndarray = np.zeros((PTM_dim, PTM_dim), dtype=np.complex64)

    if num_of_qubits == 0:
        return matrix

    num_of_qubits -= 1
    dimension: int = 1 << num_of_qubits

    top_left: np.ndarray = matrix[0:dimension, 0:dimension]
    top_right: np.ndarray = matrix[0:dimension, dimension:]
    bottom_left: np.ndarray = matrix[dimension:, 0:dimension]
    bottom_right: np.ndarray = matrix[dimension:, dimension:]
    cmw_1: np.ndarray = top_left + bottom_right
    cmw_x: np.ndarray = bottom_left + top_right
    cmw_y: np.ndarray = -1.0j * (bottom_left - top_right)
    cmw_z: np.ndarray = top_left - bottom_right
    cmws: list[np.ndarray] = [cmw_1, cmw_x, cmw_y, cmw_z]

    for pauliInd,cmw in zip(PAULI_LIST,cmws):
        if cmw.any():
            PTM_rec_anti: np.ndarray = recursive_anticommutator(cmw, num_of_qubits)
            PTM_single_anti: np.ndarray = SINGLE_ANTICOMMUTATOR[pauliInd]
            PTM += np.reshape(np.tensordot(PTM_single_anti, PTM_rec_anti, axes=0), (PTM_dim, PTM_dim))

            if pauliInd != "I":
                PTM_rec_com = recursive_commutator(cmw, num_of_qubits)
                PTM_single_com = SINGLE_COMMUTATOR[pauliInd]
                PTM += np.reshape(np.tensordot(PTM_single_com, PTM_rec_com, axes=0), (PTM_dim, PTM_dim))

    return PTM


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    num_of_qubits = 3
    data = test_matrices.rand_dense_mat(2 ** num_of_qubits)
    mat = commutator_to_ptm(data)
    print(mat)
