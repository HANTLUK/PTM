import numpy as np
import test_matrices

import logging
logger = logging.getLogger(__name__)

mulI: np.ndarray = np.eye(4)
mulX: np.ndarray = np.array([[0,1,0,0],[1,0,0,0],[0,0,0,1j],[0,0,-1j,0]])
mulY: np.ndarray = np.array([[0,0,1,0],[0,0,-1j,0],[1,0,0,0],[0,1j,0,0]])
mulZ: np.ndarray = np.array([[0,0,0,1],[0,1j,0,0],[0,0,-1j,0],[1,0,0,0]])

SINGLE_LMULTIPLICATION: dict = {"I":mulI,"X":mulX,"Y":mulY,"Z":mulZ}

pauliI: list = [1.,np.add,[0,3]]
pauliX: list = [1.,np.add,[1,2]]
pauliY: list = [1.j,np.subtract,[1,2]]
pauliZ: list = [1.,np.subtract,[0,3]]
PAULIS: dict = {"I":pauliI,"X":pauliX,"Y":pauliY,"Z":pauliZ}

PAULI_LIST: list = ["I","X","Y","Z"]


def multiplication_to_ptm(matrix: np.ndarray):
    matrix_dim: int = matrix.shape[0]
    num_of_qubits: int = matrix_dim.bit_length() - 1
    PTM_dim: int = 1 << (2*num_of_qubits)
    PTM: np.ndarray = recursive_multiplication(matrix, num_of_qubits)
    return PTM

def recursive_multiplication(matrix: np.ndarray, num_of_qubits: int, factor: float = 1.):
    num_of_qubits -= 1
    dimension: int = 1 << num_of_qubits
    if num_of_qubits == 0:
        return np.array([matrix[0][0]*factor])

    PTM_dim: int = 1 << (2*num_of_qubits)
    PTM: np.ndarray = np.zeros((PTM_dim,PTM_dim),dtype=np.complex64)
    top_left: np.ndarray = matrix[0:dimension, 0:dimension]
    top_right: np.ndarray = matrix[0:dimension, dimension:]
    bottom_left: np.ndarray = matrix[dimension:, 0:dimension]
    bottom_right: np.ndarray = matrix[dimension:, dimension:]
    cmw_1: np.ndarray = top_left + bottom_right
    cmw_x: np.ndarray = bottom_left + top_right
    cmw_y: np.ndarray = -1.0j * (bottom_left - top_right)
    cmw_z: np.ndarray = top_left - bottom_right
    cmws: list[np.ndarray] = [cmw_1, cmw_x, cmw_y, cmw_z]

    for pauliInd, cmw in zip(PAULI_LIST, cmws):
        pauli: list = PAULIS[pauliInd]
        if cmw.any():
            PTM_rec_mul: np.ndarray = recursive_multiplication(cmw, num_of_qubits, factor * pauli[0])
            PTM_single_mul: np.ndarray = SINGLE_LMULTIPLICATION[pauliInd]
            PTM += np.reshape(np.tensordot(PTM_single_mul,PTM_rec_mul,axes=0),(PTM_dim,PTM_dim))
    return PTM


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    num_of_qubits: int = 7
    data: np.ndarray = test_matrices.rand_dense_mat(2 ** num_of_qubits)
    mat: np.ndarray = multiplication_to_ptm(data)
    logger.info(f"{mat}")
