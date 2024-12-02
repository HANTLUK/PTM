import numpy as np
import numba as nb
import test_matrices
import sys

import logging
logger = logging.getLogger(__name__)

np.set_printoptions(suppress=True, linewidth=sys.maxsize, threshold=sys.maxsize)

PAULI_LIST: list = ["I","X","Y","Z"]

conjII = np.eye(4)
conjIX = np.array([[0,1,0,0],[1,0,0,0],[0,0,0,1j],[0,0,-1j,0]])
conjIY = np.array([[0,0,1,0],[0,0,0,-1j],[1,0,0,0],[0,1j,0,0]])
conjIZ = np.array([[0,0,0,1],[0,0,1j,0],[0,-1j,0,0],[1,0,0,0]])
conjXI = np.array([[0,1,0,0],[1,0,0,0],[0,0,0,-1j],[0,0,1j,0]])
conjXX = np.eye(4)
conjXY = np.array([[0,0,0,-1j],[0,1,0,0],[0,0,1,0],[1j,0,0,0]])
conjXZ = np.array([[0,0,1j,0],[0,0,0,1],[-1j,0,0,0],[0,1,0,0]])
conjYI = np.array([[0,0,1,0],[0,0,0,1j],[1,0,0,0],[0,-1j,0,0]])
conjYX = np.array([[0,0,0,1j],[0,1,0,0],[0,0,1,0],[-1j,0,0,0]])
conjYY = np.eye(4)
conjYZ = np.array([[0,-1j,0,0],[1j,0,0,0],[0,0,0,1],[0,0,1,0]])
conjZI = np.array([[0,0,0,1],[0,-1j,0,0],[0,0,1j,0],[1,0,0,0]])
conjZX = np.array([[0,0,-1j,0],[0,0,0,1],[1j,0,0,0],[0,1,0,0]])
conjZY = np.array([[0,1j,0,0],[-1j,0,0,0],[0,0,0,1],[0,0,1,0]])
conjZZ = np.eye(4)

SINGLE_CONJUGATION: dict = {"II":conjII,
                        "IX":conjIX,
                        "IY":conjIY,
                        "IZ":conjIZ,
                        "XI":conjXI,
                        "XX":conjXX,
                        "XY":conjXY,
                        "XZ":conjXZ,
                        "YI":conjYI,
                        "YX":conjYX,
                        "YY":conjYY,
                        "YZ":conjYZ,
                        "ZI":conjZI,
                        "ZX":conjZX,
                        "ZY":conjZY,
                        "ZZ":conjZZ}

@nb.jit(nopython=True)
def sc_any(array):
    for x in array.flat:
        if x:
            return True
    return False

def sandwich_to_ptm(A: np.ndarray, B: np.ndarray = None):
    if B is None:
        B = A
    matrix_dim: int = A.shape[0]
    num_of_qubits: int = matrix_dim.bit_length() - 1
    PTM_dim: int = 1 << (2*num_of_qubits)
    PTM: np.ndarray = recursive_sandwich(A, B, num_of_qubits)
    return PTM

def recursive_sandwich(A: np.ndarray, B: np.ndarray, num_of_qubits: int):
    PTM_dim: int = 1 << (2*num_of_qubits)
    PTM: np.ndarray = np.zeros((PTM_dim,PTM_dim), dtype=np.complex64)
    if num_of_qubits == 0:
        return A*np.conj(B)
    num_of_qubits -= 1
    dimension: int = 1 << num_of_qubits
    matrix: np.ndarray = A
    top_left: np.ndarray = matrix[0:dimension, 0:dimension]
    top_right: np.ndarray = matrix[0:dimension, dimension:]
    bottom_left: np.ndarray = matrix[dimension:, 0:dimension]
    bottom_right: np.ndarray = matrix[dimension:, dimension:]
    cmw_1: np.ndarray = top_left + bottom_right
    cmw_x: np.ndarray = bottom_left + top_right
    cmw_y: np.ndarray = -1.0j * (bottom_left - top_right)
    cmw_z: np.ndarray = top_left - bottom_right
    cmws_A: list[np.ndarray] = [cmw_1, cmw_x, cmw_y, cmw_z]

    matrix = B
    top_left: np.ndarray = matrix[0:dimension, 0:dimension]
    top_right: np.ndarray = matrix[0:dimension, dimension:]
    bottom_left: np.ndarray = matrix[dimension:, 0:dimension]
    bottom_right: np.ndarray = matrix[dimension:, dimension:]
    cmw_1: np.ndarray = top_left + bottom_right
    cmw_x: np.ndarray = bottom_left + top_right
    cmw_y: np.ndarray = -1.0j * (bottom_left - top_right)
    cmw_z: np.ndarray = top_left - bottom_right
    cmws_B: list[np.ndarray] = [cmw_1, cmw_x, cmw_y, cmw_z]

    for cmw_A, indexA in zip(cmws_A, PAULI_LIST):
        if sc_any(cmw_A):
            for cmw_B, indexB in zip(cmws_B, PAULI_LIST):
                if sc_any(cmw_B):
                    PTMsc: np.ndarray = SINGLE_CONJUGATION[indexA+indexB]
                    PTMrc: np.ndarray = recursive_sandwich(cmw_A,cmw_B, num_of_qubits)
                    PTM += np.real(np.kron(PTMsc,PTMrc))
    return PTM

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    num_of_qubits: int = 3
    data1: np.ndarray = test_matrices.rand_dense_mat(2 ** num_of_qubits)
    mat: np.ndarray = sandwich_to_ptm(data1)
    print(f"{mat} \n {mat.shape}")
