import numpy as np
import math

import TestMatrices as TM

import sys

np.set_printoptions(suppress=True, linewidth=sys.maxsize, threshold=sys.maxsize)

def tpd(vector: np.ndarray, num_of_qubits: int):
    if num_of_qubits == 0:
        return np.array(vector[0])
    else:
        matrix_dim: int = 1 << num_of_qubits
        vector_dim: int = 1 << (2*num_of_qubits)
        matrix: np.ndarray = vector.reshape(matrix_dim, matrix_dim)
        pauli_weights: np.ndarray = np.zeros(vector_dim, dtype=np.complex64)
        recursive_tpd(matrix, num_of_qubits, pauli_weights)
        return pauli_weights

def recursive_tpd(matrix: np.ndarray, num_of_qubits: int, pauli_weights: np.ndarray, index_outer: int = 0):
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
    
    if num_of_qubits == 0:
        index: int = 0
        for cmw in cmws:
            if cmw.any():
                pauli_weights[index_outer+index] = cmw[0][0]
            index += 1
    else:
        index: int = 0
        for cmw in cmws:
            if cmw.any():
                recursive_tpd(cmw, num_of_qubits,pauli_weights,index_outer+index)
            index += halved_dim ** 2

    return

def itpd(vector: np.ndarray, num_of_qubits: int):
    if num_of_qubits == 0:
        return np.array(vector[0])
    else:
        matrix_dim: int = 1 << num_of_qubits
        vector_dim: int = 1 << (2*num_of_qubits)
        matrix: np.ndarray = vector.reshape(matrix_dim, matrix_dim)
        can_weights: np.ndarray = np.zeros(vector_dim, dtype=np.complex64)
        recursive_itpd(matrix, num_of_qubits, can_weights)
        return can_weights


def recursive_itpd(matrix: np.ndarray, num_of_qubits: int, can_weights: np.ndarray, index_outer: int = 0):
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
    
    if num_of_qubits == 0:
        index: int = 0
        for cmw in cmws:
            if cmw.any():
                can_weights[index_outer+index] = cmw[0][0]
            index += 1
    else:
        index: int = 0
        for cmw in cmws:
            if cmw.any():
                recursive_itpd(cmw, num_of_qubits,can_weights,index_outer+index)
            index += halved_dim ** 2

    return

if __name__ == "__main__":
    qDim = 1
    mat = TM.diagRandom(2 ** qDim)
    print(mat, mat.shape)
    res = itpd(mat, qDim)
    print(res)
