import numpy as np
import math
import test_matrices
from qiskit.quantum_info.operators.channel import Choi, PTM
import sys

import logging
logger = logging.getLogger(__name__)

np.set_printoptions(suppress=True, linewidth=sys.maxsize, threshold=sys.maxsize)

def choi_to_ptm(matrix: np.ndarray):
	matrix_dim: int = matrix.shape[0]
	num_of_qubits: int = math.ceil(np.log(matrix_dim)/np.log(2))
	half_qubits: int = int(num_of_qubits/2)
	PTM: np.ndarray = np.zeros((matrix_dim,matrix_dim),dtype=np.complex64)
	recursive_Choi(matrix, num_of_qubits, PTM, half_qubits)
	return 2**(-half_qubits)*PTM

def recursive_Choi(matrix: np.ndarray, num_of_qubits: int, PTM: np.ndarray, half_qubits: int, index_col: int=0, index_row: int=0):
	num_of_qubits -= 1
	halved_dim: int = 1 << num_of_qubits
	top_left: np.ndarray = matrix[0:halved_dim, 0:halved_dim]
	top_right: np.ndarray = matrix[0:halved_dim, halved_dim:]
	bottom_left: np.ndarray = matrix[halved_dim:, 0:halved_dim]
	bottom_right: np.ndarray = matrix[halved_dim:, halved_dim:]

	if num_of_qubits >= half_qubits:
		cmw_1: np.ndarray = top_left + bottom_right
		cmw_2: np.ndarray = bottom_left + top_right
		cmw_3: np.ndarray = 1.0j * (bottom_left - top_right)
		cmw_4: np.ndarray = top_left - bottom_right
		cmws: list[np.ndarray] = [cmw_1, cmw_2, cmw_3, cmw_4]
		
		col_dim: int = 1 << (num_of_qubits - half_qubits)
		index: int = 0
		for cmw in cmws:
			# if cmw.any():
			recursive_Choi(cmw, num_of_qubits,PTM,half_qubits,index_col+index)
			index += col_dim**2

	elif num_of_qubits > 0:
		cmw_1: np.ndarray = (top_left + bottom_right)
		cmw_x: np.ndarray = (bottom_left + top_right)
		cmw_y: np.ndarray = -1.0j * (bottom_left - top_right)
		cmw_z: np.ndarray = (top_left - bottom_right)
		cmws: list[np.ndarray] = [cmw_1, cmw_x, cmw_y, cmw_z]
		

		index: int = 0
		for cmw in cmws:
			recursive_Choi(cmw, num_of_qubits,PTM,half_qubits,index_col,index_row+index)
			index += halved_dim**2
	
	elif num_of_qubits == 0:
		cmw_1: np.ndarray = (top_left + bottom_right)
		cmw_x: np.ndarray = (bottom_left + top_right)
		cmw_y: np.ndarray = -1.j * (bottom_left - top_right)
		cmw_z: np.ndarray = (top_left - bottom_right)
		cmws: list[np.ndarray] = [cmw_1, cmw_x, cmw_y, cmw_z]

		for i,cmw in enumerate(cmws):
			PTM[index_row+i][index_col] = cmw[0][0]
	return

if __name__ == "__main__":
    num_of_qubits = 2
    data = test_matrices.rand_dense_mat(4 ** num_of_qubits)
    mat = choi_to_ptm(data)
    choi_mat = Choi(data)
    mat2 = PTM(choi_mat).data
    print("mat\n", np.array_str(mat, precision=1))
    print("mat2\n", np.array_str(mat2, precision=1))
    print("dif\n", np.array_str(mat2 - mat, precision=1))
