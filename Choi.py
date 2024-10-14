import numpy as np
import math

from qiskit.quantum_info import SparsePauliOp
from qiskit.quantum_info.operators import Operator, Pauli

import scipy.sparse as SS

import TestMatrices as TM

from PTM_utils import matrix_slice, matrix_embedding

from qiskit.quantum_info.operators.channel import Choi, PTM

import sys
np.set_printoptions(suppress=True,linewidth=sys.maxsize,threshold=sys.maxsize)

import cProfile

def PTM_Choi(matrix: np.ndarray):
	matrix_dim: int = matrix.shape[0]
	num_of_qubits: int = math.ceil(np.log(matrix_dim)/np.log(2))
	half_qubits: int = int(num_of_qubits/2)
	PTM: np.ndarray = np.zeros((matrix_dim,matrix_dim),dtype=np.complex64)
	recursive_Choi(matrix, num_of_qubits, PTM, half_qubits)
	return PTM

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
		cmw_1: np.ndarray = 0.5 * (top_left + bottom_right)
		cmw_x: np.ndarray = 0.5 * (bottom_left + top_right)
		cmw_y: np.ndarray = -0.5j * (bottom_left - top_right)
		cmw_z: np.ndarray = 0.5 * (top_left - bottom_right)
		cmws: list[np.ndarray] = [cmw_1, cmw_x, cmw_y, cmw_z]
		

		index: int = 0
		for cmw in cmws:
			# if cmw.any():
			recursive_Choi(cmw, num_of_qubits,PTM,half_qubits,index_col,index_row+index)
			index += halved_dim**2
	
	elif num_of_qubits == 0:
		cmw_1: np.ndarray = 0.5 * (top_left + bottom_right)
		cmw_x: np.ndarray = 0.5 * (bottom_left + top_right)
		cmw_y: np.ndarray = -0.5j * (bottom_left - top_right)
		cmw_z: np.ndarray = 0.5 * (top_left - bottom_right)
		cmws: list[np.ndarray] = [cmw_1, cmw_x, cmw_y, cmw_z]

		index: int = 0
		for cmw in cmws:
			# if cmw.any():
			PTM[index_row+index][index_col] = cmw[0][0]
			index += 1
	return

if __name__ == "__main__":
	qDim = 6
	data = TM.denseRandom(4**qDim)
	cProfile.run("PTM_Choi(data)")

	ChoiMat = Choi(data)
	mat2 = PTM(ChoiMat).data
	# print("dif\n",np.array_str(mat2-mat, precision=1))

	
