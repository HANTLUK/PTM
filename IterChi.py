import numpy as np

import logging
logger = logging.getLogger(__name__)

from qiskit.quantum_info import SparsePauliOp
from qiskit.quantum_info.operators import Operator, Pauli
from qiskit.quantum_info.operators.channel import Chi, PTM

import scipy.sparse as SS

import TestMatrices as TM

from PTM_utils import matrix_slice, matrix_embedding

import sys
np.set_printoptions(suppress=True,linewidth=sys.maxsize,threshold=sys.maxsize)

pauliList: list = ["I","X","Y","Z"]

conjII = [1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1]
conjIX = [0,1,0,0,1,0,0,0,0,0,0,1j,0,0,-1j,0]
conjIY = [0,0,1,0,0,0,0,-1j,1,0,0,0,0,1j,0,0]
conjIZ = [0,0,0,1,0,0,1j,0,0,-1j,0,0,1,0,0,0]

conjXI = [0,1,0,0,1,0,0,0,0,0,0,-1j,0,0,1j,0]
conjXX = [1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,-1]
conjXY = [0,0,0,-1j,0,0,1,0,0,1,0,0,1j,0,0,0]
conjXZ = [0,0,1j,0,0,0,0,1,-1j,0,0,0,0,1,0,0]

conjYI = [0,0,1,0,0,0,0,1j,1,0,0,0,0,-1j,0,0]
conjYX = [0,0,0,1j,0,0,1,0,0,1,0,0,-1j,0,0,0]
conjYY = [1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,-1]
conjYZ = [0,-1j,0,0,1j,0,0,0,0,0,0,1,0,0,1,0]

conjZI = [0,0,0,1,0,0,-1j,0,0,1j,0,0,1,0,0,0]
conjZX = [0,0,-1j,0,0,0,0,1,1j,0,0,0,0,1,0,0]
conjZY = [0,1j,0,0,-1j,0,0,0,0,0,0,1,0,0,1,0]
conjZZ = [1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1]

single_conjugation = [conjII,conjIX,conjIY,conjIZ,conjXI,conjXX,conjXY,conjXZ,conjYI,conjYX,conjYY,conjYZ,conjZI,conjZX,conjZY,conjZZ]

def PTM_Chi(matrix: np.ndarray):
	matrix_dim: int = matrix.shape[0]
	num_of_qubits: int = int((matrix_dim.bit_length() - 1)/2)
	output_dim: int = 1 << int(2*num_of_qubits)
	PTM: np.ndarray = np.zeros((output_dim,output_dim),dtype=np.complex64)
	recursive_Chi(matrix, num_of_qubits, PTM)
	return PTM
	
def recursive_Chi(matrix: np.ndarray, num_of_qubits: int, output: np.ndarray, row_index: int = 0, col_index: int = 0):
	num_of_qubits -= 1
	minor_dim: int = 1 << int(2*num_of_qubits)
	slices: list[np.ndarray] = [matrix[i*minor_dim:(i+1)*minor_dim,j*minor_dim:(j+1)*minor_dim] for j in range(4) for i in range(4)]
	cmwis: list[list[np.ndarray]] = [[sum([slices[l]*single_conjugation[4*i+j][l] for l in range(16)]) for j in range(4)] for i in range(4)]

	if num_of_qubits == 0:
		indexi: int = 0
		for i, cmwi in enumerate(cmwis):
			indexj: int = 0
			for j, cmwij in enumerate(cmwi):
				output[row_index+indexi][col_index+indexj] =  cmwij[0][0]
				indexj += 1
			indexi += 1

	else:
		indexi: int = 0
		for i, cmwi in enumerate(cmwis):
			indexj: int = 0
			for j, cmwij in enumerate(cmwi):
				recursive_Chi(cmwij, num_of_qubits, output, row_index+indexi, col_index+indexj)
				indexj += minor_dim
			indexi += minor_dim
	
	return

if __name__ == "__main__":
	logging.basicConfig(level=logging.DEBUG)
	qDim = 3
	data1 = TM.denseRandom(4**qDim)
	mat = PTM_Chi(data1)
	Chi = Chi(data1)
	mat2 = PTM(Chi).data
	logger.info(f"{np.array_str(4**(qDim/2)*mat2-mat, precision=1)}")
