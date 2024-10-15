import numpy as np

import logging
logger = logging.getLogger(__name__)

from qiskit.quantum_info import SparsePauliOp
from qiskit.quantum_info.operators import Operator, Pauli
from qiskit.quantum_info.operators.channel import Chi, PTM

import scipy.sparse as SS

import test_matrices as TM

from PTM_utils import matrix_slice, matrix_embedding

import sys
np.set_printoptions(suppress=True,linewidth=sys.maxsize,threshold=sys.maxsize)

PAULI_LIST: list = ["I","X","Y","Z"]

conjII: list = np.array([1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1])
conjIX: list = np.array([0,1,0,0,1,0,0,0,0,0,0,1j,0,0,-1j,0])
conjIY: list = np.array([0,0,1,0,0,0,0,-1j,1,0,0,0,0,1j,0,0])
conjIZ: list = np.array([0,0,0,1,0,0,1j,0,0,-1j,0,0,1,0,0,0])

conjXI: list = np.array([0,1,0,0,1,0,0,0,0,0,0,-1j,0,0,1j,0])
conjXX: list = np.array([1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,-1])
conjXY: list = np.array([0,0,0,-1j,0,0,1,0,0,1,0,0,1j,0,0,0])
conjXZ: list = np.array([0,0,1j,0,0,0,0,1,-1j,0,0,0,0,1,0,0])

conjYI: list = np.array([0,0,1,0,0,0,0,1j,1,0,0,0,0,-1j,0,0])
conjYX: list = np.array([0,0,0,1j,0,0,1,0,0,1,0,0,-1j,0,0,0])
conjYY: list = np.array([1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,-1])
conjYZ: list = np.array([0,-1j,0,0,1j,0,0,0,0,0,0,1,0,0,1,0])

conjZI: list = np.array([0,0,0,1,0,0,-1j,0,0,1j,0,0,1,0,0,0])
conjZX: list = np.array([0,0,-1j,0,0,0,0,1,1j,0,0,0,0,1,0,0])
conjZY: list = np.array([0,1j,0,0,-1j,0,0,0,0,0,0,1,0,0,1,0])
conjZZ: list = np.array([1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1])

SINGLE_CONJUGATION: list[np.ndarray] = [conjII,conjIX,conjIY,conjIZ,conjXI,conjXX,conjXY,conjXZ,conjYI,conjYX,conjYY,conjYZ,conjZI,conjZX,conjZY,conjZZ]

def chi_to_ptm(matrix: np.ndarray):
	matrix_dim: int = matrix.shape[0]
	num_of_qubits: int = int((matrix_dim.bit_length() - 1)/2)
	output_dim: int = 1 << int(2*num_of_qubits)
	PTM: np.ndarray = np.zeros((output_dim,output_dim),dtype=np.complex64)
	recursive_Chi(matrix, num_of_qubits, PTM)
	return PTM
	
def recursive_Chi(matrix: np.ndarray, num_of_qubits: int, output: np.ndarray, row_index: int = 0, col_index: int = 0):
	num_of_qubits -= 1
	minor_dim: int = 1 << int(2*num_of_qubits)
	slices: np.ndarray = np.array([matrix[i*minor_dim:(i+1)*minor_dim,j*minor_dim:(j+1)*minor_dim] for j in range(4) for i in range(4)])
	cmwis: list[list[np.ndarray]] = [[np.einsum("i,ijk->jk",SINGLE_CONJUGATION[4*i+j],slices) for j in range(4)] for i in range(4)]

	if num_of_qubits == 0:
		output[row_index:row_index+4][col_index:col_index+4] = cmwis[:][:][0][0]

	else:
		indexi: int = 0
		for cmwi in cmwis:
			indexj: int = 0
			for cmwij in cmwi:
				recursive_Chi(cmwij, num_of_qubits, output, row_index+indexi, col_index+indexj)
				indexj += minor_dim
			indexi += minor_dim
	
	return

if __name__ == "__main__":
	logging.basicConfig(level=logging.DEBUG)
	qDim = 3
	data1 = TM.denseRandom(4**qDim)
	mat = chi_to_ptm(data1)
	Chi = Chi(data1)
	mat2 = PTM(Chi).data
	logger.info(f"{np.array_str(4**(qDim/2)*mat2-mat, precision=1)}")
