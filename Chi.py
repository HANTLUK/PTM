import numpy as np

from qiskit.quantum_info import SparsePauliOp
from qiskit.quantum_info.operators import Operator, Pauli
from qiskit.quantum_info.operators.channel import Chi, PTM

import scipy.sparse as SS

import TestMatrices as TM

from PTM_utils import matrix_slice, matrix_embedding

import sys
np.set_printoptions(suppress=True,linewidth=sys.maxsize,threshold=sys.maxsize)

pauliList = ["I","X","Y","Z"]

conjII = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
conjIX = np.array([[0,1,0,0],[1,0,0,0],[0,0,0,1j],[0,0,-1j,0]])
conjIY = np.array([[0,0,1,0],[0,0,0,-1j],[1,0,0,0],[0,1j,0,0]])
conjIZ = np.array([[0,0,0,1],[0,0,1j,0],[0,-1j,0,0],[1,0,0,0]])

conjXI = np.array([[0,1,0,0],[1,0,0,0],[0,0,0,-1j],[0,0,1j,0]])
conjXX = np.array([[1,0,0,0],[0,1,0,0],[0,0,-1,0],[0,0,0,-1]])
conjXY = np.array([[0,0,0,-1j],[0,0,1,0],[0,1,0,0],[1j,0,0,0]])
conjXZ = np.array([[0,0,1j,0],[0,0,0,1],[-1j,0,0,0],[0,1,0,0]])

conjYI = np.array([[0,0,1,0],[0,0,0,1j],[1,0,0,0],[0,-1j,0,0]])
conjYX = np.array([[0,0,0,1j],[0,0,1,0],[0,1,0,0],[-1j,0,0,0]])
conjYY = np.array([[1,0,0,0],[0,-1,0,0],[0,0,1,0],[0,0,0,-1]])
conjYZ = np.array([[0,-1j,0,0],[1j,0,0,0],[0,0,0,1],[0,0,1,0]])

conjZI = np.array([[0,0,0,1],[0,0,-1j,0],[0,1j,0,0],[1,0,0,0]])
conjZX = np.array([[0,0,-1j,0],[0,0,0,1],[1j,0,0,0],[0,1,0,0]])
conjZY = np.array([[0,1j,0,0],[-1j,0,0,0],[0,0,0,1],[0,0,1,0]])
conjZZ = np.array([[1,0,0,0],[0,-1,0,0],[0,0,-1,0],[0,0,0,1]])

single_conjugation = {"II":conjII,
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


def PTM_Chi(matrix):
	debug = False
	"""
	input: Chi matrix
	"""
	dim = matrix.shape[0]
	qDim = int((dim.bit_length() - 1)/2)
	if debug: print("Dim",qDim)
	sliceDim = 4**(qDim - 1)
	PTMdim = 1 << (2*qDim)
	PTM = np.zeros((PTMdim,PTMdim),dtype=np.complex64)
	if qDim == 0:
		return matrix[0][0]
	else:
		for i,pauliIndI in enumerate(pauliList):
			for j,pauliIndJ in enumerate(pauliList):
				if debug: print("i,j",i,j)
				rest = matrix[i*sliceDim:(i+1)*sliceDim,j*sliceDim:(j+1)*sliceDim]
				PTMsc = single_conjugation[pauliIndI+pauliIndJ]
				if debug: print(PTMsc)
				PTMrc = PTM_Chi(rest)
				PTM += np.kron(PTMsc,PTMrc)
				if debug: print(np.array_str(PTM, precision=1))
				del PTMrc,PTMsc
	return PTM

if __name__ == "__main__":
	qDim = 5
	data1 = TM.denseRandom(4**qDim)
	mat = PTM_Chi(data1)
	Chi = Chi(data1)
	mat2 = PTM(Chi).data
	print(np.array_str(4**(qDim/2)*mat2-mat, precision=1))