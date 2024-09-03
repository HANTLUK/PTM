import numpy as np

from qiskit.quantum_info import SparsePauliOp
from qiskit.quantum_info.operators import Operator, Pauli

import scipy.sparse as SS

import TestMatrices as TM

from PTM_utils import matrix_slice, matrix_embedding

import sys
np.set_printoptions(suppress=True,linewidth=sys.maxsize,threshold=sys.maxsize)

# Datastructure
pauliI = [1.,np.add,[0,3]]
pauliX = [1.,np.add,[1,2]]
pauliY = [1.j,np.subtract,[1,2]]
pauliZ = [1.,np.subtract,[0,3]]
paulis = {"I":pauliI,"X":pauliX,"Y":pauliY,"Z":pauliZ}
pauliList = ["I","X","Y","Z"]

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


def PTM_conjugation(A,B=None):
	debug = False
	"""
	input: two matrices A,B output PTM(A . . B^dagger)
	"""
	dim = A.shape[0]
	qDim = dim.bit_length() - 1
	PTMdim = 1 << (2*qDim)
	PTM = np.zeros((PTMdim,PTMdim),dtype=np.complex64)
	if B is None:
		B = A
	if qDim == 0:
		return np.array([[A[0][0]*np.conj(B[0][0])]],dtype=np.complex64)
	else:
		SlicesA = matrix_slice(A)
		SlicesB = matrix_slice(B)
		for i in range(4):
			pauliIndI = pauliList[i]
			pauli = paulis[pauliIndI]
			restA = pauli[1](SlicesA[pauli[2][0]],SlicesA[pauli[2][1]])
			if restA.any():
				# Diagonal ones
				restB = pauli[1](SlicesB[pauli[2][0]],SlicesB[pauli[2][1]])
				if restB.any():
					PTM += np.reshape(np.tensordot(np.eye(4),PTM_conjugation(restA,restB),axes=0),(PTMdim,PTMdim))
				# Off-diagonal
				for j in range(i+1,4):
					pauliIndJ = pauliList[j]
					pauliJ = paulis[pauliIndJ]
					restB = pauliJ[1](SlicesB[pauliJ[2][0]],SlicesB[pauliJ[2][1]])
					if restB.any():
						PTMsc = single_conjugation[pauliIndI+pauliIndJ]
						PTMrc = PTM_conjugation(restA,restB)
						PTM += np.real(np.reshape(np.tensordot(PTMsc,PTMrc,axes=0),(PTMdim,PTMdim)))
						del PTMrc,PTMsc
	return PTM

if __name__ == "__main__":
	qDim = 1
	data1 = TM.denseRandom(2**qDim)
	mat = PTM_conjugation(data1)
	print(mat)
