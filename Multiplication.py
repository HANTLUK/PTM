import numpy as np

from qiskit.quantum_info import SparsePauliOp
from qiskit.quantum_info.operators import Operator, Pauli

import scipy.sparse as SS

import TestMatrices as TM

from PTM_utils import matrix_slice, matrix_embedding

mulI = np.eye(4)
mulX = np.array([[0,1,0,0],[1,0,0,0],[0,0,0,1j],[0,0,-1j,0]])
mulY = np.array([[0,0,1,0],[0,0,-1j,0],[1,0,0,0],[0,1j,0,0]])
mulZ = np.array([[0,0,0,1],[0,1j,0,0],[0,0,-1j,0],[1,0,0,0]])
single_lmultiplication = {"I":mulI,"X":mulX,"Y":mulY,"Z":mulZ}

# Datastructure
pauliI = [1.,np.add,[0,3]]
pauliX = [1.,np.add,[1,2]]
pauliY = [1.j,np.subtract,[1,2]]
pauliZ = [1.,np.subtract,[0,3]]
paulis = {"I":pauliI,"X":pauliX,"Y":pauliY,"Z":pauliZ}
pauliList = ["I","X","Y","Z"]


def PTM_multiplication(matrix, factor=1.):
	debug = False
	"""
	recursive anticommutators
	"""
	dim = matrix.shape[0]
	qDim = dim.bit_length() - 1
	PTMdim = 1 << (2*qDim)
	PTM = np.zeros((PTMdim,PTMdim),dtype=np.complex64)
	if qDim == 0:
		return np.array([matrix[0][0]*factor])
	else:
		Slices = matrix_slice(matrix)
		for pauliInd in pauliList:
			pauli = paulis[pauliInd]
			rest = pauli[1](Slices[pauli[2][0]], Slices[pauli[2][1]])
			# Use the tensorproduct multiplication rule.
			if rest.any():
				PTMrm = PTM_multiplication(rest,factor*pauli[0])
				if debug: print("PTMrm",PTMrm)
				PTMsm = single_lmultiplication[pauliInd]
				if debug: print("PTMsm",PTMsm)
				if debug: print("Kron",np.kron(PTMsm,PTMrm))
				PTM += np.reshape(np.tensordot(PTMsm,PTMrm,axes=0),(PTMdim,PTMdim))
				del PTMsm,PTMrm

	return PTM

if __name__ == "__main__":
	qDim = 7
	data = TM.denseRandom(2**qDim)
	PTM_multiplication(data)
