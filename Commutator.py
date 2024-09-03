import numpy as np

from qiskit.quantum_info import SparsePauliOp
from qiskit.quantum_info.operators import Operator, Pauli

import scipy.sparse as SS

import TestMatrices as TM

from PTM_utils import matrix_slice, matrix_embedding

import sys
np.set_printoptions(suppress=True,linewidth=sys.maxsize,threshold=sys.maxsize)

comX = np.array([[0,0,0,0],[0,0,0,0],[0,0,0,-1j],[0,0,1j,0]])
comY = np.array([[0,0,0,0],[0,0,0,1j],[0,0,0,0],[0,-1j,0,0]])
comZ = np.array([[0,0,0,0],[0,0,-1j,0],[0,1j,0,0],[0,0,0,0]])

single_commutator = {"X":comX,"Y":comY,"Z":comZ}

antiI = np.eye(4)
antiX = np.array([[0,1,0,0],[1,0,0,0],[0,0,0,0],[0,0,0,0]])
antiY = np.array([[0,0,1,0],[0,0,0,0],[1,0,0,0],[0,0,0,0]])
antiZ = np.array([[0,0,0,1],[0,0,0,0],[0,0,0,0],[1,0,0,0]])
single_anticommutator = {"I":antiI,"X":antiX,"Y":antiY,"Z":antiZ}

# Datastructure
pauliI = [1.,np.add,[0,3]]
pauliX = [1.,np.add,[1,2]]
pauliY = [1.j,np.subtract,[1,2]]
pauliZ = [1.,np.subtract,[0,3]]
paulis = {"I":pauliI,"X":pauliX,"Y":pauliY,"Z":pauliZ}
pauliList = ["I","X","Y","Z"]

def PTM_commutator(matrix, factor=1.):
	debug = False
	"""
	recursive commutators
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
			# Use the tensorproduct commutator rule.
			if rest.any():
				if pauliInd != "I":
					PTMra = PTM_anticommutator(rest,factor*2.j*pauli[0])
					PTMsc = single_commutator[pauliInd]
					PTM += np.reshape(np.tensordot(PTMsc,PTMra,axes=0),(PTMdim,PTMdim))
					del PTMra,PTMsc

				PTMrc = PTM_commutator(rest,factor*2*pauli[0])
				PTMsa = single_anticommutator[pauliInd]
				PTM += np.reshape(np.tensordot(PTMsa,PTMrc,axes=0),(PTMdim,PTMdim))
				del PTMrc, PTMsa
	return PTM

def PTM_anticommutator(matrix, factor=1.):
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
			# Use the tensorproduct anticommutator rule.
			if rest.any():
				PTMra = PTM_anticommutator(rest,factor*2.*pauli[0])
				PTMsa = single_anticommutator[pauliInd]
				PTM += np.reshape(np.tensordot(PTMsa,PTMra,axes=0),(PTMdim,PTMdim))
				del PTMra,PTMsa
				
				if pauliInd != "I":
					PTMrc = PTM_commutator(rest,factor*2.j*pauli[0])
					PTMsc = single_commutator[pauliInd]
					PTM += np.reshape(np.tensordot(PTMsc,PTMrc,axes=0),(PTMdim,PTMdim))
					del PTMrc, PTMsc
	return PTM
		

if __name__ == "__main__":
	qDim = 1
	data = TM.diagRandom(2**qDim)
	mat = PTM_commutator(data)
	print(mat)
