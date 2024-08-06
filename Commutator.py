import numpy as np

from qiskit.quantum_info import SparsePauliOp
from qiskit.quantum_info.operators import Operator, Pauli

import scipy.sparse as SS

import TestMatrices as TM

def single_commutator(pauli):
	if pauli == "I":
		return []
	elif pauli == "X":
		return [[2,3],[3,2]]
	elif pauli == "Y":
		return [[3,1],[1,3]]
	elif pauli == "Z":
		return [[1,2],[2,1]]

def single_anticommutator(pauli):
	if pauli == "I":
		return [[0,0],[1,1],[2,2],[3,3]]
	elif pauli == "X":
		return [[0,1],[1,0]]
	elif pauli == "Y":
		return [[0,2],[2,0]]
	elif pauli == "Z":
		return [[0,3],[3,0]]

# Datastructure
pauliI = ["I",1.,np.add,[[0,0],[1,1]]]
pauliX = ["X",1.,np.add,[[0,1],[1,0]]]
pauliY = ["Y",1.j,np.subtract,[[0,1],[1,0]]]
pauliZ = ["Z",1.,np.subtract,[[0,0],[1,1]]]
paulis = [pauliI,pauliX,pauliY,pauliZ]

def matrix_slice(matrix,coords):
	debug = False
	dim = matrix.shape[0]
	bDim = int(dim/2)
	if debug: print("Range 1", range(coords[0]*bDim,(coords[0]+1)*bDim))
	if debug: print("Range 2", range(coords[1]*bDim,(coords[1]+1)*bDim))
	return matrix[coords[0]*bDim:(coords[0]+1)*bDim,coords[1]*bDim:(coords[1]+1)*bDim]

def matrix_embedding(sparse,coords,eDim):
	debug = False
	row,col,data = sparse
	data = data*len(coords)
	row_ind = []
	col_ind = []
	for coord in coords:
		row_ind += [rowi + coord[0]*eDim for rowi in row]
		col_ind += [coli + coord[1]*eDim for coli in col]
	return [row_ind,col_ind,data]

def PTM_commutator(matrix, factor=1., output="mat"):
	debug = False
	"""
	recursive commutators
	"""
	dim = matrix.shape[0]
	qDim = dim.bit_length() - 1
	row_ind = []
	col_ind = []
	data = []
	if qDim == 0:
		return [[0],[0],[matrix[0][0]*factor]]
	else:
		for pauli in paulis:
			rest = pauli[2](matrix_slice(matrix,pauli[3][0]), matrix_slice(matrix,pauli[3][1]))
			# Use the tensorproduct commutator rule.
			if rest.any():
				if pauli[0] != "I":
					PTMra = PTM_anticommutator(rest,factor*2.j*pauli[1],None)
					PTMsc = single_commutator(pauli[0])
					nrow,ncol,ndata = matrix_embedding(PTMra,PTMsc,4**(qDim-1))
					row_ind += nrow
					col_ind += ncol
					data += ndata
					del PTMra,PTMsc

				PTMrc = PTM_commutator(rest,factor*2*pauli[1],None)
				PTMsa = single_anticommutator(pauli[0])
				nrow,ncol,ndata = matrix_embedding(PTMrc,PTMsa,4**(qDim-1))
				row_ind += nrow
				col_ind += ncol
				data += ndata
				del PTMrc, PTMsa
	if output == "mat":
		return SS.csc_array((data, (row_ind, col_ind)), shape=(4**qDim,4**qDim)).toarray()

	return [row_ind,col_ind,data]

def PTM_anticommutator(matrix, factor=1., output="mat"):
	"""
	recursive anticommutators
	"""
	dim = matrix.shape[0]
	qDim = dim.bit_length() - 1
	row_ind = []
	col_ind = []
	data = []
	if qDim == 0:
		return [[0],[0],[matrix[0][0]*factor]]
	else:
		for pauli in paulis:
			rest = pauli[2](matrix_slice(matrix,pauli[3][0]), matrix_slice(matrix,pauli[3][1]))
			# Use the tensorproduct anticommutator rule.
			if rest.any():
				PTMra = PTM_anticommutator(rest,factor*2.*pauli[1],None)
				PTMsa = single_anticommutator(pauli[0])
				nrow,ncol,ndata = matrix_embedding(PTMra,PTMsa,4**(qDim-1))
				row_ind += nrow
				col_ind += ncol
				data += ndata
				del PTMra,PTMsa

				PTMrc = PTM_commutator(rest,factor*2.j*pauli[1],None)
				PTMsc = single_commutator(pauli[0])
				nrow,ncol,ndata = matrix_embedding(PTMrc,PTMsc,4**(qDim-1))
				row_ind += nrow
				col_ind += ncol
				data += ndata
				del PTMrc, PTMsc

	if output == "mat":
		return SS.csc_array((data, (row_ind, col_ind)), shape=(4**qDim,4**qDim)).toarray()

	return [row_ind,col_ind,data]
		

if __name__ == "__main__":
	qDim = 1
	data = TM.diagRandom(2**qDim)
	mat = PTM_commutator(data)
	print(mat)
