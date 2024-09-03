import numpy as np
import math

from qiskit.quantum_info import SparsePauliOp
from qiskit.quantum_info.operators import Operator, Pauli

import scipy.sparse as SS

import TestMatrices as TM

from PTM_utils import matrix_slice, matrix_embedding

import sys
np.set_printoptions(suppress=True,linewidth=sys.maxsize,threshold=sys.maxsize)

def PTM_SuperOp(matrix):
	debug = False
	matDim = matrix.shape[0]
	qBitDim = int(np.log(matDim)/np.log(4))
	mat = matrix.reshape(2**qBitDim,2**qBitDim,4**qBitDim)
	RowTPD = TPD_Rows(mat)
	if debug: print("Row TPD",RowTPD)
	mat = RowTPD.reshape(4**qBitDim,2**qBitDim,2**qBitDim)
	ColTPD = TPD_Cols(mat)
	return ColTPD

def TPD_Rows(matrix):
	debug = False
	matDim = matrix.shape[0]
	colDim = matrix.shape[2]
	qBitDim = math.ceil(np.log(matDim)/np.log(2))
	decomposition = []	
	
	# Output for dimension 1

	if qBitDim == 0:
		decomposition = np.array(matrix[0,0,:])
		del matrix

	# Calculates the tensor product coefficients via the sliced submatrices.
	# If one of these components is zero that coefficient is ignored.

	if qBitDim > 0:
		halfDim = int(2**(qBitDim-1))

		coeff1 = 0.5*(matrix[0:halfDim, 0:halfDim,:]
						+ matrix[halfDim:, halfDim:,:])
		coeffX = 0.5*(matrix[halfDim:, 0:halfDim,:]
						+ matrix[0:halfDim, halfDim:,:])
		coeffY = -1.j*0.5*(matrix[halfDim:, 0:halfDim,:]
						- matrix[0:halfDim, halfDim:,:])
		coeffZ = 0.5*(matrix[0:halfDim, 0:halfDim,:]
						- matrix[halfDim:, halfDim:,:])

		coefficients = {"I": coeff1, "X": coeffX, "Y": coeffY, "Z": coeffZ}
		del matrix

		# Recursion for the Submatrices

		for i,c in enumerate(coefficients):
			mat = coefficients[c]
			if mat.any() != 0:
				subDec = TPD_Rows(mat)
				if debug: print("SubDec",subDec)
			else:
				subDec = np.zeros(halfDim,colDim)
				if debug: print("SubDec",subDec)
			decomposition.append(subDec)

	return np.array(decomposition)

def TPD_Cols(matrix):
	matDim = matrix.shape[1]
	rowDim = matrix.shape[0]
	qBitDim = math.ceil(np.log(matDim)/np.log(2))
	decomposition = []
	# Output for dimension 1

	if qBitDim == 0:
		decomposition = np.array(matrix[:,0,0])
		del matrix

	# Calculates the tensor product coefficients via the sliced submatrices.
	# If one of these components is zero that coefficient is ignored.

	if qBitDim > 0:
		halfDim = int(2**(qBitDim-1))

		coeff1 = 0.5*(matrix[:,0:halfDim, 0:halfDim]
						+ matrix[:,halfDim:, halfDim:])
		coeffX = 0.5*(matrix[:,halfDim:, 0:halfDim]
						+ matrix[:,0:halfDim, halfDim:])
		coeffY = -1.j*0.5*(matrix[:,halfDim:, 0:halfDim]
						- matrix[:,0:halfDim, halfDim:])
		coeffZ = 0.5*(matrix[:,0:halfDim, 0:halfDim]
						- matrix[:,halfDim:, halfDim:])

		coefficients = {"I": coeff1, "X": coeffX, "Y": coeffY, "Z": coeffZ}
		del matrix

		# Recursion for the Submatrices

		for i,c in enumerate(coefficients):
			mat = coefficients[c]
			if mat.any() != 0:
				subDec = TPD_Cols(mat)
			else:
				subDec = np.zeros(rowDim,halfDim)
			decomposition.append(subDec)

	return np.array(decomposition)

if __name__ == "__main__":
	qDim = 6
	mat = TM.diagRandom(4**qDim)
	PTM_SuperOp(mat)
