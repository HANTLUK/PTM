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

def PTM_Choi(matrix, midDim=None):
	debug = False
	"""
		Calculates CMWs for n qubits, then TPD for the rows.
	"""

	matDim = matrix.shape[0]
	qBitDim = math.ceil(np.log(matDim)/np.log(2))
	halfDim = int(2**(qBitDim-1))

	# Flag for recursion start
	flag1 = False
	if midDim is None:
		flag1 = True
		midDim = int(2**(qBitDim/2))

	decomposition = []

	# Output for dimension 1

	if qBitDim == 0:
		decomposition = [matrix[0,0]]
		del matrix

	# Calculates the tensor product coefficients via the sliced submatrices.
	# If one of these components is zero that coefficient is ignored.

	if halfDim >= midDim:

		coeff1 = (matrix[0:halfDim, 0:halfDim]
						+ matrix[halfDim:, halfDim:])
		coeffX = (matrix[halfDim:, 0:halfDim]
						+matrix[0:halfDim, halfDim:])
		coeffY = (1.j*matrix[halfDim:, 0:halfDim]
						-1.j*matrix[0:halfDim, halfDim:])
		coeffZ = (matrix[0:halfDim, 0:halfDim]
						- matrix[halfDim:, halfDim:])

		coefficients = {"I": coeff1, "X": coeffX, "Y": coeffY, "Z": coeffZ}
		del matrix

		# Recursion for the Submatrices

		for i,c in enumerate(coefficients):
			mat = coefficients[c]
			if mat.any() != 0:
				subDec = PTM_Choi(mat,midDim)
			else:
				subDec = [[0. for _ in range(int(midDim**2))] for i in range(int(matDim**2/midDim**2))]
			if halfDim == midDim:
				decomposition.append(subDec)
			else:
				decomposition.extend(subDec)

	elif qBitDim > 0:

		coeff1 = 0.5*(matrix[0:halfDim, 0:halfDim]
						+ matrix[halfDim:, halfDim:])
		coeffX = 0.5*(matrix[halfDim:, 0:halfDim]
						+ matrix[0:halfDim, halfDim:])
		coeffY = -1.j*0.5*(matrix[halfDim:, 0:halfDim]
						- matrix[0:halfDim, halfDim:])
		coeffZ = 0.5*(matrix[0:halfDim, 0:halfDim]
						- matrix[halfDim:, halfDim:])

		coefficients = {"I": coeff1, "X": coeffX, "Y": coeffY, "Z": coeffZ}
		del matrix

		# Recursion for the Submatrices

		for i,c in enumerate(coefficients):
			mat = coefficients[c]
			if mat.any() != 0:
				subDec = PTM_Choi(mat,midDim)
			else:
				subDec = [0. for _ in range(int(halfDim**2))]
			if halfDim == midDim:
				decomposition.append(subDec)
			else:
				decomposition.extend(subDec)	

	if flag1:
		return np.array(decomposition).transpose()

	return decomposition

if __name__ == "__main__":
	qDim = 1
	data = TM.denseRandom(4**qDim)
	mat = PTM_Choi(data)

	ChoiMat = Choi(data)
	mat2 = PTM(ChoiMat).data
	print("mat\n",np.array_str(mat, precision=1))
	print("mat2\n",np.array_str(mat2, precision=1))
	print("dif\n",np.array_str(mat2-mat, precision=1))
