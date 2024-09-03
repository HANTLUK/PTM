import numpy as np
import math

import TestMatrices as TM

import sys
np.set_printoptions(suppress=True,linewidth=sys.maxsize,threshold=sys.maxsize)

def TPD(matrix,flag=None):
	"""
		Computes the Pauli decomposition of a square matrix.

		Iteratively splits tensor factors off and decomposes those smaller
		matrices. This is done using submatrices of the original matrix.
		The Pauli strings are generated in each step.

		Args:
			matrix: Matrix to be decomposed
	"""

	matDim = matrix.shape[0]
	qBitDim = math.ceil(np.log(matDim)/np.log(2))
	if flag is None:
		flag = True

	decomposition = []	
	# Output for dimension 1

	if qBitDim == 0:
		decomposition = [matrix[0,0]]
		del matrix

	# Calculates the tensor product coefficients via the sliced submatrices.
	# If one of these components is zero that coefficient is ignored.

	if qBitDim > 0:
		halfDim = int(2**(qBitDim-1))

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
				subDec = TPD(mat,flag)
			else:
				subDec = np.zeros(halfDim)
			decomposition.extend(subDec)
	if flag:
		return np.array(decomposition)
	return decomposition

if __name__ == "__main__":
	qDim = 10
	mat = TM.diagRandom(2**qDim)
	print(TPD(mat))
