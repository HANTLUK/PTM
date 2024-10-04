import numpy as np
import math

from qiskit.quantum_info import SparsePauliOp
from qiskit.quantum_info.operators import Operator, Pauli

import scipy.sparse as SS

import TestMatrices as TM

from PTM_utils import matrix_slice, matrix_embedding

from TPD import TPD,iTPD

from qiskit.quantum_info.operators.channel import SuperOp, PTM

import sys
np.set_printoptions(suppress=True,linewidth=sys.maxsize,threshold=sys.maxsize)

def PTM_SuperOp(matrix):
	debug = False
	matDim = matrix.shape[0]
	
	PTM = np.zeros((matDim,matDim),dtype=np.complex64)
	for i in range(matDim):
		if debug: print(f"\n Row {i}, TPD")	
		if debug: print(np.array_str(matrix[i,:], precision=1))
		PTM[i,:] = TPD(matrix[i,:])
		if debug: print(np.array_str(PTM[i,:], precision=1))
	W = np.zeros((matDim,matDim),dtype=np.complex64)
	for i in range(matDim):
		if debug: print(f"\n Col {i}, iTPD")
		if debug: print(np.array_str(PTM[:,i], precision=3))
		W[:,i] = iTPD(PTM[:,i])
		if debug: print(np.array_str(W[:,i], precision=3))
	return W
        

if __name__ == "__main__":
	qDim = 2
	data = TM.denseRandom(4**qDim)
	print(np.array_str(data, precision=1))
	mat = PTM_SuperOp(data)
	
	SOp = SuperOp(data)
	mat2 = PTM(SOp).data
	print("dif\n",np.array_str(mat2-mat, precision=1))
