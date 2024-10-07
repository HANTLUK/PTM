import numpy as np
import math
from concurrent.futures import ThreadPoolExecutor

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
	# matDim = matrix.shape[0]
	qBitDim = int(np.log(matrix.shape[0]) / np.log(4))
	W = np.apply_along_axis(TPD, 1, matrix, qBitDim)
	PTM = np.apply_along_axis(iTPD, 0, W, qBitDim)
	return PTM
        

if __name__ == "__main__":
	qDim = 2
	data = TM.denseRandom(4**qDim)
	print(np.array_str(data, precision=1))
	mat = PTM_SuperOp(data)

	SOp = SuperOp(data)
	mat2 = PTM(SOp).data
	dif = mat2-mat
	print("dif\n",np.array_str(dif, precision=1))
