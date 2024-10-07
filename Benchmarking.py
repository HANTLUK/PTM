import numpy as np

import timeit
import time

import matplotlib.pyplot as plt


from Commutator import PTM_commutator,PTM_anticommutator
from Multiplication import PTM_multiplication
from Conjugation import PTM_conjugation
from Chi import PTM_Chi
from SuperOp import PTM_SuperOp
# from Choi import PTM_Choi
from PTM_utils import matrix_slice
import TestMatrices as TM

maps = {"Commutator" : PTM_commutator, "Anticommutator" : PTM_anticommutator, "Multiplication" : PTM_multiplication, "Conjugation" : PTM_conjugation,"Slicing":matrix_slice,"Chi":PTM_Chi, "SuperOp": PTM_SuperOp}
matrizes = {"Diag":TM.diagRandom,"Dense":TM.denseRandom,"Id":TM.identityMatrix}
maxQdims = {"Commutator" : 7, "Anticommutator" : 7, "Multiplication" : 7, "Conjugation" : 5, "Slicing": 14, "Chi" : 7, "SuperOp": 7}

def benchmark(matrixType,mapType,maxQdim,SOp=False):
	output = True
	times = []
	base = 2
	if SOp: base = 4
	for qDim in range(0,maxQdim):
		data = matrizes[matrixType](base**qDim)
		start_time = timeit.default_timer()
		maps[mapType](data)
		times.append(timeit.default_timer() - start_time)
		if output and qDim > 0: print(f"{qDim} & {times[-1]} & {times[-1]/times[-2]}")
	return times

def benchmarkPlot(matrixType,mapType,maxQdim,SOp=False):
	times = benchmark(matrixType,mapType,maxQdim,SOp)
	fig, ax = plt.subplots()
	ax.scatter(range(1,maxQdim),times[1:])
	ax.set_yscale('log')
	ax.set_ylim([10**(-6), 10**2])
	ax.set_title(f"PTM {mapType} {matrixType}")
	plt.savefig(f"PlotsPTM/PTM{mapType}{matrixType}.png",dpi=150)

if __name__ == "__main__":
	matrixTypes = ["Diag","Dense"]
	mapTypes = []
#	 ["Commutator","Anticommutator","Multiplication","Conjugation","Slicing"]
	SOpmapTypes = ["SuperOp"] #
	for matrixType in matrixTypes:
		print("\n"+matrixType)
		for mapType in mapTypes:
			maxQdim = maxQdims[mapType]
			print("\n"+mapType)
			benchmarkPlot(matrixType,mapType,maxQdim)
		for mapType in SOpmapTypes:
			maxQdim = maxQdims[mapType]
			print("\n"+mapType)
			benchmarkPlot(matrixType,mapType,maxQdim,True)
