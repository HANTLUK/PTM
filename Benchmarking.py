import numpy as np

import timeit
import time

import matplotlib.pyplot as plt


from Commutator import PTM_commutator,PTM_anticommutator
from Multiplication import PTM_multiplication
from Conjugation import PTM_conjugation
from PTM_utils import matrix_slice
import TestMatrices as TM

maps = {"Commutator" : PTM_commutator, "Anticommutator" : PTM_anticommutator, "Multiplication" : PTM_multiplication, "Conjugation" : PTM_conjugation,"Slicing":matrix_slice}
matrizes = {"Diag":TM.diagRandom,"Dense":TM.denseRandom,"Id":TM.identityMatrix}
maxQdims = {"Commutator" : 7, "Anticommutator" : 7, "Multiplication" : 7, "Conjugation" : 5, "Slicing": 14}

def benchmark(matrixType,mapType,maxQdim):
	output = True
	times = []
	for qDim in range(1,maxQdim):
		data = matrizes[matrixType](2**qDim)
		start_time = timeit.default_timer()
		maps[mapType](data)
		times.append(timeit.default_timer() - start_time)
		if output: print("Dim:",qDim,"\t\t Time:",timeit.default_timer() - start_time)
	return times

def benchmarkPlot(matrixType,mapType,maxQdim):
	times = benchmark(matrixType,mapType,maxQdim)
	fig, ax = plt.subplots()
	ax.plot(range(1,maxQdim),times)
	ax.set_yscale('log')
	ax.set_xlim([0,maxQdim])
	ax.set_ylim([10**(-6), 10**2])
	ax.set_title(f"PTM {mapType} {matrixType}")
	plt.savefig(f"PlotsPTM/PTM{mapType}{matrixType}.png",dpi=150)

if __name__ == "__main__":
	matrixTypes = ["Diag","Dense"]
	mapTypes = ["Commutator","Anticommutator","Multiplication","Conjugation","Slicing"] # 
	for matrixType in matrixTypes:
		print("\n"+matrixType)
		for mapType in mapTypes:
			maxQdim = maxQdims[mapType]
			print("\n"+mapType)
			benchmarkPlot(matrixType,mapType,maxQdim)
