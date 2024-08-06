import numpy as np

import timeit
import time

import matplotlib.pyplot as plt


from Commutator import PTM_commutator
import TestMatrices as TM

if __name__ == "__main__":
	maxQdim = 7
	times = []
	for qDim in range(1,maxQdim):
		data = TM.diagRandom(2**qDim)
		start_time = timeit.default_timer()
		PTM_commutator(data)
		times.append(timeit.default_timer() - start_time)
		print("Dim",qDim,timeit.default_timer() - start_time)
	
	fig, ax = plt.subplots()
	ax.scatter(range(1,maxQdim),times)
	ax.set_yscale('log')
	ax.set_title(f"PTM Commutator Diag")
	plt.savefig(f"PlotsPTM/PTMCommutatorDiag.png",dpi=150)
