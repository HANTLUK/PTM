import numpy as np

import timeit
import time

import matplotlib.pyplot as plt

import TBD as TBD
import ExampleMatrizes as EM
import transformations as trafo
	

Reps = ["Choi","SuperOp","Operator","Kraus","Chi","PTM","Stinespring"]
to = {"Choi" : trafo._to_choi, "SuperOp" : trafo._to_superop, "Kraus" : trafo._to_kraus, "PTM" : trafo._to_ptm, "Stinespring" : trafo._to_stinespring, "Chi" : trafo._to_chi}

Test = ["PTM","Choi"]

def numKrauss(qDim):
	# Smaller than qDim**2-1
	return qDim

if __name__ == "__main__":
	maxQdim = 7
	times = []
	From = Test[1]
	To = Test[0]
	for qDim in range(1,maxQdim):
		if From == "SuperOp" or From == "Choi":
			data = EM.denseRandom(4**qDim)
		if From == "Kraus":
			k = numKrauss(qDim)
			data = [[EM.denseRandom(2**qDim) for i in range(k)],None]
		start_time = timeit.default_timer()
		to[To](From,data,2**qDim,2**qDim)
		times.append(timeit.default_timer() - start_time)
	
	fig, ax = plt.subplots()
	ax.scatter(range(1,maxQdim),times)
	ax.set_yscale('log')
	ax.set_title(f"{From} to {To}")
	plt.savefig(f"PlotsCHReps/{From}to{To}.png",dpi=150)
	
	

