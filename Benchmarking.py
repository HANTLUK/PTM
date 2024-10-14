import numpy as np

import timeit
import time

import logging
logger = logging.getLogger(__name__)

import matplotlib.pyplot as plt

from Commutator import PTM_commutator, PTM_anticommutator
from Multiplication import PTM_multiplication
from Conjugation import PTM_conjugation
from Chi import PTM_Chi
from SuperOp import PTM_SuperOp, PTM_SuperOp_old
from Choi import PTM_Choi
from PTM_utils import matrix_slice
import TestMatrices as TM

import qiskit.quantum_info.operators.channel as CH # Choi, Chi, PTM, SuperOp

maps = {"Commutator": PTM_commutator,
		"Anticommutator": PTM_anticommutator,
		"Multiplication": PTM_multiplication,
		"Conjugation": PTM_conjugation,
		"Slicing": matrix_slice,
		"Chi": PTM_Chi,
		"SuperOp_old": PTM_SuperOp_old,
		"SuperOp": PTM_SuperOp,
		"Choi": PTM_Choi}

mapsQK = {"Chi": CH.Chi,
		"Choi": CH.Choi,
		"SuperOp": CH.SuperOp}

matrizes = {"Diag": TM.diagRandom,
			"Dense": TM.denseRandom,
			"Id": TM.identityMatrix}

maxQdims = {"Commutator": 7,
			"Anticommutator": 7,
			"Multiplication": 7,
			"Conjugation": 5,
			"Slicing": 14,
			"Chi": 6,
			"SuperOp_old": 7,
            "SuperOp": 8,
			"Choi": 7}


def benchmark(matrixType, mapType, maxQdim, SOp=False, QK=False):
	output = True
	times = []
	base = 2
	if SOp:
		base = 4
	if QK:
		logger.info(f"QK")
		logger.info(f"qubits \t time \t growth")
		for qDim in range(1, maxQdim):
			data = matrizes[matrixType](base ** qDim)
			data = mapsQK[mapType](data)
			start_time = timeit.default_timer()
			CH.PTM(data)
			times.append(timeit.default_timer() - start_time)
			if qDim > 1:
				logger.info(f"{qDim} \t {times[-1]:.2f} \t {(times[-1] / times[-2]):.2f}")
	else:
		logger.info(f"TPD")
		logger.info(f"qubits \t time \t growth")
		for qDim in range(1, maxQdim):
			data = matrizes[matrixType](base ** qDim)
			start_time = timeit.default_timer()
			maps[mapType](data)
			times.append(timeit.default_timer() - start_time)
			if qDim > 1:
				logger.info(f"{qDim} \t {times[-1]:.2f} \t {(times[-1] / times[-2]):.2f}")
	return times


def benchmarkPlot(matrixType, mapType, maxQdim, SOp=False):
	fig, ax = plt.subplots()
	ax.set_ylim([10 ** (-6), 10 ** 2])
	if SOp and mapType in mapsQK.keys():
		timesQK = benchmark(matrixType, mapType, maxQdim, SOp, True)
		ax.scatter(range(1, maxQdim), timesQK)
		ax.set_ylim([10 ** (-6), 10 ** 3])
	times = benchmark(matrixType, mapType, maxQdim, SOp)
	ax.scatter(range(1, maxQdim), times)
	ax.set_yscale('log')
	ax.set_title(f"PTM {mapType} {matrixType}")
	plt.savefig(f"PlotsPTM/PTM{mapType}{matrixType}.png", dpi=150)


if __name__ == "__main__":
	logging.basicConfig(level=logging.INFO)
	matrixTypes = ["Diag", "Dense"]
	mapTypes = []
	#	 ["Commutator","Anticommutator","Multiplication","Conjugation","Slicing"]
	SOpmapTypes = ["Chi"]  # SuperOp_old, "Chi"
	for matrixType in matrixTypes:
		logger.info(f"{matrixType}")
		for mapType in mapTypes:
			maxQdim = maxQdims[mapType]
			logger.info(f"{mapType}")
			benchmarkPlot(matrixType, mapType, maxQdim)
		for mapType in SOpmapTypes:
			maxQdim = maxQdims[mapType]
			logger.info(f"{mapType}")
			benchmarkPlot(matrixType, mapType, maxQdim, True)
