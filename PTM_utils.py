import numpy as np

def matrix_slice(matrix):
	debug = False
	dim = matrix.shape[0]
	bDim = int(dim/2)
	coords = [[0,0],[0,1],[1,0],[1,1]]
	slices = [matrix[coord[0]*bDim:(coord[0]+1)*bDim,coord[1]*bDim:(coord[1]+1)*bDim] for coord in coords]
	return slices

def matrix_embedding(sparse,sparseEmb,eDim,factor=None):
	debug = False
	if debug: print(sparseEmb)
	comps,coords = sparseEmb
	row,col,dat = sparse
	if factor is not None:
		comps = [comp*factor for comp in comps]
	row_ind = [rowi + coord[0]*eDim for rowi in row for coord in coords]
	col_ind = [coli + coord[1]*eDim for coli in col for coord in coords]
	data = [dati*comp for dati in dat for comp in comps]
	return [row_ind,col_ind,data]
