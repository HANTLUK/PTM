import numpy as np
import scipy.sparse as sp


def diagRandom(dim):
	return np.add(np.zeros([dim,dim]),np.diag(np.add(np.random.rand(dim),1.j*np.random.rand(dim))))
def denseRandom(dim):
	return np.add(np.random.rand(dim,dim),1.j*np.random.rand(dim,dim))
def identityMatrix(dim):
	return np.identity(dim)
def sparseRandom(dim):
	return sp.random(dim,dim,density=1./dim,format="csr")
def symmRandom(dim):
	rand = randomMatrix(dim)
	return np.add(rand,np.transpose(rand))
def hermRandom(dim):
	rand = randomMatrix(dim)
	return np.add(rand,rand.conjugate().transpose())
def oneTerm(dim):
	rand = randomMatrix(2)
	return np.kron(rand,np.eye(int(dim/2)))
