import numpy as np
import scipy.sparse as sp


def rand_diag_mat(dim):
    return np.add(np.zeros([dim, dim]), np.diag(np.add(np.random.rand(dim), 1.j * np.random.rand(dim))))


def rand_dense_mat(dim):
    return np.add(np.random.rand(dim, dim), 1.j * np.random.rand(dim, dim))


def identity_mat(dim):
    return np.identity(dim)


def rand_sparse_mat(dim):
    return sp.random(dim, dim, density=1. / dim, format="csr").todense()


def rand_symmetric_mat(dim):
    rand = np.random.rand(dim, dim)
    return np.add(rand, np.transpose(rand))


def rand_hermitian_mat(dim):
    rand = rand_dense_mat(dim)
    return np.add(rand, rand.conjugate().transpose())


def rand_one_factor_mat(dim):
    rand = rand_dense_mat(2)
    return np.kron(rand, np.eye(int(dim / 2)))
