import timeit
import matplotlib.pyplot as plt
from Commutator import commutator_to_ptm, anticommutator_to_ptm
from Multiplication import multiplication_to_ptm
from Conjugation import sandwich_to_ptm
from Chi import chi_to_ptm
from SuperOp import can_to_ptm, can_to_ptm_old, iter_can_to_ptm
from Choi import choi_to_ptm
from PTM_utils import matrix_slice
import test_matrices
import qiskit.quantum_info.operators.channel as ch  # Choi, Chi, PTM, SuperOp

MAPS = {"Commutator": commutator_to_ptm,
        "Anticommutator": anticommutator_to_ptm,
        "Multiplication": multiplication_to_ptm,
        "Conjugation": sandwich_to_ptm,
        "Slicing": matrix_slice,
        "Chi to PTM": chi_to_ptm,
        "canonical to PTM": can_to_ptm,
        "canonical to PTM (old)": can_to_ptm_old,
        "canonical to PTM (iterative)": iter_can_to_ptm,
        "Choi to PTM": choi_to_ptm}

MAPSQK = {"Chi to PTM": ch.Chi,
          "Choi to PTM": ch.Choi,
          "canonical to PTM": ch.SuperOp,
          "canonical to PTM (parallel)": ch.SuperOp}

MATRICES = {"random diagonal": test_matrices.rand_diag_mat,
            "random dense": test_matrices.rand_dense_mat,
            "identity": test_matrices.identity_mat}

MAX_NUM_OF_QUBITS = {"Commutator": 7,
                     "Anticommutator": 7,
                     "Multiplication": 7,
                     "Conjugation": 5,
                     "Slicing": 14,
                     "Chi to PTM": 7,
                     "canonical to PTM": 7,
                     "canonical to PTM (old)": 7,
                     "canonical to PTM (iterative)": 7,
                     "Choi to PTM": 6}


def benchmark(matrix_type, map_type, max_num_of_qubits, SOp=False, QK=False, verbose=False):
    times = []
    base = 4 if SOp else 2
    if QK:
        print("|   |-> \033[4mQiskit\033[0m")
    else:
        print("|   |-> \033[4mTPD-based\033[0m")

    for num_of_qubits in range(max_num_of_qubits + 1):
        data = MATRICES[matrix_type](base) if num_of_qubits == 0 else MATRICES[matrix_type](base ** num_of_qubits)
        if QK:
            start_time = timeit.default_timer()
            data = MAPSQK[map_type](data)
            ch.PTM(data)
            times.append(timeit.default_timer() - start_time)
        else:
            start_time = timeit.default_timer()
            MAPS[map_type](data)
            times.append(timeit.default_timer() - start_time)
        if verbose:
            if num_of_qubits == 1:
                print(f"|   |   {num_of_qubits}: {times[-1]}")
            elif num_of_qubits > 1:
                print(f"|   |   {num_of_qubits}: {times[-1]} (x{times[-1] / times[-2]})")
    return times


def benchmark_can_to_ptm(matrix_type, SOp=False):
    benchmark(matrix_type, "canonical to PTM (old)", MAX_NUM_OF_QUBITS["canonical to PTM (old)"], SOp, verbose=True)
    benchmark(matrix_type, "canonical to PTM (iterative)", MAX_NUM_OF_QUBITS["canonical to PTM (iterative)"], SOp, verbose=True)
    benchmark(matrix_type, "canonical to PTM", MAX_NUM_OF_QUBITS["canonical to PTM"], SOp, verbose=True)

def benchmarkPlot(matrix_type, map_type, max_num_of_qubits, SOp=False):
    # fig, ax = plt.subplots()
    # ax.set_ylim([10 ** (-6), 10 ** 2])
    times = benchmark(matrix_type, map_type, max_num_of_qubits, SOp, verbose=True)
    # ax.scatter(range(2, max_num_of_qubits + 1), times[1:])
    # ax.set_yscale('log')
    # ax.set_title(f"PTM {map_type} {matrix_type}")
    if SOp and map_type in MAPSQK.keys():
        timesQK = benchmark(matrix_type, map_type, max_num_of_qubits, SOp, True, verbose=True)
        # ax.scatter(range(2, max_num_of_qubits + 1), timesQK)
        # ax.set_ylim([10 ** (-6), 10 ** 3])
    # plt.savefig(f"PlotsPTM/PTM{map_type}{matrix_type}.png", dpi=150)


if __name__ == "__main__":
    matrix_types = ["random dense"]
    # "random diagonal"
    # map_types = []
    #	 ["Commutator","Anticommutator","Multiplication","Conjugation","Slicing"]
    # SOpmapTypes = ["canonical to PTM (parallel)"]  # SuperOp_new, "Chi"
    for mat_type in matrix_types:
        benchmark_can_to_ptm(mat_type, True)
    #     print(f"\033[4m{mat_type}\033[0m")
    #     for map_type in map_types:
    #         print(f"|-> \033[4m{map_type}\033[0m")
    #         benchmarkPlot(mat_type, map_type, MAX_NUM_OF_QUBITS[map_type])
    #     for map_type in SOpmapTypes:
    #         print(f"|-> \033[4m{map_type}\033[0m")
    #         benchmarkPlot(mat_type, map_type, MAX_NUM_OF_QUBITS[map_type], True)
