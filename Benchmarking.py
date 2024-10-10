import timeit
import matplotlib.pyplot as plt
from Commutator import commutator_to_ptm, anticommutator_to_ptm
from Multiplication import multiplication_to_ptm
from Conjugation import sandwich_to_ptm
from Chi import chi_to_ptm
from SuperOp import can_to_ptm, can_to_ptm_old
from Choi import choi_to_ptm
from PTM_utils import matrix_slice
import test_matrices
import qiskit.quantum_info.operators.channel as ch  # Choi, Chi, PTM, SuperOp

maps = {"Commutator": commutator_to_ptm,
        "Anticommutator": anticommutator_to_ptm,
        "Multiplication": multiplication_to_ptm,
        "Conjugation": sandwich_to_ptm,
        "Slicing": matrix_slice,
        "Chi": chi_to_ptm,
        "SuperOp_old": can_to_ptm_old,
        "SuperOp": can_to_ptm,
        "Choi": choi_to_ptm}

mapsQK = {"Chi": ch.Chi,
          "Choi": ch.Choi,
          "SuperOp": ch.SuperOp}

matrices = {"random diagonal": test_matrices.rand_diag_mat,
            "random dense": test_matrices.rand_dense_mat,
            "identity": test_matrices.identity_mat}

max_num_of_qubits = {"Commutator": 7,
                     "Anticommutator": 7,
                     "Multiplication": 7,
                     "Conjugation": 5,
                     "Slicing": 14,
                     "Chi": 7,
                     "SuperOp_old": 4,
                     "SuperOp": 4,
                     "Choi": 6}


def benchmark(matrix_type, map_type, max_num_of_qubits, SOp=False, QK=False):
    output = True
    times = []
    base = 2
    if SOp:
        base = 4
    if QK:
        print("|   |-> \033[4mQiskit\033[0m")
        for num_of_qubits in range(1, max_num_of_qubits + 1):
            data = matrices[matrix_type](base ** num_of_qubits)
            data = mapsQK[map_type](data)
            start_time = timeit.default_timer()
            ch.PTM(data)
            times.append(timeit.default_timer() - start_time)
            if output and num_of_qubits > 1:
                print(f"|   |   {num_of_qubits}: {times[-1]} (x{times[-1] / times[-2]})")
    else:
        print("|   |-> \033[4mTPD-based\033[0m")
        for num_of_qubits in range(0, max_num_of_qubits + 1):
            data = matrices[matrix_type](base ** num_of_qubits)
            start_time = timeit.default_timer()
            maps[map_type](data)
            times.append(timeit.default_timer() - start_time)
            if output and num_of_qubits > 1:
                print(f"|   |   {num_of_qubits}: {times[-1]} (x{times[-1] / times[-2]})")
    return times


def benchmarkPlot(matrix_type, map_type, max_num_of_qubits, SOp=False):
    fig, ax = plt.subplots()
    ax.set_ylim([10 ** (-6), 10 ** 2])
    if SOp and map_type in mapsQK.keys():
        timesQK = benchmark(matrix_type, map_type, max_num_of_qubits, SOp, True)
        ax.scatter(range(1, max_num_of_qubits + 1), timesQK)
        ax.set_ylim([10 ** (-6), 10 ** 3])
    times = benchmark(matrix_type, map_type, max_num_of_qubits, SOp)
    ax.scatter(range(1, max_num_of_qubits + 1), times[1:])
    ax.set_yscale('log')
    ax.set_title(f"PTM {map_type} {matrix_type}")
    plt.savefig(f"PlotsPTM/PTM{map_type}{matrix_type}.png", dpi=150)


if __name__ == "__main__":
    matrix_types = ["random diagonal", "random dense"]
    map_types = []
    #	 ["Commutator","Anticommutator","Multiplication","Conjugation","Slicing"]
    SOpmapTypes = ["SuperOp"]  # SuperOp_old, "Chi"
    for mat_type in matrix_types:
        print(f"\033[4m{mat_type}\033[0m")
        for map_type in map_types:
            print(f"|-> \033[4m{map_type}\033[0m")
            benchmarkPlot(mat_type, map_type, max_num_of_qubits[map_type])
        for map_type in SOpmapTypes:
            print(f"|-> \033[4m{map_type}\033[0m")
            benchmarkPlot(mat_type, map_type, max_num_of_qubits[map_type], True)
