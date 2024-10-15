import timeit
import time

import logging
logger = logging.getLogger(__name__)

import matplotlib.pyplot as plt
from Commutator import commutator_to_ptm, anticommutator_to_ptm
from Multiplication import multiplication_to_ptm
from Conjugation import sandwich_to_ptm
from IterChi import chi_to_ptm
from SuperOp import can_to_ptm
from Choi import choi_to_ptm
from PTM_utils import matrix_slice
import test_matrices
import qiskit.quantum_info.operators.channel as ch  # Choi, Chi, PTM, SuperOp

MAPS = {"C": ("commutator to PTM", commutator_to_ptm, None, 7),
        "AC": ("anticommutator to PTM", anticommutator_to_ptm, None, 7),
        "L": ("left multiplication to PTM", multiplication_to_ptm, None, 7),
        "M": ("sandwich multiplication to PTM", sandwich_to_ptm, None, 7),
        "Slicing": ("matrix slicing", matrix_slice, None, 7),
        "Can": ("canonical to PTM", can_to_ptm, ch.SuperOp, 7),
        "Chi": ("Chi to PTM", chi_to_ptm, ch.Chi, 7),
        "Choi": ("Choi to PTM", choi_to_ptm, ch.Choi, 7)
        }

MATRICES = {"randDiag": ("random diagonal", test_matrices.rand_diag_mat),
            "randDense": ("random dense", test_matrices.rand_dense_mat),
            "id": ("identity", test_matrices.identity_mat)
            }


def benchmark(matrix_type: str, map_type: str, from_qiskit=False, verbose=False):
    times = []
    base = 4 if MAPS[map_type][2] else 2
    str_from_qiskit: str = "Qiskit" if from_qiskit else "TPD-based"
    print(f"|   |-> \033[4m {str_from_qiskit}\033[0m")

    warm_up_data = MATRICES[matrix_type][1](base ** 2)
    if from_qiskit and MAPS[map_type][2]:
        warm_up_data = MAPS[map_type][2](warm_up_data)
        ch.PTM(warm_up_data)
    else:
        MAPS[map_type][1](warm_up_data)

    for num_of_qubits in range(1, MAPS[map_type][3] + 1):
        data = MATRICES[matrix_type][1](base ** num_of_qubits)
        if from_qiskit and MAPS[map_type][2]:
            start_time = time.perf_counter()
            data = MAPS[map_type][2](data)
            ch.PTM(data)
            times.append(time.perf_counter() - start_time)
        else:
            start_time = time.perf_counter()
            MAPS[map_type][1](data)
            times.append(time.perf_counter() - start_time)
        if verbose:
            str_multiplicator: str = "" if num_of_qubits == 1 else f" (x{times[-1] / times[-2]})"
            print(f"|   |   {num_of_qubits}: {times[-1]}{str_multiplicator}")
    return times


def benchmark_plot(matrix_type: str, map_type: str):
    name_file: str = f"{map_type}-PTM_{matrix_type}"
    fig, ax = plt.subplots()
    ax.set_ylim([10 ** (-6), 10 ** 2])
    times = benchmark(matrix_type, map_type, verbose=True)
    with open(f"TabsPTM/{name_file}.dat", "w") as myfile:
        out_string: str = "".join([f"{i+1}\t {time:.5f}\n" for i,time in enumerate(times)])
        myfile.write(f"{out_string}")
    ax.scatter(range(1, MAPS[map_type][3] + 1), times)
    ax.set_yscale('log')
    ax.set_title(f"{MAPS[map_type][0]}: {MATRICES[matrix_type][0]}")
    if MAPS[map_type][2]:
        times_qiskit = benchmark(matrix_type, map_type, True, verbose=True)
        ax.scatter(range(1, MAPS[map_type][3] + 1), times_qiskit)
        ax.set_ylim([10 ** (-6), 10 ** 3])
        with open(f"TabsPTM/{name_file}.dat", "a") as myfile:
            out_string: str = "".join([f"{i+1}\t {time:.5f}\n" for i,time in enumerate(times)])
            myfile.write(f"{out_string}") 
    plt.savefig(f"PlotsPTM/{name_file}.png", dpi=150)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    MATRIX_TYPES = ["randDense"]
    MAP_TYPES = ["Can"]
    for MAT_TYPE in MATRIX_TYPES:
        print(f"\033[4m{MATRICES[MAT_TYPE][0]}\033[0m")
        for MAP_TYPE in MAP_TYPES:
            print(f"|-> \033[4m{MAPS[MAP_TYPE][0]}\033[0m")
            benchmark_plot(MAT_TYPE, MAP_TYPE)
