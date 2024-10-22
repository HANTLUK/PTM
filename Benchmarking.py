import timeit
import time
import numpy as np

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
        "M": ("sandwich multiplication to PTM", sandwich_to_ptm, None, 6),
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
    times: list = []
    mul: list = []
    base: int = 4 if MAPS[map_type][2] else 2
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
            if num_of_qubits > 1:
                mul.append(times[-1]/times[-2])
        else:
            start_time = time.perf_counter()
            MAPS[map_type][1](data)
            times.append(time.perf_counter() - start_time)
            if num_of_qubits > 1:
                mul.append(times[-1]/times[-2])
        if verbose:
            str_multiplicator: str = "" if num_of_qubits == 1 else f" (x{times[-1] / times[-2]})"
            print(f"|   |   {num_of_qubits}: {times[-1]}{str_multiplicator}")
    avg_mul: float = (times[-1]/times[0])**(1./(MAPS[map_type][3]-1))
    max_mul: float = max(mul)
    act_mul: float = mul[-1]
    return [times, [avg_mul, max_mul, act_mul]]


def benchmark_plot(matrix_type: str, map_type: str):
    name_file: str = f"{map_type}-PTM_{matrix_type}"
    fig, ax = plt.subplots()
    ax.set_ylim([10 ** (-6), 10 ** 4])
    times, mul = benchmark(matrix_type, map_type, verbose=True)
    avg_mul, max_mul, act_mul = mul
    ran_ex: list = range(len(times),len(times)+3)
    ex_times: np.ndarray = np.array([times[-1]*avg_mul**(i-len(times)) for i in ran_ex])
    max_ex_times: np.ndarray = np.array([times[-1]*max_mul**(i-len(times)) for i in ran_ex])
    act_ex_times: np.ndarray = np.array([times[-1]*act_mul**(i-len(times)) for i in ran_ex])
    with open(f"TabsPTM/{name_file}_TPD.dat", "w") as myfile:
        out_string: str = "".join([f"{i+1}\t {time:.5f}\n" for i,time in enumerate(times)])
        myfile.write(f"{out_string}")
    with open(f"TabsPTM/{name_file}_TPD_ex.dat", "w") as myfile:
        out_string: str = "".join([f"{i+len(times)}\t {time:.5f} \t {max_time:.5f} \t {act_time:.5f} \n" for i,(time,max_time,act_time) in enumerate(zip(ex_times,max_ex_times,act_ex_times))])
        myfile.write(f"{out_string}")
    ax.scatter(range(1, MAPS[map_type][3] + 1), times)
    ax.set_yscale('log')
    ax.set_title(f"{MAPS[map_type][0]}: {MATRICES[matrix_type][0]}")
    if MAPS[map_type][2]:
        times_qiskit, mul_qiskit = benchmark(matrix_type, map_type, True, verbose=True)
        ax.scatter(range(1, MAPS[map_type][3] + 1), times_qiskit)
        avg_mul, max_mul, act_mul = mul_qiskit
        ran_ex: list = range(len(times_qiskit),len(times_qiskit)+3)
        ex_times: np.ndarray = np.array([times_qiskit[-1]*avg_mul**(i-len(times_qiskit)) for i in ran_ex])
        max_ex_times: np.ndarray = np.array([times_qiskit[-1]*max_mul**(i-len(times_qiskit)) for i in ran_ex])
        act_ex_times: np.ndarray = np.array([times_qiskit[-1]*act_mul**(i-len(times_qiskit)) for i in ran_ex])
        ax.set_ylim([10 ** (-6), 10 ** 3])
        with open(f"TabsPTM/{name_file}_QK.dat", "w") as myfile:
            out_string: str = "".join([f"{i+1}\t {time:.5f}\n" for i,time in enumerate(times_qiskit)])
            myfile.write(f"{out_string}") 
        with open(f"TabsPTM/{name_file}_QK_ex.dat", "w") as myfile:
            out_string: str = "".join([f"{i+len(times_qiskit)}\t {time:.5f} \t {max_time:.5f} \t {act_time:.5f} \n" for i,(time,max_time,act_time) in enumerate(zip(ex_times,max_ex_times,act_ex_times))])
            myfile.write(f"{out_string}") 
    plt.savefig(f"PlotsPTM/{name_file}.png", dpi=150)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    MATRIX_TYPES = ["randDense", "randDiag"]
    MAP_TYPES = ["M"]
    for MAT_TYPE in MATRIX_TYPES:
        print(f"\033[4m{MATRICES[MAT_TYPE][0]}\033[0m")
        for MAP_TYPE in MAP_TYPES:
            print(f"|-> \033[4m{MAPS[MAP_TYPE][0]}\033[0m")
            benchmark_plot(MAT_TYPE, MAP_TYPE)
