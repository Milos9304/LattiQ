# General imports
import numpy as np
from numpy import linalg as lg
import matplotlib.pyplot as plt
import glob

# Pre-defined ansatz circuit, operator class and visualization tools
from qiskit.circuit.library import QAOAAnsatz
from qiskit_algorithms.minimum_eigensolvers.diagonal_estimator import _DiagonalEstimator

from qiskit.quantum_info import SparsePauliOp
from qiskit.visualization import plot_distribution

from qiskit_algorithms.minimum_eigensolvers import QAOA
from qiskit_algorithms.optimizers import COBYLA
from qiskit.primitives import Sampler
from qiskit.result import QuasiDistribution
from qiskit.quantum_info import Pauli, Statevector
from qiskit_algorithms.utils import validate_initial_point, validate_bounds



from qiskit_algorithms.utils.set_batching import _set_default_batchsize

#from qiskit.primitives import Estimator
from qiskit_aer.primitives import Estimator # import change!!!

# SciPy minimizer routine
from scipy.optimize import minimize
from fileinput import filename
import math


#%config InlineBackend.figure_format='retina'
def bitfield(n: int, L: int) -> list[int]:
    result = np.binary_repr(n, L)
    return [int(digit) for digit in result]  # [2:] to chop off the "0b" part
# Problem to Hamiltonian operator
def objective_value(x: np.ndarray, g: np.ndarray) -> float:
    print(x.shape)
    print(g.shape)
    return np.matmul(np.matmul(x.transpose(),g),x)
def _compare_measurements(candidate, current_best):
    """Compare two best measurements. Returns True if the candidate is better than current value.

    This compares the following two criteria, in this precedence:

        1. The smaller objective value is better
        2. The higher probability for the objective value is better

    """
    if candidate["value"] < current_best["value"]:
        return True
    elif candidate["value"] == current_best["value"]:
        return candidate["probability"] > current_best["probability"]
    return False
def sample_most_likely(state_vector: QuasiDistribution | Statevector) -> np.ndarray:
    """Compute the most likely binary string from state vector.
    Args:
        state_vector: State vector or quasi-distribution.

    Returns:
        Binary string as an array of ints.
    """
    if isinstance(state_vector, QuasiDistribution):
        values = list(state_vector.values())
    else:
        values = state_vector
    n = int(np.log2(len(values)))
    k = np.argmax(np.abs(values))
    x = bitfield(k, n)
    x.reverse()
    return np.asarray(x)
def cost_func(params, ansatz, hamiltonian, estimator):
    """Return estimate of energy from estimator

    Parameters:
        params (ndarray): Array of ansatz parameters
        ansatz (QuantumCircuit): Parameterized ansatz circuit
        hamiltonian (SparsePauliOp): Operator representation of Hamiltonian
        estimator (Estimator): Estimator primitive instance

    Returns:
        float: Energy estimate
    """
    cost = estimator.run(ansatz, hamiltonian, parameter_values=params, shots=1000).result().values[0]
    return cost
def store_best_measurement(best):
    for best_i in best:
        if best_measurement["best"] is None or _compare_measurements(
            best_i, best_measurement["best"]
        ):
            best_measurement["best"] = best_i
files=glob.glob("*.ising")
for file_name in files:
    #print(file_name)
    f=open(file_name)
    fg=open(file_name[0:-5]+"gram")
    G=fg.readlines()
    G=np.array(list(map(lambda y:list(map(int,y.split())),map(lambda x: x[:-1] if x[-1] == '\n' else x,G))))
    
    hamil_list=[]
    for line in f.readlines():
        c, s = line.split(" ")
        if s[-1] == '\n':
            s=s[0:-1]
        hamil_list.append((s,c))

hamiltonian = SparsePauliOp.from_list(hamil_list)
evals, evects = lg.eig(hamiltonian)
minimum_indices = np.where(evals == min(evals))
#Lambda = np.diag(Eigenvalues)
print(evals)
#print(Eigenvectors)
#print(Lambda)
#print(hamiltonian)
global p
p=1

# QAOA ansatz circuit
#ansatz = QAOAAnsatz(hamiltonian, reps=2)
# Draw
#ansatz.decompose(reps=3).draw("mpl", filename="f")

#estimator = Estimator(run_options= {"method": "statevector"})
#sampler = Sampler(session=session, options={"shots": int(1e4)})

solutions=list(map(lambda x: bin(x)[2:], list(minimum_indices[0]))) #["0010","1001"]
for i in range(len(solutions)):
    if(len(solutions[i]) < hamiltonian.num_qubits):
        solutions[i] = "0"*(hamiltonian.num_qubits-len(solutions[i]))+solutions[i]
num_samples = 10

#hamiltonian_m=hamiltonian.to_matrix()
#v=np.zeros(16)
#v[9]=1
#print(np.matmul(hamiltonian_m, v))

#x0 = [1.02207,-0.313715,1.02207,-0.313715,1.02207,-0.313715,1.02207,-0.313715,1.02207,-0.313715,1.02207,-0.313715,1.02207,-0.313715,1.02207,-0.313715]#, 0.0193662, 2.10431]#2 * np.pi * np.random.rand(ansatz.num_parameters)
final_stats = []

np.random.seed(0)
for sample_i in range(num_samples):
    
    x0=2 * np.pi * np.random.rand(p*2) - np.pi  
    sampler = Sampler()
    optimizer = COBYLA()
    ansatz = QAOAAnsatz(hamiltonian, reps=p)#, mixer_operator=None).decompose()
    ansatz.parameter_bounds = [(-math.pi,math.pi)]*(2*p)

    if len(ansatz.clbits) > 0:
        ansatz.remove_final_measurements()
    ansatz.measure_all()
    initial_point = validate_initial_point(x0, ansatz)
    bounds = validate_bounds(ansatz)
    best_measurement = {"best": None}
    estimator = _DiagonalEstimator(
        sampler=sampler, callback=store_best_measurement, aggregation=None)
    def evaluate_energy(parameters: np.ndarray) -> np.ndarray | float:
        global eval_count
        # handle broadcasting: ensure parameters is of shape [array, array, ...]
        parameters = np.reshape(parameters, (-1, 2*p)).tolist()
        batch_size = len(parameters)
        
        estimator_result = estimator.run(
            batch_size * [ansatz], batch_size * [hamiltonian], parameters
        ).result()
        values = estimator_result.values
        
        result = values if len(values) > 1 else values[0]
        return np.real(result)
    eval_count = 0
    was_updated = _set_default_batchsize(optimizer)
    optimizer_result = optimizer.minimize(
                    fun=evaluate_energy, x0=initial_point, bounds=bounds
                )
    if was_updated:
        optimizer.set_max_evals_grouped(None)
    
    print(optimizer_result)
    hit_rate = 0
    final_state = sampler.run([ansatz], [optimizer_result.x]).result().quasi_dists[0]
    new_row=np.zeros(len(final_state.items()))
    for key, value in final_state.items():
        new_row[key] = value
    final_stats.append(new_row)
    
    for w in sorted(final_state, key=final_state.get, reverse=True):
        b=str(bin(w))[2:]
        b="0"*(ansatz.num_qubits-len(b))+b
        if b in solutions:
            hit_rate += final_state[w]
            print(b,":",final_state[w])
    print("HR:", hit_rate)


#plt.bar(range(len(final_state_list)),final_state_list)
bars=plt.bar(range(len(final_stats[0])), np.mean(final_stats, axis=0), yerr=np.std(final_stats, axis=0), ecolor='black', capsize=10)#, fmt='-o')
for sol in list(minimum_indices[0]):
    bars[sol].set_color('red')

for i in range(len(bars)):
    yval = bars[i].get_height()
    plt.text(bars[i].get_x(), yval + .005, int(np.real(evals[i])),fontsize='x-large')
plt.title("p="+str(p)+" num_samples="+str(num_samples))
plt.show()

    #optimizer_result = qaoa.optimizer(fun=evaluate_energy, x0=initial_point, bounds=bounds)
    #result = qaoa.compute_minimum_eigenvalue(hamiltonian)
    #x = sample_most_likely(result.eigenstate)
    #print(final_state)
    #print(bounds)
    #print(initial_point)
    #print(result)
    #print(f'Objective value computed by QAOA is {objective_value(x, G)}')