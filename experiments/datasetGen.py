from fpylll import FPLLL, GSO, LLL, IntegerMatrix, Enumeration, EnumerationError, Pruning
import csv
import random
from multiprocessing import Pool
import sys

seed = 2021

LLL_reduce = True

# Number of CPUs available
# Number of instances should be multiple of number of cpus!
cpus = 1  # 64
num_instances = 2  # 128

if len(sys.argv) != 3 or sys.argv[1].isnumeric() == False or sys.argv[2].isnumeric() == False:
    print("Rank limits not provided. The correct usage is 'python datasetGen.py rank_min rank_max'")
    sys.exit()

ranks = range(int(sys.argv[1]), int(sys.argv[2])+1)
q = 65537
full_dim = 180

out_name = "qary_" + str(ranks[0]) + "_" + str(ranks[-1]) + ".csv"

"Runs classical enumeration and return the found vector"
def run_enum(M, pruned=False, pruned_dist=0):
    enum_obj = Enumeration(M)
    dist, exp = M.get_r_exp(0, 0)

    try:
        if pruned:
            dist, v = enum_obj.enumerate(0, M.d, pruned_dist, 1)[0]
        else:
            dist, v = enum_obj.enumerate(0, M.d, dist, exp)[0]  # , pruning=pruning.coefficients)[0]

    except EnumerationError:
        return None, None

    return dist, v


def generate_basis(d, q, _ranks):
    A = IntegerMatrix.random(d, "qary", k=d // 2, q=q)
    if LLL_reduce:
        return LLL.reduction(A)[:_ranks[-1]]
    return A[:_ranks[-1]]


def instance_calc(instance_seed):
    FPLLL.set_random_seed(instance_seed)
    basis = generate_basis(full_dim, q, ranks)
    dists = []
    vs = []
    j = 0
    for rank in ranks:
        A = basis[:rank]
        M = GSO.Mat(A)
        M.update_gso()

        if j == 0:
            dist, v = run_enum(M)
        else:
            dist, v = run_enum(M, True, dist)

        dists.append(dist)
        vs.append(v)
        j += 1
    return dists, vs, basis


if __name__ == "__main__":
    random.seed(seed)

    pool = Pool(cpus)

    if num_instances % cpus != 0:
        print("ERROR. Number of instances should be multiple of number of cpus")

    with open(out_name[:-4] + "_matrices.csv", 'w', newline='') as matricesfile:
        writer_matrices = csv.writer(matricesfile, delimiter=';')
        with open(out_name, 'w', newline='') as csvfile:
            csvfile.write("num_ranks " + str(ranks[-1]-ranks[0]+1)+"\n")
            csvfile.write("rank_min " + str(ranks[0])+"\n")
            csvfile.write("num_instances " + str(num_instances)+"\n")
            csvfile.write("dim " + str(full_dim)+"\n")
            writer = csv.writer(csvfile, delimiter=';')

            for i in range(num_instances // cpus):
                res = pool.map(instance_calc, [random.randint(0, 99999999) for i in range(cpus)])
                for dists, vs, A_reduced in res:
                    writer.writerow(zip(dists, vs))
                    writer_matrices.writerow(A_reduced)
                csvfile.flush()
                matricesfile.flush()
