import math
import numpy as np
from itertools import combinations as comb


def get_perm_array_ind(num_control, num_copd):
    """
    Get all (here:5005) combinations of possible permutations, build array,
    in each row: new permutation
    e.g., indices 0,..5: control;       6,...,14: copd patients
    permute indices
    num_perm = 5005

    Args:
        num_control: int
            number of patients in Group 1
        num_copd: int
            number of patients in Group 2

    Returns:
        a_ind: nd.array (5005,15)
            for each permutation (total 5005) new ordering of indices (for
            patients)

    """

    num_perm = math.factorial(num_control + num_copd) / (
            math.factorial(num_control) * math.factorial(num_copd))
    num_perm = int(num_perm)

    li = 0
    a_ind = np.zeros([num_perm, num_control + num_copd], dtype=int)
    pat_combination = range(0, num_control + num_copd)  # [0,1,2,3,....,14]
    # from pat_combination (all patient indices) get 6 permuted values
    comb1 = comb(pat_combination, num_control)
    for i in list(comb1):
        # print(i)
        a_ind[li, 0:num_control] = np.asarray(i)  # i:tuple
        li = li + 1

    # fill the remaining columns (for copd group) with remaining indices,
    # which werent used in the control permutation
    for i_row in range(0, num_perm):
        c = num_control
        for item in pat_combination:
            if not item in a_ind[i_row]:
                a_ind[i_row, c] = item
                c = c + 1
    return a_ind