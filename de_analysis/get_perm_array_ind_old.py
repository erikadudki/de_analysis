import math
import numpy as np
from itertools import combinations as comb
import itertools
# from de_analysis import *
# import de_analysis


def get_perm_array_ind(num_control, num_copd, modus = 'usual'):
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
        modus: 'usual'|'compare_clusters'
            'usual': compare group of patients condition1 vs different
                group of patients condition2
            'compare_clusters': compare same group of patients but with
                different cluster/celltype  annotations
                (e.g. [Pat1_cluster1,Pat2_cluster1,Pat3_cluster1] vs.
                [Pat1_cluster2,Pat2_cluster2,Pat3_cluster2])

    Returns:
        a_ind: nd.array (5005,15)
            for each permutation (total 5005) new ordering of indices (for
            patients)
        num_perm: int
            number of  permutations

    """
    if modus == 'usual':
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

    elif modus == 'compare_clusters':
        # Comparison of same patients in different clusters ( 1_cl1 vs 1_cl2,...)
        # here: num_control = num_copd
        pat_combination = range(0, num_control + num_copd)
        all_list = []
        for i_pat in range(0, num_control):
            all_list.append([i_pat, i_pat + num_control])
            # all_list = [[1, 12], [2, 22], [3,32],[4,42]]

        # to compute all possible permutations
        all_poss_perm = list(itertools.product(*all_list))
        # get only first half of permutations, other half is mirror of the ones before
        num_perm = int(len(all_poss_perm) / 2)
        all_need_perm = all_poss_perm[0:num_perm]

        # fill permutation indices for one group side into list with all patient indices
        a_ind = np.zeros([num_perm, num_control + num_copd], dtype=int)
        for li in range(num_perm):
            a_ind[li, 0:num_control] = all_need_perm[li]

    # fill the remaining columns (for copd group) with remaining indices,
    # which werent used in the control permutation
    for i_row in range(0, num_perm):
        c = num_control
        for item in pat_combination:
            if not item in a_ind[i_row]:
                a_ind[i_row, c] = item
                c = c + 1
    # remove first and last entries which are e.g. ( 0123,4567 ) & (4567,0123)
    col_to_keep = np.shape(a_ind)[0]-1
    a_sub_ind = a_ind[1:col_to_keep,:]
    num_perm = num_perm - 2

    return a_sub_ind, num_perm

#
# num_control=8
# a_ind, num_perm =get_perm_array_ind(num_control,num_control, modus = 'compare_clusters')
# print(num_perm)
# print(2**num_control/2)