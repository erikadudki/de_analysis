import math
import numpy as np
from itertools import combinations as comb
import itertools
# from de_analysis import *
# import de_analysis


def get_perm_array_ind(num_control, num_copd, modus = 'towsided'):
    """
    Get all (here:5005) combinations of possible permutations, build array,
    in each row: new permutation
    e.g., indices 0,..5: control;       6,...,14: copd patients
    permute indices
    remove constellation of main group ([0 1 2 3 4 5])
    num_perm = 5005

    Args:
        num_control: int
            number of patients in Group 1
        num_copd: int
            number of patients in Group 2
        modus: 'onesided'|'twosided'|'compare_clusters'
            'onesided': if number of patients in both groups are the same,
                the index-permutations are cut in half, because of symmetry of
                permutations (e.g. main index:[0 1 2 3]; permutations [0 2 1 3],
                [0 3 1 2],[1 2 3 0],[1 3 2 0], last two are repetition)
            'twosided': get all indices of all permutations, also mirrored
                indices
            'compare_clusters': compare same group of patients but with
                different cluster/celltype  annotations
                Comparison of same patients in different clusters ( 1_cl1 vs 1_cl2,...)
                (e.g. [Pat1_cluster1,Pat2_cluster1,Pat3_cluster1] vs.
                [Pat1_cluster2,Pat2_cluster2,Pat3_cluster2])

    Returns:
        a_ind: nd.array (5005,15)
            for each permutation new ordering of indices (for
            patients)
        num_perm: int
            number of  permutations

    """
    if (modus == 'onesided') or (modus == 'twosided'):
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
        # here: num_control = num_copd !!!
        if num_control == num_copd:
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
        else: #if  num_control != num_copd
            # TODO
            pat_combination = range(0, num_control + num_copd)
            all_list = []

            # construct the permutations first for the smaller group
            min_num_g1_g2 = min(num_control, num_copd)
            for i_pat in range(0, min_num_g1_g2):
                all_list.append([i_pat, i_pat + min_num_g1_g2])

            # for i_pat in range(0, num_control):
            #     all_list.append([i_pat, i_pat + num_control])
                # all_list = [[1, 12], [2, 22], [3,32],[4,42]]

            # to compute all possible permutations
            all_poss_perm = list(itertools.product(*all_list))
            # get only first half of permutations, other half is mirror of the ones before
            num_perm = len(all_poss_perm)  # int(len(all_poss_perm) / 2)
            all_need_perm = all_poss_perm#[0:num_perm]

            # fill permutation indices for one group side into list with all patient indices
            a_ind = np.zeros([num_perm, num_control + num_copd], dtype=int)
            for li in range(num_perm):
                a_ind[li, 0:min_num_g1_g2] = all_need_perm[li]



    else:
        raise ValueError('Given permutation modus not available. Choose between: '
                         'onesided, twosided or compare_clusters.')

    # fill the remaining columns (for copd group) with remaining indices,
    # which werent used in the control permutation
    for i_row in range(0, num_perm):
        if num_control != num_copd:
            c = min_num_g1_g2
        else:
            c = num_control
        for item in pat_combination:
            if not item in a_ind[i_row]:
                a_ind[i_row, c] = item
                c = c + 1
    # flip the index-results, if group2 is smaller then group1, such that
    # order of indexes is correct (permuted indices group1, then group2)
    if num_copd < num_control:
        a_ind = np.fliplr(a_ind)
    # remove first and last entries which are e.g. ( 0123,4567 ) & (4567,0123)
    # if patient groups same number of patients
    if num_control == num_copd:
        # remove only first entry which are e.g. remove( 012,345 ), last entry (045,123) is needed!
        # if patient groups same number of patients
        if modus == 'compare_clusters':
            col_to_keep = np.shape(a_ind)[0]    #keep until all entries
            a_sub_ind = a_ind[1:col_to_keep, :] #remove first entry
            num_perm = num_perm - 1
        else:
            # remove first and last entries which are e.g. ( 0123,4567 ) & (4567,0123)
            # if patient groups same number of patients
            col_to_keep = np.shape(a_ind)[0]-1
            a_sub_ind = a_ind[1:col_to_keep,:]
            num_perm = num_perm - 2
    else:
        # remove first row (first constellation of main group) (123-4567)
        a_sub_ind = a_ind[1:num_perm,:]
        num_perm = num_perm - 1

    if modus == 'onesided':
        if num_control == num_copd:
            # get only first half of permutations, other half is mirror of the ones before
            num_perm = int(num_perm/2)
            a_sub_ind = a_sub_ind[0:num_perm]

    return a_sub_ind, num_perm

#
# num_control=2
# a_ind, num_perm =  get_perm_array_ind(4,4, modus = 'compare_clusters')
# print(num_perm)
# print(2**num_control/2)
# print(a_ind)
