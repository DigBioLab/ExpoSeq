cimport numpy as np
import numpy as np
cimport cython


cpdef np.ndarray[double, ndim = 2] loop_morosita_horn(nested_seq,
                                                        nested_clone_fraction):
    cdef int first_cluster, second_cluster, seq_index, index_seq2, start, axis
    cdef double xi, yi, product_xy, denmuerator_morosita_horn, numerator_morosita_horn, morosita_horn_index, sum_xy,
    cdef np.ndarray [double, ndim = 2]heatmap_morosita_horn
    cdef np.ndarray  clone_frac_one, clone_frac_two
    cdef np.ndarray [double, ndim = 1] all_xy, all_xi, all_yi, lambda_x, lambda_y
    cdef str seq, seq2
    cdef list sample_one, sample_two

    axis = nested_clone_fraction.shape[0]
    heatmap_morosita_horn = np.zeros([axis, axis])
    start = 0
    for first_cluster in range(0, axis):
        clone_frac_one = np.array(nested_clone_fraction[first_cluster])
        sample_one = nested_seq[first_cluster]
        for second_cluster in range(start, axis):
            all_xy = np.zeros([1])
            all_xi = np.zeros([1])
            all_yi = np.zeros([1])
            sample_two = nested_seq[second_cluster]
            clone_frac_two = np.array(nested_clone_fraction[second_cluster])
            for seq_index in range(0, len(sample_one)):
                seq = sample_one[seq_index]
                length_seq = len(seq)
                if seq in sample_two: # oherwise if seq2 > seq1 they will be counted twice  and seq not in used_seq
                    index_seq2 = sample_two.index(seq)
                    seq2 = sample_two[index_seq2]
                    if length_seq == len(seq2):
                        xi = clone_frac_one[seq_index]
                        all_xi = np.append(all_xi, xi)
                        yi = clone_frac_two[index_seq2]
                        all_yi = np.append(all_yi, yi)
                        product_xy = xi*yi
                        all_xy = np.append(all_xy, product_xy)
                    else:
                        pass
                else:
                    pass
            sum_xy = np.sum(all_xy)
            numerator_morosita_horn = 2*sum_xy
            lambda_x = all_xi**2
            lambda_y = all_yi**2
            denumerator_morosita_horn = (np.sum(lambda_x) + np.sum(lambda_y))
            morosita_horn_index = numerator_morosita_horn/denumerator_morosita_horn
            heatmap_morosita_horn[start, second_cluster] = morosita_horn_index
            heatmap_morosita_horn[second_cluster, start] = morosita_horn_index
        start += 1
    return heatmap_morosita_horn