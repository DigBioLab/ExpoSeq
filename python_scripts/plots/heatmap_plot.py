import numpy as np
import seaborn as sns



heatmap_axis = input("How many experiments do you want to plot? Max. No. {max_experiments}")
heatmap_axis = int(heatmap_axis)
heatmap_absolute = np.zeros([heatmap_axis, heatmap_axis])
heatmap_absolute_jaccard = np.zeros([heatmap_axis, heatmap_axis])
heatmap_absolute_sorensen = np.zeros([heatmap_axis, heatmap_axis])
heatmap_absolute = heatmap_absolute.astype(np.uint32)
heatmap_morosita_horn = np.zeros([heatmap_axis, heatmap_axis])
matches_collected = []
sum_matches = 0

nested_seq = sequencing_report_all.nSeqCDR3

start = 0
index_x = 0
lib_matches = 0
for first_cluster in range(0, 37):
    for second_cluster in range(start, 37):
        for seq in nested_seq[first_cluster]:
            for seq2 in nested_seq[second_cluster]:
                if seq == seq2:
                    lib_matches += 1
                else:
                    pass
        a = lib_matches
        b = len(nested_seq[first_cluster]) - lib_matches
        c = len(nested_seq[second_cluster]) - lib_matches
        heatmap_absolute_jaccard[start, second_cluster] = lib_matches/(lib_matches+b+c)
        heatmap_absolute_jaccard[second_cluster, start] = lib_matches / (lib_matches + b + c)
        heatmap_absolute_sorensen[start, second_cluster] = 2*lib_matches/(2*lib_matches+b+c)
        heatmap_absolute_sorensen[second_cluster, start] = 2 * lib_matches / (2 * lib_matches + b + c)
        heatmap_absolute[start, second_cluster] = lib_matches
        heatmap_absolute[second_cluster, start] = lib_matches

       # index_x += 1
        lib_matches = 0
    start += 1


nested_clone_fraction = sequencing_report_all.cloneFraction
array_clone_counts = np.array(nested_clone_fraction)
nested_seq = sequencing_report_all.nSeqCDR3
#summed_clone_counts = tidy_data.sum_nest_data(nested_column_name = "cloneCount",
                                #              axis = 0)


start = 0
index_x = 0
lib_matches = 0
for first_cluster in range(0, 37):
    for second_cluster in range(start, 37):
        sum_xy = []
        all_xi = []
        all_yi = []
        sequences1 = []
        sequences2 = []
        for seq_index in range(0, len(nested_seq[first_cluster])):
            seq = nested_seq[first_cluster][seq_index]
            for seq2_index in range(0, len(nested_seq[second_cluster])):
                seq2 = nested_seq[second_cluster][seq2_index]
                if seq == seq2 and (seq not in sequences1 or seq2 not in sequences2): # oherwise if seq2 > seq1 they will be counted twice
                    xi = nested_clone_fraction[first_cluster][seq_index]
                    all_xi.append(xi)
                    yi = nested_clone_fraction[second_cluster][seq2_index]
                    all_yi.append(yi)
                    product_xy = xi*yi
                    sum_xy.append(product_xy)
                    sequences1.append(seq)
                    sequences2.append(seq2)
                else:
                    pass
        sum_xy = np.sum(np.array(sum_xy))
        all_xi = np.array(all_xi).T
        all_yi = np.array(all_yi).T
        numerator_morosita_horn = 2*sum_xy
        lambda_x = np.dot(all_xi, all_xi)
        lambda_y = np.dot(all_yi, all_yi)
        denumerator_morosita_horn = (np.sum(lambda_y) + np.sum(lambda_y))
        morosita_horn_index = numerator_morosita_horn/denumerator_morosita_horn
        heatmap_morosita_horn[start, second_cluster] = morosita_horn_index
        heatmap_morosita_horn[second_cluster, start] = morosita_horn_index

    start += 1




## normalizing
clone_counts = sequencing_report_all.cloneCount


x = np.array(list(zip(summed_counts)))
y = np.array(summed_counts)
b = np.broadcast(x, y)
out = np.empty(b.shape)
out.flat = [u+v for (u,v) in b]

normalized_heatmap = heatmap_absolute/(out/2)
sns.heatmap(normalized_heatmap, cmap="Blues", annot=True, annot_kws={"size":4}, fmt='.2f')