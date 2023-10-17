### Heatmap based on Morosita-Horn Index

The heatmap shows a matrix which values are calculated based on the morosita horn index. This index captures the degree of identity between two samples which are composed of multiple values of different sizes. This is especially useful for comparing two samples which are composed of different numbers of clones and have a different number of aligned reads. In the context of this heatmap, identity means that two sequences have the same sequence length and the exact same amino acid sequence. 

The heatmap has several possible applications during the quality control of the sequencing, which are:
1. identification of cross contamination between samples. This can be seen if a high degree of identity can be observed, although the samples were panned against different antigens.
2. Validation of the panning process. This can be seen if a high degree of identity can be observed between samples which were panned against the same antigen in different rounds.


*Weaknesses*:
- the morosita horn becomes less accuracte if the number of elements (sequences) is very low in one of the samples

*Alternatives*:

The user can analyse the identity between the samples based on:
- Jaccard Index
- Sorensen-Dice Index