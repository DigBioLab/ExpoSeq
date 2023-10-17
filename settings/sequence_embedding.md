### Clustering based on t-SNE

Addressing the limitations of using the levenshtein distance for clustering hcdr3 sequences, a plot which shows the sequences and their whole sequence relation was created to enable a more extensive analysis of the sequences and their binding. To be able to have the sequences as vectors the SGT embedding was applied.

This embedding does not capture any chemical properties of the amino acid strand but instead tries to recognize the characteristic relative position of letters within a sequence which enables an identification of patterns between sequences of different lengths. This addresses the variability of the cdr3 sequences and enables to capture the similarity of sequences. Therefore, the sgt embedder was initialized using the package sgt where the parameter length sensitive was set to True and the parameter kappa was assigned with 1 . After the initialization, the embedding creates  an output of 400 dimensions per sequence which were reduced to 80 dimensions by using principal component analysis (PCA).

Subsequently t-SNE was implemented to show a representative arrangement of the 80 dimensions in the two dimensional space because a further dimension reduction to only two dimensions would insufficiently describe the data globally by using PCA. That is why t-SNE was introduced to be able to show a representative arrangement of the points without losing too much of its information. 

The plot can be used in variety of ways, such as: 
1. Identify similar CDR3 regions
2. learn about the differences between samples
3. identify specific sequence differences for samples which were panned against different antigens.