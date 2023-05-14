# Aluminum_genes_classification

- First use the interactome downloaded from Mishra et al., (2022) in edge list format
- Add a count matrix and run the GRAPH_AD2_OMICAS.R script to add the correlations of log2(X+1) to the edges
- Run PECAN.sh to obtain an embedded graph using Pecanpy
- Run PULEARN_ALUMINUM.py with the obtained embedded graph. This code runs support vector machine, Random forest and linear discriminant analyisis using PULEARN. 
> This code helps to obtain fixed metric for PULEARN since the regular metrics do not work well with PULEARNING
