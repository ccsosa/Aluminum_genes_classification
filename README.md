# Aluminum_genes_classification
Research schema used: 

![Code schema](https://raw.githubusercontent.com/ccsosa/Aluminum_genes_classification/master/DIAGRAM.png)



- First use the interactome downloaded from Mishra et al., (2022) in edge list format
- Add a count matrix and run the GRAPH_AD2_OMICAS.R script to add the correlations of log2(X+1) to the edges
- Run PECAN.sh to obtain an embedded graph using Pecanpy Python module and a undirected weighted graph using node2vec+ algorithm
- Run PULEARN_ALUMINUM.py with the obtained embedded graph. This code runs support vector machine, Random forest and linear discriminant algorithms using PULEARN with 10000 estimators. 
> This code helps to obtain fixed metrics for PULEARN since the regular metrics do not work well with PULEARNING

- Please check summary files to obtain a brief summary of results
- Please check results folder to obtain detailed files including prediction per gene
