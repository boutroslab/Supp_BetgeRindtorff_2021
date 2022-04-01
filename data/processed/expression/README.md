### How was gene expression data processed within the promise project? How was it fed into MOFA?
The processing of the gene expression data within the promise project was a multi-step process: 
1) supp_BetgeRindtorff_2021/code/data/make_expression.R to read all microarray raw data, annotate probes and perform RMA pre-processing. The Current script takes two command line arguments that are specific to the promise project and should not be relevant for follow-up studies. 
2) supp_BetgeRindtorff_2021/code/data/tidy_expression.R (this is probably where you got stuck, because I just realised the file was not pushed because of an issue with .gitignore) The script turns the promise_expr object through a series of transformations including a) adding more metadata, b) removing corrupted samples and c) calculating the coefficient of variation and keeping the genes within the top 10% of varying genes. 
3) supp_BetgeRindtorff_2021/code/models/mofa/tidy_mofa.R is the script that takes the filtered and top 10% varying gene expression data and formats it into a MOFA compatible format. It then feeds it into the MOFA model, leaving the scale_input argument at TRUE (default for MOFA). 

### Why did you filter the genes by their coefficient of variation? 
The coefficient of variation is the standard deviation of a distribution scaled by its mean. This allows us to select for genes that *vary* among samples a lot while accounting for the fact that some genes have different average expression levels. For example, a ribosomal RNA sample is more highly expressed than a differentiation marker. By dividing the standard deviation with the average expression value, we can try to account for these differences. 
Generally, MOFA works best if the various data matrixes fed into the model have about the same dimensionality. To reduce the number of features for the gene expression data dimension, it was necessary to select a subset of genes. 

