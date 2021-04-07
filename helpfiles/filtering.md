## Filtering for low expression

***

Prior to calculation of the haemolysis metric, the difference between the geometric mean of signature and background miRNA,
strict filtering is performed. This filtering helps to ensure miRNA with very low counts are removed.
For a miRNA to be considered biologically relevant it must be expressed at some minimun level.
Here filtering is modeled on the gene expression filtering recommended in the _edgeR_ [Users Guide](http://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf).

***

Briefly, this number should be equal to the number of samples in your smallest group.
For example, if your experiment contains 3 samples from control and 4 samples from a treatment you select 3 as the number of samples to use for filtering. 

