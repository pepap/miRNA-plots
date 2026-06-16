# miRNA-plots
MA-plots, CrossHair-plots, fidelity-heatmaps and fidelity-boxplots for smallRNA-seq data

## <ins>miRNAs' MA-plot</ins>
The expected input is the result table of the differential expression analysis from [DESeq2](https://www.bioconductor.org/packages//2.13/bioc/html/DESeq2.html) library. The output table is merged with the [miRBase annotation table]().
The [miRBase 22.1.](https://www.mirbase.org/download/mmu.gff3) set of miRNAs was used for the annotation. The annotation of the mirtrons was taken from [Ladewig et al.](https://doi.org/10.1101/gr.133553.111).
