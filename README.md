# diffUTR

## Streamlining differential exon and 3' UTR usage analysis

The _diffUTR_ R package streamlines differential exon usage (DEU) analyses, and leverages existing DEU tools and alternative poly-adenylation site databases to enable differential 3' UTR usage analysis ([Gerber et al., 2021](https://doi.org/10.1186/s12859-021-04114-7)) . 

<img src="inst/docs/figure1.svg" alt="diffUTR scheme" style="margin-top: 15px; margin-bottom: 15px;" />

Popular bin-based DEU methods are provided by the [limma](https://bioconductor.org/packages/release/bioc/html/limma.html), [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) and in particular [DEXSeq](https://bioconductor.org/packages/release/bioc/html/DEXSeq.html) packages. However, their usage is not straightforward for non-experienced users, and their results often difficult to interpret. We therefore developed a simple workflow (Figure 1A), usable with any of the three methods but standardizing inputs and outputs. In particular, bin annotation and quantification, as well as different usage results, are all stored in a [RangedSummarizedExperiment](https://www.bioconductor.org/packages/devel/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html), which facilitates data storage and exploration, and enables advanced plotting functions irrespective of the underlying method. _diffUTR_ also provides an improved version of the `limma::diffSplice` method.

In addition, _diffUTR_ allows the extension of the DEU framework to UTR usage. A chief difficulty here is that most UTR variants are not catalogued in standard transcript annotations, limiting the utility of standard transcript-level quantification based on reference transcripts. However, based on databases of poly-adenylation (APA) sites such as [polyASite](https://polyasite.unibas.ch), _diffUTR_ can use alternative APA sites to further segment and extend UTR bins, as illustrated in Figure 1B.

In this way, _diffUTR_ outperforms alternative methods for detecting UTR changes from standard transcriptomics (see the [paper](https://doi.org/10.1186/s12859-021-04114-7) for more details) :

<img src="inst/docs/benchmark.png" alt="Differential UTR usage benchmark" width="600px" style="margin-top: 10px; margin-bottom: 15px;" />

**Note, however, that for all methods the FDR is considerably higher than the nominal one given by the method. For this reason, we urge users to use more stringent thresholds to avoid spurious results.**

Finally, _diffUTR_ provides a number of plotting utilities (see the vignette for more details), compatible with the results of any of the three underlying statistical methods.

## Installation

```{r}
BiocManager::install("ETHZ-INS/diffUTR")
```

If this fails because you don't have the latest R version, you can use:
```{r}
BiocManager::install("ETHZ-INS/diffUTR", ref="R36")
```

See the vignette for more details!
