# BLENDM

**BLEND-M: Cellular Deconvolution with Personalized DNA Methylation References**

## Overview

`BLENDM` is an R package for reference-based cellular deconvolution of bulk DNA methylation (DNAm) data. Unlike existing methods that assume a shared cell-type-specific (CTS) reference across all samples, BLEND-M learns a **personalized** CTS reference profile for each bulk sample by blending available purified reference data. It also explicitly models **heteroscedasticity** across marker CpGs via a two-step inverse-variance weighted non-negative least squares estimator, improving robustness to noisy markers.

## Key Features

- **Personalized references**: learns sample-specific CTS DNAm profiles rather than assuming a shared population reference
- **Heteroscedasticity modeling**: down-weights high-variance CpG sites for more robust fraction estimates
- **Two-step estimator**: computationally efficient and provably optimal (asymptotically equivalent to the oracle estimator with known variances)
- **Marker CpG selection**: built-in two-step procedure (ANOVA on M-values + effect size filter) for identifying informative CpGs from purified reference panels

## Ready-to-use reference list

Ready-to-use. Covers the marker CpGs selected using marker selection steps in BLEND-M manuscript. 

Two Choices. Essentially the same panel, just different resolutions.

- If just need general cell fractions: https://drive.google.com/file/d/1-Jwp7i7eAHRLVgB2v-XnRQ_pHXJ0KJ9O/view?usp=drive_link
- (Mostly useful for children blood DNAm deconvolution) If further need cord cell type fractions and adult blood cell type fractions: https://drive.google.com/file/d/1wrnUvJGrV7BIb0nnfP5KxwvtE7KARQSA/view?usp=drive_link

## Installation

```r
devtools::install_github("Penghuihuang2000/BLENDM")
```

Dependencies: `nnls`, `rlist`

## Usage

### Step 1: Select marker CpGs from a purified reference panel

```r
library(BLENDM)

# beta_ref:          G x M matrix of beta values from purified cell-type samples
# cell_type_labels:  character vector of cell type labels, length M

markers <- select_marker_cpg(
  beta_val    = beta_ref,
  cell_type   = cell_type_labels,
  pval_cutoff = 0.05,
  diff_cutoff = 10
)
```

### Step 2: Deconvolve bulk samples

```r
# reference.list: named list of matrices, one per cell type (G x M_t each)
# bulk_beta:      G x N matrix of bulk sample beta values

fractions <- BLENDM(
  mixture_sample = bulk_beta[markers, ],
  reference.list = lapply(reference.list, function(x) x[markers, ])
)

# fractions: N x T matrix (samples x cell types), rows sum to 1
```


## Citation

If you use BLENDM, please cite:

> Huang P, Peters DG, McKennan C. BLEND-M: cellular deconvolution with personalized DNA methylation references. 2025.

## License

MIT
