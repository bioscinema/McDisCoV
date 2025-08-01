# McDisCoV

**McDisCoV**: *Outcome-linked Microbiome Discovery via Compositional lOgit-normal modeling under Varying support*

McDisCoV is an R package for differential abundance testing of microbiome compositional data using a robust logit-normal model. It leverages phylogenetic structure and performs feature selection under varying levels of support. McDisCoV is especially designed for datasets where features are sparse and zero-inflated, such as in 16S rRNA and shotgun metagenomic data.

---

## 📦 Installation

To install McDisCoV from GitHub:

```r
# You need devtools or remotes to install from GitHub
install.packages("devtools")  # or install.packages("remotes")
devtools::install_github("bioscinema/McDisCoV")
```

You can implement McDisCov following the code:
```r
# Load required packages
library(McDisCoV)
library(phyloseq)

# Run McDisCoV differential abundance analysis
res <- differential_abundance(
  physeq                = gp,
  group_var             = "Group",       # must match sample_data column
  imputation            = TRUE,
  zero_ratio_threshold  = 0.9,
  zero_impute           = 0.5,
  Wc                    = 0.7,
  max_iter              = 20,
  tol                   = 1e-6,
  min_nonzero           = 5
)

# View top differential taxa
head(res$wald_results)

# View grouped taxa (Group A, B, C)
table(res$group_results$group)
```
