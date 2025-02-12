# supacow-share

This repository contains R scripts to generate most of the figures and other results presented in Kobel et al. 2025 ([preprint](https://doi.org/10.1101/2024.12.05.626740))

The scripts are written as R Markdown documents that can be compiled with [knitr](https://yihui.org/knitr/) to produce the figures and tables in the manuscript. The script `RCT_definition.Rmd` provides the computational definition for the two rumen community types (RCTs) compared in the manuscrip, based on [Dirichlet Multinomial Mixtures](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0030126) (DMM) clustering of 16S rRNA gene amplicon data. The rest of the scripts are named according to the figure/table they correspond to.

Input data for the scripts is available through the Norwegian National Infrastructure for Research Data (NIRD) at XXXXX.

Packages were managed with `renv`; see [their documentation](https://rstudio.github.io/renv/articles/renv.html) for how to use this.
