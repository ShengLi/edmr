edmr
====

Comprehensive DMR analysis based on bimodal normal distribution model and cost function for regional methylation analysis.

Installation
---------
```R
install.packages( c("data.table", "mixtools" )
source("http://bioconductor.org/biocLite.R")
biocLite(c("GenomicRanges","IRanges"))
# install from github
library(devtools)
install_github("methylKit", username = "al2na",build_vignettes=FALSE)
install_github("edmr", username = "ShengLi",build_vignettes=FALSE)
```

Citing edmr
---------
Li S, Garrett-Bakelman FE, Akalin A, Zumbo P, Levine R, To BL, Lewis ID, Brown AL, D'Andrea RJ, Melnick A, Mason CE. An optimized algorithm for detecting and annotating regional differential methylation. BMC Bioinformatics. 2013;14 Suppl 5:S10.

