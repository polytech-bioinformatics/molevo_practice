## Install devtools package, if necessary
install.packages("devtools")
## Install iMKT package from GitHub repository
devtools::install_github("sergihervas/iMKT")

## Load iMKT library
library(iMKT)

daf <- read.delim("jingwei_dteissieri_dyakuba_aligned.daf")
div <- read.delim("jingwei_dteissieri_dyakuba_aligned.div")
standardMKT(daf, div)

