# polyPK
:exclamation: This is a fork from read-only mirror of the CRAN R package repository.  polyPK — The Pharmacokinetics (PK) of Multi-Component Drugs Using a Metabolomics Approach  

- fix issue to make sure it can be used on R 4.2.X

## Install
install.packages(c("BiocManager", "remotes", "xlsx"))

BiocManager::install(c("impute", "pcaMethods", "ropls", "mixOmics"))

remotes::install_github("plyush1993/polyPK", force = T, upgrade = F, INSTALL_opts="--no-multiarch")
