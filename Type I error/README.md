# Family-wise error rate

To emperically evaluate FWER of weighted p-value approach and sFDR, we simulate GWAS z-scores from N(0,1) for variants in 1000 Genome Project. For each simulation, we obatin CADD-score and integrate. The FWER is defined as (#simulations that have at least one findings/ #simulations).

```console
wget http://tcag.ca/documents/tools/omni25_indep.tar.gz
tar -xzvf 
```
```r
library(bigsnpr)
rds <- snp_readBed("indep.bed")
```
