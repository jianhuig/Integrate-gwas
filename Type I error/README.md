# Family-wise error rate

To emperically evaluate FWER of weighted p-value approach and sFDR, we simulate GWAS z-scores from N(0,1) for variants in 1000 Genome Project. For each simulation, we obatin CADD-score and integrate. The FWER is defined as (#simulations that have at least one findings/ #simulations). We also compare the FWER with recently purposed [FINDOR](https://github.com/gkichaev/FINDOR) method.

# Data Preparation

We obtain pre-cleaned 1000 Genome Project PLINK file for indepdent samples from [The Center for Applied Genomics at Sickkids](http://tcag.ca/tools/1000genomes.html)

```shell
wget http://tcag.ca/documents/tools/omni25_indep.tar.gz
tar -xzvf 
```


```r
library(bigsnpr)
library(dplyr)

# read PLINK file into R
snp_readBed("indep.bed")
bigSNP <- snp_attach("indep.rds")

n = nrow(bigSNP$fam) # Sample Size
sum_stat <- data.frame(SNP = bigSNP$map$marker.ID, A1 = bigSNP$map$allele1, A2 = bigSNP$map$allele2) %>% mutate( N = n)
load("cadd.RData")
id <- intersect(sum_stat$SNP, cadd$rsid) # only 545028 SNPs left
sum_stat <- sum_stat %>% filter(SNP %in% id)

# simulate 1000 sets
for ( i in 1:1000){
sum_stat <- sum_stat %>% mutate(Z = rnorm(nrow(sum_stat)))
gz1 <- gzfile(paste0("summary_statistics/1kg_sim",i,".tsv.gz"), "w")
write.table(sum_stat, gz1, row.names=F, quote=F, sep= " ")
close(gz1)
}
```
# FINDOR

You will need [ldsc](https://github.com/bulik/ldsc) and [findor](https://github.com/gkichaev/FINDOR).
```shell
for i in summary_statistics/*.gz
do
	sbatch --export=file=$i ldsc.sh
done
```
