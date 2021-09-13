# Family-wise error rate

To emperically evaluate FWER of weighted p-value approach and sFDR, we simulate GWAS z-scores from N(0,1) for variants in 1000 Genome Project. For each simulation, we obatin CADD-score and integrate. The FWER is defined as (#simulations that have at least one findings/ #simulations). We also compare the FWER with recently purposed [FINDOR](https://github.com/gkichaev/FINDOR) method.

# Data Preparation

We obtain the genotypes from UKBB data.

```shell
# Download .bim file
wget  -nd  biobank.ctsu.ox.ac.uk/ukb/ukb/auxdata/ukb_snp_bim.tar
tar -xvf ukb_snp_bim.tar

# Download .bed file per chromosome

```


```r
library(bigsnpr)
library(dplyr)
library(doParallel)

# read PLINK file into R
snp_readBed("indep.bed")
bigSNP <- snp_attach("indep.rds")

n = nrow(bigSNP$fam) # Sample Size
sum_stat <- data.frame(SNP = bigSNP$map$marker.ID, A1 = bigSNP$map$allele1, A2 = bigSNP$map$allele2) %>% mutate( N = n)
load("cadd.RData")
id <- intersect(sum_stat$SNP, cadd$rsid) # only 545028 SNPs left
sum_stat <- sum_stat %>% filter(SNP %in% id)

# simulate 1000 sets
cl <- parallel::makeCluster(detectCores()-1)
doParallel::registerDoParallel(cl)
foreach(i = 1:1000,.packages=c("dplyr")) %dopar% {
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
# Type I error evaluation
```r
library(data.table)
library(dplyr)
cl <- parallel::makeCluster(detectCores()-1)
doParallel::registerDoParallel(cl)
load("~/cadd.RData")
files <- list.files(path="~/summary_statistics/" ,pattern = "[.]reweighted$",full.names=TRUE)
result <- foreach(.packages = c("dplyr","qvalue"), i = files,.combine=rbind ,.export='fread',.errorhandling='pass') %dopar%{
    gwas <- fread(i) %>% rename(rsid = SNP)
    gwas <- gwas %>% left_join(cadd, by = 'rsid') %>% select(P, P_weighted, PHRED)
    gwas <- gwas %>% mutate(v = pnorm(PHRED-2)) %>% mutate(w = v/mean(v)) %>% mutate(pwfdr = P / w) %>% select(-v,-w)
    size = 2
    gwas$group <- with(gwas, cut(PHRED, breaks=quantile(PHRED, probs=c(0,0.05,1)), na.rm=TRUE, include.lowest=TRUE, labels=c(1:2)))
    gwas <- gwas %>% group_by(group) %>% mutate(q=qvalue(P)$qval)
    alpha = rep(0, size)
    m <- sapply(levels(gwas$group),function(x) nrow(gwas %>% filter(group == x) %>% filter(q < 0.05)))
    if(sum(m)==0){wei = rep(1,size) # unweighted case
    }else{
    	 index <- which(as.numeric(m)>0)
    	 alpha[index] <- sapply(levels(gwas$group)[index],function(x) max(gwas %>% filter(group == x) %>% filter(q < 0.05) %>% ungroup() %>% select(P)))
    	 alpha_bar <- sum(m*alpha)/sum(m)
    	 wei <- alpha/alpha_bar
    }
    gwas <- gwas %>% mutate(p.sfdr = P/wei[as.numeric(group)]) %>% ungroup()
	p.bon <- 0.05/nrow(gwas)
	c(sum(gwas$P_weighted < p.bon), sum(gwas$pwfdr < p.bon), sum(gwas$p.sfdr < p.bon))
	rm(gwas)
}
```
