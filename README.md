# Integration of Functional Annotation Score into GWAS

Here we use childhood asthma as an example. Following are detailed code and procedure to integrate functional annotation score (Eigen) into GWAS. 

## 1. Download GWAS summary statistics

This GWAS result is perfomed by Nealelab based on UkBiobank data. See http://www.nealelab.is/uk-biobank.
```console
$~/ wget https://www.dropbox.com/s/kp9bollwekaco0s/20002_1111.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O 20002_1111.gwas.imputed_v3.both_sexes.tsv.bgz
$~/ rename 20002_1111.gwas.imputed_v3.both_sexes.tsv.bgz 20002_1111.gwas.imputed_v3.both_sexes.tsv.gz
$~/ gunzip 20002_1111.gwas.imputed_v3.both_sexes.tsv.gz
```
## 2. QC control

To accomodate availabity of meta function score Eigen, we excluded InDel variants and SNPs on X-chromosome. In addition, we excluded SNPs with minor allele 
frequency less than 5%.
