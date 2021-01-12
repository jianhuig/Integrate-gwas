# Integration of Functional Annotation Score into GWAS

Here we use childhood asthma as an example. Following are detailed code and procedure to integrate functional annotation score (CADD/Eigen) into GWAS. 

## 1. Download GWAS summary statistics

This GWAS result is perfomed by Nealelab based on UkBiobank data. See http://www.nealelab.is/uk-biobank.
```shell
wget https://www.dropbox.com/s/kp9bollwekaco0s/20002_1111.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O 20002_1111.gwas.imputed_v3.both_sexes.tsv.bgz
rename 20002_1111.gwas.imputed_v3.both_sexes.tsv.bgz 20002_1111.gwas.imputed_v3.both_sexes.tsv.gz
gunzip 20002_1111.gwas.imputed_v3.both_sexes.tsv.gz
```
## 2. QC control

To accomodate availabity of meta function score Eigen, we excluded InDel variants and SNPs on X-chromosome. In addition, we excluded SNPs with minor allele 
frequency less than 5%.

## 3. Obatin CADD/Eigen score

We use ANNOVAR to efficiently obatin Eigen scores on each SNP. ANNOVAR can be freely downloaded at http://download.openbioinformatics.org/annovar_download_form.php.
Unfortunately Eigen is not included in ANNOVAR by defualt, we need to download it manually from ANNOVAR database.
```shell
tar xvfz annovar.latest.tar.gz 
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar -eigen humandb/
perl annotate_variation.pl -filter asthma.variants.avinput humandb/ -dbtype eigen -build hg19
```
CADD score can be easily computed by running scripts from https://cadd.gs.washington.edu/download.

## 4. Integration Methods

Now your gwas dataset should contain at least the following 7 columns "chr","pos","ref","alt","tstat","pval","Eigen/CADD".

Function.R contains four integration methods: inverse-variance meta-analysis, Fisher's method, straitified FDR and weighted p-value approach. For meta-analysis and Fisher's method to run, we must obtain p-values for functional annotations scores. However, the functional annoations are 

i) not necessarily normal

ii) negative values indicate "less deleterious" and therefore should have large p-values.

