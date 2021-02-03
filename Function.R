library(dplyr)
library(qvalue)

pvalue <- function(x){return(2*pnorm(-abs(x)))}

INT <- function(anno,tstat,k=3/8){
	sorted <-  sort(anno, decreasing = T,index.return=T)
	n = length(anno)
	p.sorted <- (seq(1:n)-k)/(n-2*k+1) # default k = 3/8, Bloom offset
	inv <- rep(0,n)
	inv[sorted$ix] <- abs(qnorm(p.sorted/2))
	sign = tstat * inv # manually fix the sign for meta-analysis
	inv[which(sign<0)]<- -inv[which(sign<0)]
	return(inv)
}

# Traditional Inverse-variance Meta Analysis
Meta <- function(z1, z2){
	meta.z <- (z1+z2)/sqrt(var(z1)+var(z2))
	return(pvalue(meta.z))
}

# Fisher's combination method
Fisher <- function(pval, anno){
	fisher.t <- -2*(log(pval)+log(pvalue(anno)))
	return(pchisq(gwas$fisher,df=4,lower.tail = F))
}

# sFDR
sFDR <- function(pval, anno, c = 0.05, gamma = 0.05, nbins=2){
	gwas <- data.frame(pval=pval,anno=anno)
    gwas$group <- with(gwas, cut(anno, breaks=quantile(anno, probs=seq(0,1, length.out=(nbins+1)), na.rm=TRUE), include.lowest=TRUE, labels=c(1:nbins)))
    gwas <- gwas %>% group_by(group) %>% mutate(q=qvalue(pval)$qval)
    alpha = rep(0, nbins)
    m <- sapply(levels(gwas$group),function(x) nrow(gwas %>% filter(group == x) %>% filter(q < 0.05)))
    if(sum(m)==0){w = rep(1,size) # unweighted case
    }else{
   	 index <- which(as.numeric(m)>0)
   	 alpha[index] <- sapply(levels(gwas$group)[index],function(x) max(gwas %>% filter(group == x) %>% filter(q < 0.05) %>% ungroup() %>% select(pval)))
   	 alpha_bar <- sum(m*alpha)/sum(m)
   	 w <- alpha/alpha_bar
    }
    gwas <- gwas %>% mutate(p.sfdr = pval/w[as.numeric(group)])
	return(gwas$p.sfdr)
}

# wFDR
wFDR <- function(pval, anno){
	w = length(pval)*(exp(abs(anno)))/ sum(exp(abs(anno)))
	return(pval/w)
}