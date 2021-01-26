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
sFDR <- function(pval, anno, c = 0.05, gamma = 0.05, number_of_bins = 2){
	p.sfdr <- pval
	n <- length(pval)
	index1 <- which(anno < c)
	index2 <- which(anno >= c)
	m1 <- length(index1)
	m2 <- length(index2)
	q1 <- qvalue(pval[index1])$qvalue
	q2 <- qvalue(pval[index2])$qvalue
	if ((sum(q1 < gamma]) == 0) & ((sum(q2 < gamma]) == 0){
	  p.sfdr <- pval # unweighted case
	  } else if(sum(q1 < gamma]) == 0){
	  w1 = 0
	  w2 = n/m2
	  p.sfdr[index1] <- pval[index1]/w1
	  p.sfdr[index2] <- pval[index2]/w2 
	  } else if(sum(q2 < gamma]) == 0){
	  w2 = 0
	  w1 = n/m1
	  p.sfdr[index1] <- pval[index1]/w1
	  p.sfdr[index2] <- pval[index2]/w2 
	  } else{
	  alpha1 = max(p1[which(q1<gamma)])
	  alpha2 = max(p2[which(q2<gamma)])
	  alpha_bar = (m1*alpha1+m2*alpha2)/(m1+m2)
	  w1 = alpha1/alpha_bar
	  w2 = alpha2/alpha_bar
	  p.sfdr[index1] <- pval[index1]/w1
	  p.sfdr[index2] <- pval[index2]/w2
	  }
   return(p.sfdr)	
}

# wFDR
wFDR <- function(pval, anno){
	w = length(pval)*(exp(abs(anno)))/ sum(exp(abs(anno)))
	return(pval/w)
}
