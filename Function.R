
IVT <- function(anno,tstat,k=3/8){
	sorted <-  sort(anno, decreasing = T,index.return=T)
	n = length(anno)
	p.sorted <- (seq(1:n)-k)/(n-2*k+1) # default k = 3/8, Bloom offset
	inv <- rep(0,n)
	inv[sorted$ix] <- abs(qnorm(p.sorted/2))
	sign = tstat * inv # manually fix the sign for meta-analysis
	inv[which(sign<0)]<- -inv[which(sign<0)]
	return(2*pnorm(-abs(inv)))
}

