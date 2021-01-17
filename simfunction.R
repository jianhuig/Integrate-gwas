library(dplyr)
# define p value function
pvalue <- function(x) 2*pnorm(-abs(x))

# Estimate pi0 using lambda
compute_pi0 <- function(m, pvalue,lambda=0.5){
  min(sum(pvalue>lambda)/(m*(1-lambda)),1)
}

# GWAS Simulation
GWAS_simulation <- function(mu1, mu2, sd = 1, n = 100, random=FALSE, n.random=100) {
  # GWAS1
  z1 <-
    cbind(matrix(
      rnorm(n * rep, mean = mu1, sd = 1),
      byrow = T,
      nrow = rep
    ),
    matrix(
      rnorm((10000-n) * rep, mean = 0, sd = 1),
      byrow = T,
      nrow = rep
    ))
  # Randomly shuffle
  #z1 <- t(apply(z1, 1, sample))
  if (random){
    z1[,(101-n.random):10^4] <- t(apply(z1[,(101-n.random):10^4], 1, sample))
  }
  
  # GWAS2
  z2 <-
    cbind(matrix(
      rnorm(100 * rep, mean = mu2, sd = 1),
      byrow = T,
      nrow = rep
    ),
    matrix(
      rnorm(9900 * rep, mean = 0, sd = 1),
      byrow = T,
      nrow = rep
    ))
  # Pvalue
  p1 <- t(apply(z1, 1, pvalue))
  p2 <- t(apply(z2, 1, pvalue))
  return (list(z1, z2, p1, p2))
}

rej <- function(x,n){
  return(sum(x<=n)/ntrue)
}

rank.def <- function(x){
  x.ordered <- order(x)
  return(rank(x.ordered, ties.method = "first"))
}

Baseline <- function(x){
  q <-qvalue(x, pi0 = 0.99)$qvalue # compute q-value
  q.rank <- rank.def(q)
  return(c(rej(x[1:100],5e-6),rej(q.rank[1:100], 100),
         rej(q[1:100], 0.05),rej(q[1:100], 0.2),mean(q.rank[1:ntrue])))
}

Fisher <- function(x1,x2){
  chi <- -2*(log(x1)+log(x2))
  p <- pchisq(chi,df=4,lower.tail = F)
  q <- qvalue(p, pi0=0.99)$qvalue
  q.rank <- rank.def(q)
  return(c(rej(p[1:100],5e-6),rej(q.rank[1:100], 100),
           rej(q[1:100], 0.05), rej(q[1:100], 0.2),mean(q.rank[1:ntrue])))
}  


sFDR<- function(x, z2,i,q.level=0.05){
  priority <- which(x < sfdr.c) # high-priorty index
  strata1 <- z2[i, priority]
  strata2 <- z2[i, which(x >= sfdr.c)]
  pstrata1 <- pvalue(strata1) # p-value of strata1
  pstrata2 <- pvalue(strata2) # p-value of strata2
  m1 = length(pstrata1)
  m2 = length(pstrata2)
  # Estimate pi0 using lambda = 0.5
  pi0.strata1 <-
    compute_pi0(m = m1, pvalue = pstrata1)
  pi0.strata2 <-
    compute_pi0(m = m2, pvalue = pstrata2)
  # compute q-value
  q.strat1 <- qvalue(pstrata1, pi0 = pi0.strata1)$qvalue
  q.strat2 <- qvalue(pstrata2, pi0 = pi0.strata2)$qvalue
  if ((length(pstrata1[q.strat1 < q.level]) == 0) & (length(pstrata2[q.strat2 < q.level]) == 0)){
    p = pvalue(z2[i,])
  } else if((length(pstrata1[q.strat1 < q.level]) == 0)){
    alpha2 = max(pstrata2[q.strat2 < q.level])
    alpha1 = 0
    w2 = 10000/m2
    w1 = 0
    p = pvalue(z2[i,])
    p[priority] = p[priority]/w1
    p[which(x >= sfdr.c)] = p[which(x >= sfdr.c)]/w2
  } else if((length(pstrata2[q.strat2 < q.level]) == 0)){
    alpha1 = max(pstrata1[q.strat1 < q.level])
    alpha2 = 0
    w1 = 10000/m1
    w2 = 0
    p = pvalue(z2[i,])
    p[priority] = p[priority]/w1
    p[which(x >= sfdr.c)] = p[which(x >= sfdr.c)]/w2
  }
  else{
    alpha1 = max(pstrata1[q.strat1 < q.level])
    alpha2 = max(pstrata2[q.strat2 < q.level])
    alpha_bar = (m1*alpha1+m2*alpha2)/(m1+m2)
    w1 = alpha1/alpha_bar
    w2 = alpha2/alpha_bar
    p = pvalue(z2[i,])
    p[priority] = p[priority]/w1
    p[which(x >= sfdr.c)] = p[which(x >= sfdr.c)]/w2
  }
  q<- rep(1, nsnp)
  q[priority] <- q.strat1
  q[which(x >= sfdr.c)] <- q.strat2
  q.rank <- rank.def(q)
  return(c(rej(p[1:100],5e-6),rej(q.rank[1:100], 100),
           rej(q[1:100], 0.05), rej(q[1:100], 0.2),
           mean(q.rank[1:ntrue])))
}

wFDR <- function(x, z1,i){
  w = nsnp * (exp(z1[i,]) / sum(exp(z1[i,]))) # compute weight use z from GWAS1
  #w = nsnp * (abs(exp(z1[i,])) / sum(abs(exp(z1[i,]))))
  p <- x / w
  index <- rank(p)
  p_ordered <- sort(p)
  q_ordered <- p_ordered
  q_ordered[10000] <- 0.99 * q_ordered[10000]
  for ( i in 9999:1){
    q_ordered[i] <- min((9900/i)*q_ordered[i], q_ordered[i+1])
  }
  q <- q_ordered[index]
  q.rank <- rank.def(q)
  return(c(rej(p[1:100],5e-6),rej(q.rank[1:100], 100),
           rej(q[1:100], 0.05),rej(q[1:100], 0.2),
           mean(q.rank[1:ntrue])))
}

inv <- function(x1,x2){
  z <- (x1+x2)/sqrt(2)
  p <- pvalue(z)
  q <- qvalue(p)$qvalue
  q.rank <- rank.def(q)
  return(c(rej(p[1:100],5e-6),rej(q.rank[1:100], 100),
           rej(q[1:100], 0.05),rej(q[1:100], 0.2),
           mean(q.rank[1:ntrue])))
}
  
combine <- function(GWAS){
  z1 <- GWAS[[1]]
  z2 <- GWAS[[2]]
  p1 <- GWAS[[3]]
  p2 <- GWAS[[4]]
  baseline <- NULL
  fisher <- NULL
  strat <- NULL
  wei <- NULL
  meta <- NULL
  for (i in 1:rep) {
    # Baseline 
    baseline <- rbind(baseline, Baseline(p2[i,]))
    
    # Fisher
    fisher <- rbind(fisher, Fisher(p1[i,],p2[i,]))
    
    # SFDR
    strat <- rbind(strat, sFDR(p1[i,], z2,i))
    
    #WFDR
    wei <- rbind(wei, wFDR(p2[i,],z1,i))
    
    # Inverse Variance Meta
    meta <- rbind(meta, inv(z1[i,],z2[i,]))
  }
  result <- cbind(baseline, fisher, strat, wei, meta)
  return(apply(result, 2, mean))
}

illustration.plot <- function(mu1, mu2){
  par(mfrow= c(2,1))
  z2 <- c(rnorm(ntrue, mu2,1),rnorm(10^4-ntrue))
  plot(-log10(2*pnorm(-abs(z2))),ylab=expression(-log[10](P)),xlab="SNPs", main  = "Primary Dataset",ylim=c(0,10),col='grey')
  points(-log10(2*pnorm(-abs(z2[1:100]))),col='red')
  z1 <- c(rnorm(ntrue, mu1, 1),rnorm(10^4-ntrue, 0, 1))
  plot(-log10(2*pnorm(-abs(z1))), ylab=expression(-log[10](P)), xlab="SNPs", main= "Other Information",ylim=c(0,10),col='grey')
  points(-log10(2*pnorm(-abs(z1[0:ntrue]))))
}

power.plot <- function(result, mu2, n){
  #title <- paste0("Power at mu1 = ",mu1)
  sub <- paste0("Number of true SNPs = ", ntrue)
  # reference line
  plot(mu2,result[,1+n],type='o',col='red',xlab=expression(mu),ylab = "power",lwd=1, ylim=c(0,1), main =title, sub = sub, pch=16, bg='red')
  # Meta
  points(mu2,result[,21+n], col='black',pch=17)
  lines(mu2,result[,21+n], col='black',lty=3,lwd=1)
  # Fisher
  points(mu2,result[,6+n],pch=25, bg='black')
  lines(mu2,result[,6+n], col='black',lty=3,lwd=1)
  #SFDR 
  points(mu2,result[,11+n], col='blue',pch=22, bg='blue')
  lines(mu2,result[,11+n], col='blue',lty=3,lwd=1)
  points(mu2,result[,16+n], col='blue',pch=23, bg='blue')
  lines(mu2,result[,16+n], col='blue',lty=3,lwd=1)
 
  legend("topleft",legend=c("Baseline","Inverse-variance Meta","Fisher's Method","Weighted p-value","Stratified FDR"), col=c("red","black","black","blue","blue"), lty = c(1,3,3,3,3), ncol=1, pch=c(16,17,25,23,22),pt.bg = c("red","black","black","blue","blue"))
}

rank.plot <- function(result, mu2){
  sub <- paste0("Number of true SNPs = ", ntrue)
  # reference line
  plot(mu2,result[,5],type='o',col='red',xlab=expression(mu),ylab = "Ranking",lwd=1, ylim=rev(c(0,max(result[,25]))), main=title, sub = sub, pch=16, bg='red')
  # Meta
  points(mu2,result[,25], col='black',pch=17)
  lines(mu2,result[,25], col='black',lty=3,lwd=1)
  # Fisher
  points(mu2,result[,10],pch=25, bg='black')
  lines(mu2,result[,10], col='black',lty=3,lwd=1)
  #SFDR 
  points(mu2,result[,15], col='blue',pch=22, bg='blue')
  lines(mu2,result[,15], col='blue',lty=3,lwd=1)
  points(mu2,result[,20], col='blue',pch=23, bg='blue')
  lines(mu2,result[,20], col='blue',lty=3,lwd=1)
  
  legend("topleft",legend=c("Baseline","Inverse-variance Meta","Fisher's Method","Weighted p-value","Stratified FDR"), col=c("red","black","black","blue","blue"), lty = c(1,3,3,3,3), ncol=1, pch=c(16,17,25,23,22),pt.bg = c("red","black","black","blue","blue"))
}

gamma.plot <- function(result, mu2, n){
  sub <- paste0("Number of true SNPs = ", ntrue)
  # reference line
  plot(mu2,result[,3+n],type='o',col='red',xlab=expression(mu),ylab = "power",lwd=1, ylim=c(0,1), main =title, sub = sub, pch=16, bg='red')
  # Meta
  points(mu2,result[,23+n], col='black',pch=17)
  lines(mu2,result[,23+n], col='black',lty=3,lwd=1)
  # Fisher
  points(mu2,result[,8+n],pch=25, bg='black')
  lines(mu2,result[,8+n], col='black',lty=3,lwd=1)
  #SFDR 
  points(mu2,result[,13+n], col='blue',pch=22, bg='blue')
  lines(mu2,result[,13+n], col='blue',lty=3,lwd=1)
  points(mu2,result[,18+n], col='blue',pch=23, bg='blue')
  lines(mu2,result[,18+n], col='blue',lty=3,lwd=1)
  
  legend("topleft",legend=c("Baseline","Inverse-variance Meta","Fisher's Method","Weighted p-value","Stratified FDR"), col=c("red","black","black","blue","blue"), lty = c(1,3,3,3,3), ncol=1, pch=c(16,17,25,23,22),pt.bg = c("red","black","black","blue","blue"))
}
