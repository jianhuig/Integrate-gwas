library(qvalue)
library(ggplot2)
library(ggpubr)
sfdr.c = 0.05 # sFDR threshold
pvalue <- function(x) 2*pnorm(-abs(x))
GwasSimulation <- function(mu1, mu2, sd = 1, n = 100, random=FALSE, n.random=100) {
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
# Homogeneity Case
ntrue= 100
nsnp=10^4
rep = 10^3
mu1 = 3
mu2 = 3
G <- GwasSimulation(mu1, mu2)
homo <- combine(G)

# Hetero 1 m = 100,mu=1.5
ntrue= 100
nsnp=10^4
rep = 10^3
mu1 = 1.5
mu2 = 3
G <- GwasSimulation(mu1, mu2)
h1 <-combine(G)

# Hetero m=50, mu =3
ntrue= 100
nsnp=10^4
rep = 10^3
mu1 = 3
mu2 = 3
G <- GwasSimulation(mu1, mu2, n=50)
h2 <-combine(G)

# Hetero n=50 + 1.5
ntrue= 100
nsnp=10^4
rep = 10^3
mu1 = 1.5
mu2 = 3
G <- GwasSimulation(mu1, mu2,n=50)
h3 <-combine(G)

# Hetero mu=0, mu =0
ntrue= 100
nsnp=10^4
rep = 10^3
mu1 = 0
mu2 = 3
G <- GwasSimulation(mu1, mu2)
h4 <-combine(G)

# misleading m=50+50, mu =3 
ntrue= 100
nsnp=10^4
rep = 10^3
mu1 = 3
mu2 = 3
G <- GwasSimulation(mu1, mu2,random = T,n.random = 50)
h5 <-combine(G)


# m=50+50, mu=1.5
ntrue= 100
nsnp=10^4
rep = 10^3
mu1 = 1.5
mu2 = 3
G <- GwasSimulation(mu1, mu2,random = T,n.random = 50)
h6 <-combine(G)

# Misleading, mu =3
ntrue= 100
nsnp=10^4
rep = 10^3
mu1 = 3
mu2 = 3
G <- GwasSimulation(mu1, mu2,random = T,n.random = 100)
h7 <-combine(G)

comp.result <- rbind(homo,h1,h2,h3,h5,h6,h7,h4)

temp <- data.frame('RE'=c(1-comp.result[,25]/comp.result[1,5],
                          1-comp.result[,10]/comp.result[1,5],
                          1-comp.result[,15]/comp.result[1,5],
                          1-comp.result[,20]/comp.result[1,5]),
                   method = c(rep("meta",8),
                               rep("Fisher",8),
                               rep("sFDR",8),
                               rep("weighted p",8)))
temp$method <- factor(temp$method)
temp$method <- factor(temp$method,levels(temp$method)[c(2,1,4,3)])
temp %>% ggplot(aes(x=rep(1:8,4),y=RE,group=method)) +
  geom_point(aes(shape=method,color=method,fill=method))+ 
  geom_line(aes(color=method,linetype=method))+
  scale_shape_manual(values=c(24, 25, 23,22))+
  scale_color_manual(values = c("#B79F00","#619CFF","#00BFC4", "#C77CFF"))+
  scale_fill_manual(values = c("#B79F00","#619CFF","#00BFC4", "#C77CFF"))+
  ylab("Relative Efficiency")+
  xlab("Simulation Senario")+
  theme_classic()+
  theme(legend.position = "bottom")+
  geom_hline(yintercept = 0,color='red')+
  geom_bracket(xmin=1,xmax=1,label="Homogenity",y.position = -0.9)+
  ylim(-1,1)+scale_x_continuous(limits = c(0,9),breaks = 1:8)+
  geom_bracket(xmin = 2,xmax = 4,label = "Partial Informative",y.position = -0.9)+
  geom_bracket(xmin = 5,xmax = 7,label = "Partial Misleading",y.position = -0.9)+
  geom_bracket(xmin = 8,xmax = 8,label = "Uninformative",y.position = -0.9)+
  scale_linetype_manual(values = c("solid","dashed","solid","dashed"))
