load("fullanalysis_cadd_random.RData")
df <- do.call(rbind.data.frame,result)
colnames(df) <- c("phenotype","m", paste0(c("p.","tp."),"meta"),
                  paste0(c("p.","tp."),"Fisher"),
                  paste0(c("p.","tp."),"sFDR"),
                  paste0(c("p.","tp."),"weighted p-value"))
df[,2:10] <- mutate_all(df[,2:10], function(x) as.numeric(as.character(x)))
df$phenotype <- gsub("\\..*","", basename(df$phenotype))
h2 <- data.table::fread("ukb31063_h2_topline.02Oct2019.tsv")
df <- df %>% left_join(h2, by='phenotype') %>% filter(!is.na(h2_sig))
m = df$m
df1 <- data.frame(h2_sig=rep(df$h2_sig,4),
                  h2=rep(df$h2_liability,4),
                  value= c(df$tp.meta/m, 
                           df$tp.Fisher/m,
                           df$`tp.weighted p-value`/m,
                           df$tp.sFDR/m),
                  method = c(rep("meta",1132),
                             rep("Fisher",1132),
                             rep("weighted p",1132),
                             rep("sFDR",1132)))
df1$method <- factor(df1$method, levels = c("meta","Fisher","weighted p","sFDR"))
p1 <- df1 %>% filter(value <=1) %>% ggplot(aes(x=method,y=value))+
  geom_boxplot(aes(fill=method)) +
  theme_bw()+
  theme(legend.position = 'none',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ylab("Recall")+
  xlab("")+
  geom_jitter(size=0.1)+
  scale_fill_manual(values=c("#B79F00","#619CFF","#00BFC4", "#C77CFF"))
df1 <- data.frame(h2_sig=rep(df$h2_sig,4),
                  h2=rep(df$h2_liability,4),
                  value= c(df$tp.meta/df$p.meta,
                           df$tp.Fisher/df$p.Fisher,
                           df$`tp.weighted p-value`/df$`p.weighted p-value`,
                           df$tp.sFDR/df$p.sFDR), 
                  method = c(rep("meta",1132),
                             rep("Fisher",1132),
                             rep("weighted p",1132),
                             rep("sFDR",1132)))
df1$method <- factor(df1$method, levels = c("meta","Fisher","weighted p","sFDR"))
p2 <- df1 %>% filter(value <= 1) %>% ggplot(aes(x = method, y = value))+ 
  geom_boxplot(aes(fill = method))+ 
  theme_bw()+ 
  theme(
  legend.position = 'none',
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)+ 
  ylab("Precision")+ 
  xlab("")+ 
  geom_jitter(size = 0.1)+ 
  scale_fill_manual(values = c("#B79F00", "#619CFF", "#00BFC4", "#C77CFF"))
ggpubr::ggarrange(p1,p2)
