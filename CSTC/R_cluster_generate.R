hcp_ph = read.csv('Desktop/hcp_ph.csv', as.is = T)
hcp_ph = hcp_ph[, 3:ncol(hcp_ph)]

hcp_ph[,'iq'][is.na(hcp_ph[,'iq'])] = mean(hcp_ph[,'iq'],na.rm=T)
hcp_ph[,'int'][is.na(hcp_ph[,'int'])] = mean(hcp_ph[,'int'],na.rm=T)
hcp_ph[,'ext'][is.na(hcp_ph[,'ext'])] = mean(hcp_ph[,'ext'],na.rm=T)
hcp_ph[,'total'][is.na(hcp_ph[,'total'])] = mean(hcp_ph[,'total'],na.rm=T)


cbic_ph = read.csv('Desktop/cbic_123_ph.csv', as.is = T)
cbic_ph = cbic_ph[, 3:ncol(cbic_ph)]
library(stats)
    
    
cbic_ph[,'iq'][is.na(cbic_ph[,'iq'])] = mean(cbic_ph[,'iq'],na.rm=T)
cbic_ph[,'int'][is.na(cbic_ph[,'int'])] = mean(cbic_ph[,'int'],na.rm=T)
cbic_ph[,'ext'][is.na(cbic_ph[,'ext'])] = mean(cbic_ph[,'ext'],na.rm=T)
cbic_ph[,'total'][is.na(cbic_ph[,'total'])] = mean(cbic_ph[,'total'],na.rm=T)

install.packages('dendextend')

library(dendextend)
library(colorspace)
library(ggplot2)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


ru_hcp = hcp_ph[hcp_ph$group==1, ]
cbic_hcp = hcp_ph[hcp_ph$group==0, ]






age_vec = seq(5,21,0.1)
new_age = data.frame(age_vec)
colnames(new_age) = 'age'

#set.seed(1)
d = double()
r_square = c()
for (i in colnames(ru_hcp)[1:444]){
  gami = mgcv::gam(ru_hcp[,i]~s(age,sp=0.6), data = ru_hcp)
  predict = predict(gami, newdata = new_age)
  d = rbind(d,predict)
  r_square = c(r_square, summary(gami)$r.sq)
}
rownames(d) = colnames(ru_hcp)[1:444]

dist_mat <- cor(t(d))
dist_mat = 1 - dist_mat
dist_mat = as.dist(dist_mat)

hclust_avg <- hclust(dist_mat, method = 'average')


#set.seed(1)
d_v = double()
r_square_v = c()
for (i in colnames(cbic_hcp)[1:444]){
  gami = mgcv::gam(cbic_hcp[,i]~s(age,sp=0.6), data = cbic_hcp)
  predict = predict(gami, newdata = new_age)
  d_v = rbind(d_v,predict)
  r_square_v = c(r_square_v, summary(gami)$r.sq)
}
rownames(d_v) = colnames(cbic_hcp)[1:444]

dist_mat_v <- cor(t(d_v))
dist_mat_v = 1 - dist_mat_v
dist_mat_v = as.dist(dist_mat_v)

hclust_avg_v <- hclust(dist_mat_v, method = 'average')

find_match = function(matric1, matric2){
  rownames(matric1) = 1:nrow(matric1)
  rownames(matric2) = 1:nrow(matric2)
  cr = as.data.frame(cor(t(matric1), t(matric2),use='complete.obs'))
  cr = cbind(cr, rep(-1, nrow(cr)))
  or = c()
  for(i in 1:nrow(cr)){
    v = cr[i, !colnames(cr) %in% or]
    m = which.max(v)
    m = colnames(v)[m]
    or = c(or,m)
    #  fff = as.data.frame(fff)
  }
  return(or)
}


setwd('~/clusters_analysis/')
dir.create('cluster_ana03_cor')
setwd('cluster_ana03_cor')
pdf('cluster_hcp.pdf')
for (c_number in c(2,3,4,5,7)){
  cut_avg <- cutree(hclust_avg, k=c_number)
  mean_c = double()
  for(i in unique(cut_avg)){
    after_scale = t(scale(t(d[cut_avg==i, ])))
    mean = apply(after_scale, 2, mean, na.rm = T)
    mean_c = rbind(mean_c, mean)
  }
  
  cut_avg <- cutree(hclust_avg_v, k=c_number)  
  mean_cv = double()
  for(i in unique(cut_avg)){
    after_scale = t(scale(t(d_v[cut_avg==i, ])))
    mean = apply(after_scale, 2, mean, na.rm = T)
    mean_cv = rbind(mean_cv, mean)
  }
  ind = as.integer(find_match(mean_cv, mean_c))
  col = rainbow_hcl(c_number,c=90)
  col_v = col[ind]
  
  
  dir.create(paste(c_number,'_cluster',sep=''))
  par(mfrow=c(2,1))
  avg_dend_obj <- as.dendrogram(hclust_avg)
  avg_col_dend <- color_branches(avg_dend_obj, k=c_number, col = col)
  plot(avg_col_dend, main = 'dendrogram for hcp_ru')
  avg_dend_obj_v <- as.dendrogram(hclust_avg_v)
  avg_col_dend_v <- color_branches(avg_dend_obj_v, k=c_number, col = col_v)
  plot(avg_col_dend_v, main = 'dendrogram for hcp_cbic')
  
  
  
  cut_avg <- cutree(hclust_avg, k=c_number)
  dir.create(paste(c_number,'_cluster/','hcp_ru',sep = ''))
  dir.create(paste(c_number,'_cluster/','hcp_cbic',sep = ''))
  plot_list=list()
  for(i in unique(cut_avg)){
    dir.create(paste(c_number,'_cluster/','hcp_ru/cluster_',i,sep = ''))
    after_scale = t(scale(t(d[cut_avg==i, ])))
    mean = apply(after_scale, 2, mean, na.rm = T)
    sd = apply(after_scale, 2, sd, na.rm = T)
    df = data.frame(age = age_vec, mean = mean, sd = sd)
    p = ggplot(df, aes(x=age, y=mean)) + 
      geom_line(col=col[i]) +
      geom_point(col=col[i])+
      geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.1,alpha = 0.1, col=col[i])+
      ylab('value')+labs(title = paste('cluster',i))
    plot_list[[2*i-1]]=p
    write(rownames(d)[cut_avg==i],paste(c_number,'_cluster/','hcp_ru/cluster_',i,'/variable.txt',sep = ''))
    write(r_square[cut_avg==i],paste(c_number,'_cluster/','hcp_ru/cluster_',i,'/r_square.txt',sep = ''))
  }
  cut_avg <- cutree(hclust_avg_v, k=c_number)
  cut_avg_v = rep(0, length(cut_avg))
  for(i in 1:c_number){
    cut_avg_v[cut_avg==i] = ind[i]
  }
  #col = rainbow_hcl(length(unique(cut_avg)),c=90)
  for(i in unique(cut_avg)){
    dir.create(paste(c_number,'_cluster/','hcp_cbic/cluster_',i,sep = ''))
    after_scale = t(scale(t(d_v[cut_avg_v==i, ])))
    mean = apply(after_scale, 2, mean, na.rm = T)
    sd = apply(after_scale, 2, sd, na.rm = T)
    df = data.frame(age = age_vec, mean = mean, sd = sd)
    p = ggplot(df, aes(x=age, y=mean)) + 
      geom_line(col=col[i]) +
      geom_point(col=col[i])+
      geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.1,alpha = 0.1, col=col[i])+
      ylab('value')+labs(title = paste('cluster',i))
    plot_list[[2*i]]=p
    write(rownames(d)[cut_avg_v==i],paste(c_number,'_cluster/','hcp_cbic/cluster_',i,'/variable.txt',sep = ''))
    write(r_square_v[cut_avg==i],paste(c_number,'_cluster/','hcp_cbic/cluster_',i,'/r_square.txt',sep = ''))
  }
  multiplot(plotlist = plot_list, cols = c_number)  
  
  
}

dev.off()

dist_mat_v <- cor(t(d_v))
colnames(dist_mat_v) = NULL
rownames(dist_mat_v) = NULL
dd = cor(cbic_hcp[1:444], use='complete.obs')
colnames(dd) = NULL
rownames(dd) = NULL


dist_mat <- cor(t(d))
colnames(dist_mat) = NULL
rownames(dist_mat) = NULL
d = cor(ru_hcp[1:444], use='complete.obs')
colnames(d) = NULL
rownames(d) = NULL

library(reshape2)


melted_cormat <- melt(dist_mat_v)
melted_cormat_1 <- melt(dd)
melted_cormat_2 <- melt(dist_mat)
melted_cormat_3 <- melt(d)
pdf('corr.pdf')

library(ggplot2)
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+
  scale_fill_gradient(low="green", high="red", limits=c(-1, 1))+
  ggtitle("correlation of GAMs result cbic_hcp ")
ggplot(data = melted_cormat_1, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+
  scale_fill_gradient(low="green", high="red", limits=c(-1, 1))+
  ggtitle("correlation of raw data cbic_hcp")
ggplot(data = melted_cormat_2, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+
  scale_fill_gradient(low="green", high="red", limits=c(-1, 1))+
  ggtitle("correlation of GAMs result ru_hcp ")
ggplot(data = melted_cormat_3, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+
  scale_fill_gradient(low="green", high="red", limits=c(-1, 1))+
  ggtitle("correlation of raw data ru_hcp")

dev.off()

ru_hcp_res = sample_n(ru_hcp, nrow(cbic_hcp))
dist_mat <- cor(ru_hcp_res[1:444], use='complete.obs')
dist_mat = 1 - dist_mat
dist_mat = as.dist(dist_mat)

hclust_avg <- hclust(dist_mat, method = 'average')

dist_mat_v <- cor(cbic_hcp[1:444], use='complete.obs')
dist_mat_v = 1 - dist_mat_v
dist_mat_v = as.dist(dist_mat_v)

hclust_avg_v <- hclust(dist_mat_v, method = 'average')






#####################################################################
####################### R-SQUARE ####################################
#####################################################################
r_1 = c()
for (i in colnames(ru_hcp)[1:444]){
  gami = mgcv::gam(ru_hcp[,i]~s(age,sp=0.6), data = ru_hcp)
  r = summary(gami)$r.sq
  r_1 = c(r_1,r) 
}
#colnames(d) = colnames(ru_hcp)[1:444]

r_2 = c()
for (i in colnames(cbic_hcp)[1:444]){
  gami = mgcv::gam(cbic_hcp[,i]~s(age,sp=0.6), data = cbic_hcp)
  r = summary(gami)$r.sq
  r_2 = c(r_2,r)
}
#rownames(d_v) = colnames(cbic_hcp)[1:444]

r_s = c(r_1, r_2)
site = c(rep('ru', 444), rep('cbic', 444))
r_data = data.frame(r_s, site)




packages <- c("ggplot2", "dplyr", "lavaan", "plyr", "cowplot", "rmarkdown", 
              "readr", "caTools", "bitops")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

setwd('~/Desktop/')
png('r_square.png')
p6 = ggplot(r_data,aes(x=site,y=r_s, fill = site, colour = site))+
  geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =1, trim = FALSE)+
  geom_point(position = position_jitter(width = .15), size = .25)+
  geom_boxplot(aes(x = as.numeric(site)+0.25, y = r_s),outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK") +
  ylab('r_square')+xlab('site')+coord_flip()+theme_cowplot()+guides(fill = FALSE, colour = FALSE) +
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  ggtitle("r_square between hcp_ru & hcp_cbic")
print(p6)
dev.off()


dist_mat_v <- cor(t(d_v))
colnames(dist_mat_v) = NULL
rownames(dist_mat_v) = NULL
dd = cor(cbic_hcp[1:444], use='complete.obs')
colnames(dd) = NULL
rownames(dd) = NULL


dist_mat <- cor(t(d))
colnames(dist_mat) = NULL
rownames(dist_mat) = NULL
d = cor(ru_hcp[1:444], use='complete.obs')
colnames(d) = NULL
rownames(d) = NULL

library(ggplot2)

uptri_v = data.frame(c(as.numeric(as.dist(dist_mat_v)), as.numeric(as.dist(dd))), c(rep('RAW', length(as.numeric(as.dist(dist_mat_v)))), rep('GAM', length(as.numeric(as.dist(dist_mat_v))))))

uptri = data.frame(c(as.numeric(as.dist(dist_mat)), as.numeric(as.dist(d))), c(rep('RAW', length(as.numeric(as.dist(dist_mat)))), rep('GAM', length(as.numeric(as.dist(dist_mat))))))

colnames(uptri) = c('Value', 'type')
colnames(uptri_v) = c('Value', 'type')

pdf('hist.pdf')
ggplot(uptri, aes(x=Value, fill= type, color=type)) +
  geom_histogram(position="identity", bins = 50, alpha = 0.5)+labs(title=paste("ru_hcp",cor(as.numeric(as.dist(dist_mat)), as.numeric(as.dist(d)))))
ggplot(uptri_v, aes(x=Value, fill= type, color=type)) +
  geom_histogram(position="identity", bins = 50, alpha = 0.5)+labs(title=paste("cbic_hcp", cor(as.numeric(as.dist(dist_mat_v)), as.numeric(as.dist(dd)))))

dev.off()


mat1 = cor(t(rbind(d,d_v)))
colnames(mat1) = NULL
rownames(mat1) = NULL
#mat2 = cor(cbind(ru_hcp[1:444],cbic_hcp[1:444]))


melted_cormat1 <- melt(mat1)
#melted_cormat_2 <- melt(dist_mat)
png('correlation.png')

library(ggplot2)
ggplot(data = melted_cormat1, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+
  scale_fill_gradient(low="green", high="red", limits=c(-1, 1))+
  ggtitle("correlation of GAMs results for ru_hcp & cbic_hcp ")
dev.off()

cross=mat1[445:888,1:444]
png('~/Desktop/CSTC/off_dig_hist.png')
library(ggplot2)         
ggplot()+geom_histogram(mapping = (aes(x=diag(cross))), bins=40, col='red', alpha = 0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()
















