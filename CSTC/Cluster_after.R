hcp_si = read.csv('~/Desktop/CSTC/after_hcp_si.csv', as.is = T)[,-1]
si_data = hcp_si[hcp_si$group == 0,]
ru_data = hcp_si[hcp_si$group == 1,]
cbic_data = hcp_si[hcp_si$group == 2,]

age_vec = seq(5,21,0.1)
new_age = data.frame(age_vec)
colnames(new_age) = 'age'

##########################################################################
#####################PREDICTION###########################################
##########################################################################


d_cbic = double()
r_square_cbic = c()
ssp_cbic=c()
for (i in colnames(cbic_data)[1:444]){
  gami = mgcv::gam(cbic_data[,i]~s(age,sp=-1), data = cbic_data)
  predict = predict(gami, newdata = new_age)
  d_cbic = rbind(d_cbic,predict)
  r_square_cbic = c(r_square_cbic, summary(gami)$r.sq)
  ssp_cbic = c(ssp_cbic,gami['sp'][[1]])
}
rownames(d_cbic) = colnames(cbic_data)[1:444]
names(r_square_cbic) = colnames(cbic_data)[1:444]
names(ssp_cbic) = colnames(cbic_data)[1:444]

ooo = rownames(d_cbic)[order(r_square_cbic, decreasing = T)]
pdf('~/Desktop/Jae_fucking_loser/after_data_cbic_auto.pdf')


for (i in ooo){
  df_pre = data.frame(age = age_vec, value = d_cbic[i,])
  df_data = data.frame(age = cbic_data[,'age'], value = cbic_data[,i])
  p2 = ggplot() + 
    geom_line(df_pre, mapping = aes(x=age, y=value), col='Red') +
    #geom_point(df, mapping = aes(x=age, y=v)) +
    geom_point(df_data, mapping = aes(x = age, y =value), alpha = 0.5)+
    ggtitle(paste('CBIC ',i,' r_square ', round(r_square_cbic[i],5), ' sp = ', round(ssp_cbic[i],3),sep=''))+ 
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  print(p2)
}



dev.off()



#####################################################
###################### CLUSTER TREE #################
#####################################################
library(dendextend)
library(colorspace)
library(ggplot2)



dist_mat_si <- cor(t(d_si))
dist_mat_si = 1 - dist_mat_si

dist_v_si = as.dist(dist_mat_si[1:169, 1:169])

dist_t_si = as.dist(dist_mat_si[383:444, 383:444])

hclust_v_si <- hclust(dist_v_si, method = 'average')
hclust_t_si <- hclust(dist_t_si, method = 'average')

dist_mat_ru <- cor(t(d_ru))
dist_mat_ru = 1 - dist_mat_ru
dist_v_ru = as.dist(dist_mat_ru[1:169, 1:169])

dist_t_ru = as.dist(dist_mat_ru[383:444, 383:444])



hclust_v_ru <- hclust(dist_v_ru, method = 'average')
hclust_t_ru <- hclust(dist_t_ru, method = 'average')

dist_mat_cbic <- cor(t(d_cbic))
dist_mat_cbic = 1 - dist_mat_cbic
dist_v_cbic = as.dist(dist_mat_cbic[1:169, 1:169])

dist_t_cbic = as.dist(dist_mat_cbic[383:444, 383:444])



hclust_v_cbic <- hclust(dist_v_cbic, method = 'average')
hclust_t_cbic <- hclust(dist_t_cbic, method = 'average')


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








setwd('~/clusters_analysis/')
dir.create('cluster_ana_after_t2')
setwd('cluster_ana_after_t2')
pdf('cluster_hcp.pdf')
for (c_number in c(2,3,4,5,6,7)){
  cut_t_si <- cutree(hclust_t_si, k=c_number)
  mean_t_si = double()
  for(i in c(1:c_number)){
    if (sum(cut_t_si==i) ==1) {
      after_scale = t(scale(d_si[383:444,][cut_t_si==i, ]))
    } else{
      after_scale = t(scale(t(d_si[383:444,][cut_t_si==i, ])))
    }
    mean = apply(after_scale, 2, mean, na.rm = T)
    mean_t_si = rbind(mean_t_si, mean)
  }
  clu_t_si = unique(cut_t_si[labels(hclust_t_si)])
  
  cut_t_ru <- cutree(hclust_t_ru, k=c_number)  
  mean_t_ru = double()
  for(i in c(1:c_number)){
    if (sum(cut_t_ru==i)==1) {
      after_scale = t(scale(d_ru[383:444,][cut_t_ru==i, ]))
    } else{
      after_scale = t(scale(t(d_ru[383:444,][cut_t_ru==i, ])))
    }
    mean = apply(after_scale, 2, mean, na.rm = T)
    mean_t_ru = rbind(mean_t_ru, mean)
  }
  
  
  ind_t_ru = as.integer(find_match(mean_t_ru, mean_t_si))
  clu_t_ru = unique(cut_t_ru[labels(hclust_t_ru)])
  
  
  cut_t_cbic <- cutree(hclust_t_cbic, k=c_number)  
  mean_t_cbic = double()
  for(i in c(1:c_number)){
    if (sum(cut_t_cbic==i)==1) {
      after_scale = t(scale(d_cbic[383:444,][cut_t_cbic==i, ]))
    } else{
      after_scale = t(scale(t(d_cbic[383:444,][cut_t_cbic==i, ])))
    }
    mean = apply(after_scale, 2, mean, na.rm = T)
    mean_t_cbic = rbind(mean_t_cbic, mean)
  }
  
  ind_t_cbic = as.integer(find_match(mean_t_cbic, mean_t_si))
  clu_t_cbic = unique(cut_t_cbic[labels(hclust_t_cbic)])
  
  
  
  col_all = rainbow_hcl(10,c=90)
  
  col_si = col_all[clu_t_si]
  col_ru = col_all[ind_t_ru[clu_t_ru]]
  col_cbic = col_all[ind_t_cbic[clu_t_cbic]]
  
  
  dir.create(paste(c_number,'_cluster',sep=''))
  par(mfrow=c(3,1))
  avg_dend_obj <- as.dendrogram(hclust_t_si)
  avg_col_dend <- color_branches(avg_dend_obj, k=c_number, col = col_si)
  plot(avg_col_dend, main = 'dendrogram for SI')
  avg_dend_obj <- as.dendrogram(hclust_t_ru)
  avg_col_dend <- color_branches(avg_dend_obj, k=c_number, col = col_ru)
  plot(avg_col_dend, main = 'dendrogram for RU')
  avg_dend_obj <- as.dendrogram(hclust_t_cbic)
  avg_col_dend <- color_branches(avg_dend_obj, k=c_number, col = col_cbic)
  plot(avg_col_dend, main = 'dendrogram for CBIC')
  
  
  
  
  dir.create(paste(c_number,'_cluster/','SI',sep = ''))
  dir.create(paste(c_number,'_cluster/','RU',sep = ''))
  dir.create(paste(c_number,'_cluster/','CBIC',sep = ''))
  plot_list=list()
  for(i in unique(cut_t_si)){
    dir.create(paste(c_number,'_cluster/','SI/cluster_',i,sep = ''))
    if (sum(cut_t_si==i)==1) {
      after_scale = t(scale(d_si[383:444,][cut_t_si==i, ]))
    } else{
      after_scale = t(scale(t(d_si[383:444,][cut_t_si==i, ])))
    }
    mean = apply(after_scale, 2, mean, na.rm = T)
    sd = apply(after_scale, 2, sd, na.rm = T)
    df = data.frame(age = age_vec, mean = mean, sd = sd)
    p = ggplot(df, aes(x=age, y=mean)) + 
      geom_line(col=col_all[i]) +
      geom_point(col=col_all[i])+
      geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.1,alpha = 0.1, col=col_all[i])+
      ylab('value')+labs(title = paste('cluster',i))
    plot_list[[3*i-2]]=p
    write(rownames(d_si[383:444,])[cut_t_si==i],paste(c_number,'_cluster/','SI/cluster_',i,'/variable.txt',sep = ''))
  }
  
  
  
  
  cut_t_ru_v = rep(0, length(cut_t_ru))
  for(i in 1:c_number){
    cut_t_ru_v[cut_t_ru==i] = ind_t_ru[i]
  }
  #col = rainbow_hcl(length(unique(cut_avg)),c=90)
  for(i in unique(cut_t_ru_v)){
    dir.create(paste(c_number,'_cluster/','RU/cluster_',i,sep = ''))
    if (sum(cut_t_ru_v==i)==1) {
      after_scale = t(scale(d_ru[383:444,][cut_t_ru_v==i, ]))
    } else{
      after_scale = t(scale(t(d_ru[383:444,][cut_t_ru_v==i, ])))
    }
    mean = apply(after_scale, 2, mean, na.rm = T)
    sd = apply(after_scale, 2, sd, na.rm = T)
    df = data.frame(age = age_vec, mean = mean, sd = sd)
    p = ggplot(df, aes(x=age, y=mean)) + 
      geom_line(col=col_all[i]) +
      geom_point(col=col_all[i])+
      geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.1,alpha = 0.1, col=col_all[i])+
      ylab('value')+labs(title = paste('cluster',i))
    plot_list[[3*i-1]]=p
    write(rownames(d_ru[383:444,])[cut_t_ru_v==i],paste(c_number,'_cluster/','RU/cluster_',i,'/variable.txt',sep = ''))
  }
  
  
  
  cut_t_cbic_v = rep(0, length(cut_t_cbic))
  for(i in 1:c_number){
    cut_t_cbic_v[cut_t_cbic==i] = ind_t_cbic[i]
  }
  #col = rainbow_hcl(length(unique(cut_avg)),c=90)
  for(i in unique(cut_t_cbic_v)){
    dir.create(paste(c_number,'_cluster/','CBIC/cluster_',i,sep = ''))
    if (sum(cut_t_cbic_v==i)==1) {
      after_scale = t(scale(d_cbic[383:444,][cut_t_cbic_v==i, ]))
    } else{
      after_scale = t(scale(t(d_cbic[383:444,][cut_t_cbic_v==i, ])))
    }
    mean = apply(after_scale, 2, mean, na.rm = T)
    sd = apply(after_scale, 2, sd, na.rm = T)
    df = data.frame(age = age_vec, mean = mean, sd = sd)
    p = ggplot(df, aes(x=age, y=mean)) + 
      geom_line(col=col_all[i]) +
      geom_point(col=col_all[i])+
      geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.1,alpha = 0.1, col=col_all[i])+
      ylab('value')+labs(title = paste('cluster',i))
    plot_list[[3*i]]=p
    write(rownames(d_cbic[383:444,])[cut_t_cbic_v==i],paste(c_number,'_cluster/','CBIC/cluster_',i,'/variable.txt',sep = ''))
  }
  
  
  
  multiplot(plotlist = plot_list, cols = c_number)  
  
  
}

dev.off()




#######################################################
############# Volume ##################################
#######################################################


setwd('~/clusters_analysis/')
dir.create('cluster_ana_after_v4')
setwd('cluster_ana_after_v4')
pdf('cluster_hcp.pdf')
for (c_number in c(2,3,4,5,6,7)){
  cut_t_si <- cutree(hclust_t_si, k=c_number)
  mean_t_si = double()
  for(i in c(1:c_number)){
    if (sum(cut_t_si==i) ==1) {
      after_scale = t(scale(d_si[383:444,][cut_t_si==i, ]))
    } else{
      after_scale = t(scale(t(d_si[383:444,][cut_t_si==i, ])))
    }
    mean = apply(after_scale, 2, mean, na.rm = T)
    mean_t_si = rbind(mean_t_si, mean)
  }
  clu_t_si = unique(cut_t_si[labels(hclust_t_si)])
  
  
  cut_v_si <- cutree(hclust_v_si, k=c_number)  
  mean_v_si = double()
  for(i in c(1:c_number)){
    if (sum(cut_v_si==i)==1) {
      after_scale = t(scale(d_si[1:169,][cut_v_si==i, ]))
    } else{
      after_scale = t(scale(t(d_si[1:169,][cut_v_si==i, ])))
    }
    mean = apply(after_scale, 2, mean, na.rm = T)
    mean_v_si = rbind(mean_v_si, mean)
  }
  
  
  ind_v_si = as.integer(find_match(mean_v_si, mean_t_si))
  clu_v_si = unique(cut_v_si[labels(hclust_v_si)])
  
  cut_v_ru <- cutree(hclust_v_ru, k=c_number)  
  mean_v_ru = double()
  for(i in c(1:c_number)){
    if (sum(cut_v_ru==i)==1) {
      after_scale = t(scale(d_ru[1:169,][cut_v_ru==i, ]))
    } else{
      after_scale = t(scale(t(d_ru[1:169,][cut_v_ru==i, ])))
    }
    mean = apply(after_scale, 2, mean, na.rm = T)
    mean_v_ru = rbind(mean_v_ru, mean)
  }
  
  
  ind_v_ru = as.integer(find_match(mean_v_ru, mean_t_si))
  clu_v_ru = unique(cut_v_ru[labels(hclust_v_ru)])
  
  
  cut_v_cbic <- cutree(hclust_v_cbic, k=c_number)  
  mean_v_cbic = double()
  for(i in c(1:c_number)){
    if (sum(cut_v_cbic==i)==1) {
      after_scale = t(scale(d_cbic[1:169,][cut_v_cbic==i, ]))
    } else{
      after_scale = t(scale(t(d_cbic[1:169,][cut_v_cbic==i, ])))
    }
    mean = apply(after_scale, 2, mean, na.rm = T)
    mean_v_cbic = rbind(mean_v_cbic, mean)
  }
  
  ind_v_cbic = as.integer(find_match(mean_v_cbic, mean_t_si))
  clu_v_cbic = unique(cut_v_cbic[labels(hclust_v_cbic)])
  
  
  
  col_all = rainbow_hcl(10,c=90)
  
  col_si = col_all[clu_t_si]
  col_si = col_all[ind_v_si[clu_v_si]]
  col_ru = col_all[ind_v_ru[clu_v_ru]]
  col_cbic = col_all[ind_v_cbic[clu_v_cbic]]
  
  
  dir.create(paste(c_number,'_cluster',sep=''))
  par(mfrow=c(3,1))
  avg_dend_obj <- as.dendrogram(hclust_v_si)
  avg_col_dend <- color_branches(avg_dend_obj, k=c_number, col = col_si)
  plot(avg_col_dend, main = 'dendrogram for SI')
  avg_dend_obj <- as.dendrogram(hclust_v_ru)
  avg_col_dend <- color_branches(avg_dend_obj, k=c_number, col = col_ru)
  plot(avg_col_dend, main = 'dendrogram for RU')
  avg_dend_obj <- as.dendrogram(hclust_v_cbic)
  avg_col_dend <- color_branches(avg_dend_obj, k=c_number, col = col_cbic)
  plot(avg_col_dend, main = 'dendrogram for CBIC')
  
  
  
  
  dir.create(paste(c_number,'_cluster/','SI',sep = ''))
  dir.create(paste(c_number,'_cluster/','RU',sep = ''))
  dir.create(paste(c_number,'_cluster/','CBIC',sep = ''))
  plot_list=list()
  for(i in unique(cut_v_si)){
    dir.create(paste(c_number,'_cluster/','SI/cluster_',i,sep = ''))
    if (sum(cut_v_si==i)==1) {
      after_scale = t(scale(d_si[1:169,][cut_v_si==i, ]))
    } else{
      after_scale = t(scale(t(d_si[1:169,][cut_v_si==i, ])))
    }
    mean = apply(after_scale, 2, mean, na.rm = T)
    sd = apply(after_scale, 2, sd, na.rm = T)
    df = data.frame(age = age_vec, mean = mean, sd = sd)
    p = ggplot(df, aes(x=age, y=mean)) + 
      geom_line(col=col_all[i]) +
      geom_point(col=col_all[i])+
      geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.1,alpha = 0.1, col=col_all[i])+
      ylab('value')+labs(title = paste('cluster',i))
    plot_list[[3*i-2]]=p
    write(rownames(d_si[1:169,])[cut_v_si==i],paste(c_number,'_cluster/','SI/cluster_',i,'/variable.txt',sep = ''))
  }
  
  
  
  
  cut_v_ru_v = rep(0, length(cut_v_ru))
  for(i in 1:c_number){
    cut_v_ru_v[cut_v_ru==i] = ind_v_ru[i]
  }
  #col = rainbow_hcl(length(unique(cut_avg)),c=90)
  for(i in unique(cut_v_ru_v)){
    dir.create(paste(c_number,'_cluster/','RU/cluster_',i,sep = ''))
    if (sum(cut_v_ru_v==i)==1) {
      after_scale = t(scale(d_ru[1:169,][cut_v_ru_v==i, ]))
    } else{
      after_scale = t(scale(t(d_ru[1:169,][cut_v_ru_v==i, ])))
    }
    mean = apply(after_scale, 2, mean, na.rm = T)
    sd = apply(after_scale, 2, sd, na.rm = T)
    df = data.frame(age = age_vec, mean = mean, sd = sd)
    p = ggplot(df, aes(x=age, y=mean)) + 
      geom_line(col=col_all[i]) +
      geom_point(col=col_all[i])+
      geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.1,alpha = 0.1, col=col_all[i])+
      ylab('value')+labs(title = paste('cluster',i))
    plot_list[[3*i-1]]=p
    write(rownames(d_ru[1:169,])[cut_v_ru_v==i],paste(c_number,'_cluster/','RU/cluster_',i,'/variable.txt',sep = ''))
  }
  
  
  
  cut_v_cbic_v = rep(0, length(cut_v_cbic))
  for(i in 1:c_number){
    cut_v_cbic_v[cut_v_cbic==i] = ind_v_cbic[i]
  }
  #col = rainbow_hcl(length(unique(cut_avg)),c=90)
  for(i in unique(cut_v_cbic_v)){
    dir.create(paste(c_number,'_cluster/','CBIC/cluster_',i,sep = ''))
    if (sum(cut_v_cbic_v==i)==1) {
      after_scale = t(scale(d_cbic[1:169,][cut_v_cbic_v==i, ]))
    } else{
      after_scale = t(scale(t(d_cbic[1:169,][cut_v_cbic_v==i, ])))
    }
    mean = apply(after_scale, 2, mean, na.rm = T)
    sd = apply(after_scale, 2, sd, na.rm = T)
    df = data.frame(age = age_vec, mean = mean, sd = sd)
    p = ggplot(df, aes(x=age, y=mean)) + 
      geom_line(col=col_all[i]) +
      geom_point(col=col_all[i])+
      geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.1,alpha = 0.1, col=col_all[i])+
      ylab('value')+labs(title = paste('cluster',i))
    plot_list[[3*i]]=p
    write(rownames(d_cbic[1:169,])[cut_v_cbic_v==i],paste(c_number,'_cluster/','CBIC/cluster_',i,'/variable.txt',sep = ''))
  }
  
  
  
  multiplot(plotlist = plot_list, cols = c_number)  
  
  
}

dev.off()




















