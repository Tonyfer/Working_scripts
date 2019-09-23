conn_l = list()
si_l = list()
ru_l = list()
cbic_l = list()
for (i in c(1:7)){
  conn = read.csv(paste0('~/Desktop/CSTC/122_conn',i,'_after.csv'), as.is = T)[,-1]
  si_l[[i]] = conn[hcp_si$group == 0,]
  ru_l[[i]] = conn[hcp_si$group == 1,]
  cbic_l[[i]] = conn[hcp_si$group == 2,]
}
age_vec = seq(5,21,0.1)
new_age = data.frame(age_vec)
colnames(new_age) = 'age'

##########################################################################
#####################PREDICTION###########################################
##########################################################################

r_square_si_l = list()
r_square_ru_l = list()
r_square_cbic_l = list()
d_si_l = list()
d_ru_l = list()
d_cbic_l = list()

for (i in c(1:7)){
  d_si = double()
  r_square_si = c()
  ssp_si=c()
  for (c in colnames(si_l[[i]])[1:(dim(si_l[[i]])[2]-6)]){
    gami = mgcv::gam(si_l[[i]][,c]~s(age,sp=-1), data = si_l[[i]])
    predict = predict(gami, newdata = new_age)
    d_si = rbind(d_si,predict)
    r_square_si = c(r_square_si, summary(gami)$r.sq)
    ssp_si = c(ssp_si,gami['sp'][[1]])
  }
  rownames(d_si) = colnames(si_l[[i]])[1:(dim(si_l[[i]])[2]-6)]
  names(r_square_si) = colnames(si_l[[i]])[1:(dim(si_l[[i]])[2]-6)]
  names(ssp_si) = colnames(si_l[[i]])[1:(dim(si_l[[i]])[2]-6)]
  
  d_si_l[[i]] = d_si
  r_square_si_l[[i]] = r_square_si
  
  #ooo = rownames(d_si)[order(r_square_si, decreasing = T)]
  #pdf(paste0('~/Desktop/Jae_fucking_loser/after_444_conn',i,'_si_auto.pdf'))


  #for (o in ooo){
   # df_pre = data.frame(age = age_vec, value = d_si[o,])
  #  df_data = data.frame(age = si_l[[i]][,'age'], value = si_l[[i]][,o])
  #   p2 = ggplot() + 
  #     geom_line(df_pre, mapping = aes(x=age, y=value), col='Red') +
  #     #geom_point(df, mapping = aes(x=age, y=v)) +
  #     geom_point(df_data, mapping = aes(x = age, y =value), alpha = 0.5)+
  #     ggtitle(paste('SI ',o,' r_square ', round(r_square_si[o],5), ' sp = ', round(ssp_si[o],3),sep=''))+ 
  #     theme_bw() + 
  #     theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  #   print(p2)
  # }
  # dev.off()
  
  
  d_ru = double()
  r_square_ru = c()
  ssp_ru=c()
  for (c in colnames(ru_l[[i]])[1:(dim(ru_l[[i]])[2]-6)]){
    gami = mgcv::gam(ru_l[[i]][,c]~s(age,sp=-1), data = ru_l[[i]])
    predict = predict(gami, newdata = new_age)
    d_ru = rbind(d_ru,predict)
    r_square_ru = c(r_square_ru, summary(gami)$r.sq)
    ssp_ru = c(ssp_ru,gami['sp'][[1]])
  }
  rownames(d_ru) = colnames(ru_l[[i]])[1:(dim(ru_l[[i]])[2]-6)]
  names(r_square_ru) = colnames(ru_l[[i]])[1:(dim(ru_l[[i]])[2]-6)]
  names(ssp_ru) = colnames(ru_l[[i]])[1:(dim(ru_l[[i]])[2]-6)]
  
  d_ru_l[[i]] = d_ru
  r_square_ru_l[[i]] = r_square_ru
  
  #ooo = rownames(d_ru)[order(r_square_ru, decreasing = T)]
  #pdf(paste0('~/Desktop/Jae_fucking_loser/after_444_conn',i,'_ru_auto.pdf'))
  
  
  # for (o in ooo){
  #   df_pre = data.frame(age = age_vec, value = d_ru[o,])
  #   df_data = data.frame(age = ru_l[[i]][,'age'], value = ru_l[[i]][,o])
  #   p2 = ggplot() + 
  #     geom_line(df_pre, mapping = aes(x=age, y=value), col='Red') +
  #     #geom_point(df, mapping = aes(x=age, y=v)) +
  #     geom_point(df_data, mapping = aes(x = age, y =value), alpha = 0.5)+
  #     ggtitle(paste('RU ',o,' r_square ', round(r_square_ru[o],5), ' sp = ', round(ssp_ru[o],3),sep=''))+ 
  #     theme_bw() + 
  #     theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  #   print(p2)
  # }
  # dev.off()
  
  
  d_cbic = double()
  r_square_cbic = c()
  ssp_cbic=c()
  for (c in colnames(cbic_l[[i]])[1:(dim(cbic_l[[i]])[2]-6)]){
    gami = mgcv::gam(cbic_l[[i]][,c]~s(age,sp=-1), data = cbic_l[[i]])
    predict = predict(gami, newdata = new_age)
    d_cbic = rbind(d_cbic,predict)
    r_square_cbic = c(r_square_cbic, summary(gami)$r.sq)
    ssp_cbic = c(ssp_cbic,gami['sp'][[1]])
  }
  rownames(d_cbic) = colnames(cbic_l[[i]])[1:(dim(cbic_l[[i]])[2]-6)]
  names(r_square_cbic) = colnames(cbic_l[[i]])[1:(dim(cbic_l[[i]])[2]-6)]
  names(ssp_cbic) = colnames(cbic_l[[i]])[1:(dim(cbic_l[[i]])[2]-6)]
  
  d_cbic_l[[i]] = d_cbic
  r_square_cbic_l[[i]] = r_square_cbic
  
  # ooo = rownames(d_cbic)[order(r_square_cbic, decreasing = T)]
  # pdf(paste0('~/Desktop/Jae_fucking_loser/after_444_conn',i,'_cbic_auto.pdf'))
  # 
  # 
  # for (o in ooo){
  #   df_pre = data.frame(age = age_vec, value = d_cbic[o,])
  #   df_data = data.frame(age = cbic_l[[i]][,'age'], value = cbic_l[[i]][,o])
  #   p2 = ggplot() + 
  #     geom_line(df_pre, mapping = aes(x=age, y=value), col='Red') +
  #     #geom_point(df, mapping = aes(x=age, y=v)) +
  #     geom_point(df_data, mapping = aes(x = age, y =value), alpha = 0.5)+
  #     ggtitle(paste('CBIC ',o,' r_square ', round(r_square_cbic[o],5), ' sp = ', round(ssp_cbic[o],3),sep=''))+ 
  #     theme_bw() + 
  #     theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  #   print(p2)
  # }
  # dev.off()
  # 
}





library(dendextend)
library(colorspace)
library(ggplot2)

hclust_c_si_l = list()
hclust_c_ru_l = list()
hclust_c_cbic_l = list()


for (i in c(1:7)){
  dist_mat_si <- cor(t(d_si_l[[i]]))
  dist_mat_si = 1 - dist_mat_si

  dist_c_si = as.dist(dist_mat_si)

#dist_t_si = as.dist(dist_mat_si[383:444, 383:444])

  hclust_c_si_l[[i]] <- hclust(dist_c_si, method = 'average')
#hclust_t_si <- hclust(dist_t_si, method = 'average')

  dist_mat_ru <- cor(t(d_ru_l[[i]]))
  dist_mat_ru = 1 - dist_mat_ru
  dist_c_ru = as.dist(dist_mat_ru)

#dist_t_ru = as.dist(dist_mat_ru[383:444, 383:444])



  hclust_c_ru_l[[i]] <- hclust(dist_c_ru, method = 'average')
#hclust_t_ru <- hclust(dist_t_ru, method = 'average')

  dist_mat_cbic <- cor(t(d_cbic_l[[i]]))
  dist_mat_cbic = 1 - dist_mat_cbic
  dist_c_cbic = as.dist(dist_mat_cbic)

#dist_t_cbic = as.dist(dist_mat_cbic[383:444, 383:444])



  hclust_c_cbic_l[[i]] <- hclust(dist_c_cbic, method = 'average')
#hclust_t_cbic <- hclust(dist_t_cbic, method = 'average')
}

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





for (net in c(1:7)){
  dir.create(paste0('~/clusters_analysis/cluster_ana_122_conn', net, '_2'))
  setwd(paste0('~/clusters_analysis/cluster_ana_122_conn', net, '_2'))
  pdf('cluster_hcp.pdf')
  for (c_number in c(2:5)){
    cut_t_si <- cutree(hclust_c_si_l[[1]], k=c_number)
    mean_t_si = double()
    for(i in c(1:c_number)){
      if (sum(cut_t_si==i) ==1) {
        after_scale = t(scale(d_si_l[[1]][cut_t_si==i, ]))
      } else{
        after_scale = t(scale(t(d_si_l[[1]][cut_t_si==i, ])))
      }
      mean = apply(after_scale, 2, mean, na.rm = T)
      mean_t_si = rbind(mean_t_si, mean)
    }
    clu_t_si = unique(cut_t_si[labels(hclust_c_si_l[[1]])])
    
    
    cut_c_si <- cutree(hclust_c_si_l[[net]], k=c_number)  
    mean_c_si = double()
    for(i in c(1:c_number)){
      if (sum(cut_c_si==i)==1) {
        after_scale = t(scale(d_si_l[[net]][cut_c_si==i, ]))
      } else{
        after_scale = t(scale(t(d_si_l[[net]][cut_c_si==i, ])))
      }
      mean = apply(after_scale, 2, mean, na.rm = T)
      mean_c_si = rbind(mean_c_si, mean)
    }
    
    
    ind_c_si = as.integer(find_match(mean_c_si, mean_t_si))
    clu_c_si = unique(cut_c_si[labels(hclust_c_si_l[[net]])])
    
    cut_c_ru <- cutree(hclust_c_ru_l[[net]], k=c_number)  
    mean_c_ru = double()
    for(i in c(1:c_number)){
      if (sum(cut_c_ru==i)==1) {
        after_scale = t(scale(d_ru_l[[net]][cut_c_ru==i, ]))
      } else{
        after_scale = t(scale(t(d_ru_l[[net]][cut_c_ru==i, ])))
      }
      mean = apply(after_scale, 2, mean, na.rm = T)
      mean_c_ru = rbind(mean_c_ru, mean)
    }
    
    
    ind_c_ru = as.integer(find_match(mean_c_ru, mean_t_si))
    clu_c_ru = unique(cut_c_ru[labels(hclust_c_ru_l[[net]])])
    
    
    cut_c_cbic <- cutree(hclust_c_cbic_l[[net]], k=c_number)  
    mean_c_cbic = double()
    for(i in c(1:c_number)){
      if (sum(cut_c_cbic==i)==1) {
        after_scale = t(scale(d_cbic_l[[net]][cut_c_cbic==i, ]))
      } else{
        after_scale = t(scale(t(d_cbic_l[[net]][cut_c_cbic==i, ])))
      }
      mean = apply(after_scale, 2, mean, na.rm = T)
      mean_c_cbic = rbind(mean_c_cbic, mean)
    }
    
    ind_c_cbic = as.integer(find_match(mean_c_cbic, mean_t_si))
    clu_c_cbic = unique(cut_c_cbic[labels(hclust_c_cbic_l[[net]])])
    
    
    
    col_all = rainbow_hcl(10,c=90)
    
    col_si = col_all[clu_t_si]
    col_si = col_all[ind_c_si[clu_c_si]]
    col_ru = col_all[ind_c_ru[clu_c_ru]]
    col_cbic = col_all[ind_c_cbic[clu_c_cbic]]
    
    
    dir.create(paste(c_number,'_cluster',sep=''))
    par(mfrow=c(3,1))
    avg_dend_obj <- as.dendrogram(hclust_c_si_l[[net]])
    avg_col_dend <- color_branches(avg_dend_obj, k=c_number, col = col_si)
    plot(avg_col_dend, main = 'dendrogram for SI')
    avg_dend_obj <- as.dendrogram(hclust_c_ru_l[[net]])
    avg_col_dend <- color_branches(avg_dend_obj, k=c_number, col = col_ru)
    plot(avg_col_dend, main = 'dendrogram for RU')
    avg_dend_obj <- as.dendrogram(hclust_c_cbic_l[[net]])
    avg_col_dend <- color_branches(avg_dend_obj, k=c_number, col = col_cbic)
    plot(avg_col_dend, main = 'dendrogram for CBIC')
    
    
    
    
    dir.create(paste(c_number,'_cluster/','SI',sep = ''))
    dir.create(paste(c_number,'_cluster/','RU',sep = ''))
    dir.create(paste(c_number,'_cluster/','CBIC',sep = ''))
    plot_list=list()
    
    
    cut_c_si_v = rep(0, length(cut_c_si))
    for(i in 1:c_number){
      cut_c_si_v[cut_c_si==i] = ind_c_si[i]
    }
    for(i in unique(cut_c_si_v)){
      dir.create(paste(c_number,'_cluster/','SI/cluster_',i,sep = ''))
      if (sum(cut_c_si_v==i)==1) {
        after_scale = t(scale(d_si_l[[net]][cut_c_si_v==i, ]))
      } else{
        after_scale = t(scale(t(d_si_l[[net]][cut_c_si_v==i, ])))
      }
      mean = apply(after_scale, 2, mean, na.rm = T)
      sd = apply(after_scale, 2, sd, na.rm = T)
      df = data.frame(age = age_vec, mean = mean, sd = sd)
      p = ggplot(df, aes(x=age, y=mean)) + 
        geom_line(col=col_all[i]) +
        geom_point(col=col_all[i])+
        geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.1,alpha = 0.1, col=col_all[i])+
        ylab('value')+labs(title = paste('cluster',i,'\n N =',  sum(cut_c_si_v==i)))+
        theme(axis.text=element_text(size=6), axis.title=element_text(size=6))
      plot_list[[3*i-2]]=p
      write(rownames(d_si_l[[net]])[cut_c_si_v==i],paste(c_number,'_cluster/','SI/cluster_',i,'/variable.txt',sep = ''))
    }
    
    
    
    
    cut_c_ru_v = rep(0, length(cut_c_ru))
    for(i in 1:c_number){
      cut_c_ru_v[cut_c_ru==i] = ind_c_ru[i]
    }
    #col = rainbow_hcl(length(unique(cut_avg)),c=90)
    for(i in unique(cut_c_ru_v)){
      dir.create(paste(c_number,'_cluster/','RU/cluster_',i,sep = ''))
      if (sum(cut_c_ru_v==i)==1) {
        after_scale = t(scale(d_ru_l[[net]][cut_c_ru_v==i, ]))
      } else{
        after_scale = t(scale(t(d_ru_l[[net]][cut_c_ru_v==i, ])))
      }
      mean = apply(after_scale, 2, mean, na.rm = T)
      sd = apply(after_scale, 2, sd, na.rm = T)
      df = data.frame(age = age_vec, mean = mean, sd = sd)
      p = ggplot(df, aes(x=age, y=mean)) + 
        geom_line(col=col_all[i]) +
        geom_point(col=col_all[i])+
        geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.1,alpha = 0.1, col=col_all[i])+
        ylab('value')+labs(title = paste('cluster',i,'\n N =',  sum(cut_c_ru_v==i)))+
        theme(axis.text=element_text(size=6), axis.title=element_text(size=6))
      plot_list[[3*i-1]]=p
      write(rownames(d_ru_l[[net]])[cut_c_ru_v==i],paste(c_number,'_cluster/','RU/cluster_',i,'/variable.txt',sep = ''))
    }
    
    
    
    cut_c_cbic_v = rep(0, length(cut_c_cbic))
    for(i in 1:c_number){
      cut_c_cbic_v[cut_c_cbic==i] = ind_c_cbic[i]
    }
    #col = rainbow_hcl(length(unique(cut_avg)),c=90)
    for(i in unique(cut_c_cbic_v)){
      dir.create(paste(c_number,'_cluster/','CBIC/cluster_',i,sep = ''))
      if (sum(cut_c_cbic_v==i)==1) {
        after_scale = t(scale(d_cbic_l[[net]][cut_c_cbic_v==i, ]))
      } else{
        after_scale = t(scale(t(d_cbic_l[[net]][cut_c_cbic_v==i, ])))
      }
      mean = apply(after_scale, 2, mean, na.rm = T)
      sd = apply(after_scale, 2, sd, na.rm = T)
      df = data.frame(age = age_vec, mean = mean, sd = sd)
      p = ggplot(df, aes(x=age, y=mean)) + 
        geom_line(col=col_all[i]) +
        geom_point(col=col_all[i])+
        geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.1,alpha = 0.1, col=col_all[i])+
        ylab('value')+labs(title = paste('cluster',i,'\n N =',  sum(cut_c_cbic_v==i)))+
        theme(axis.text=element_text(size=6), axis.title=element_text(size=6))
      plot_list[[3*i]]=p
      write(rownames(d_cbic_l[[net]])[cut_c_cbic_v==i],paste(c_number,'_cluster/','CBIC/cluster_',i,'/variable.txt',sep = ''))
    }
    
    
    
    multiplot(plotlist = plot_list, cols = c_number)  
    
    
  }
  
  dev.off()


}



















