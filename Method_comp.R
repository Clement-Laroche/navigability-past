library(udunits2)
library(sf)
library(terra)
library(gridExtra)
library(tidyr)
library(dtw)
library(raster)
library(xtable)
library(rlist)
library(osmdata) 
library(transport)
library(data.table)
library(fda)
library(fdacluster)
library(ggplot2)
library(mcp)
library(changepoint)
library(segmented)
library(patchwork)
library(dplyr)
library(FactoMineR)
library(ggmap)
library(longitudinalData)
library(rPref)
library(maotai)

################################################################
###################### INDICATORS OF PERF ######################
######################     COMPARISON     ######################
################################################################


# grouping all results indicators into the same table

load(file = "Indic_partsection_wo_constraint.Rdata") # file created by the partsection_wo_constraint.R script
tab_method1 <- tab_indic
load(file = "Indic_fullsection_wo_constraint.Rdata") # file created by the fullsection_wo_constraint.R script
tab_method2 <- tab_indic
load(file = "Indic_fullsection_w_constraint.Rdata") # file created by the fullsection_w_constraint.R script
tab_method3 <- tab_indic
rm(tab_indic)
tab_method3 <- tab_method3 %>% select(-c("Linearity_NCP","Linearity_res"))
tab_method1$meth <- 1
tab_method2$meth <- 2
tab_method3$meth <- 3
tab <- rbind(tab_method1,tab_method2,tab_method3)

# minor corrections and adding two new columns on the table
tab$Pts_order_M1[tab$Pts_order_M1=="CP1->Low_emp->Up_emp"] <- "CP1->Up_emp->Low_emp"
tab$Pts_order_M2[tab$Pts_order_M2=="Up_CP2->Low_emp->Up_emp->Low_CP2"] <- "Up_CP2->Up_emp->Low_emp->Low_CP2"
tab$Prec_up2up_M2_rel <- tab$Prec_up2up_M2/unlist(tab[["char_max_length"]])
tab$Prec_low2low_M2_rel <- tab$Prec_low2low_M2/unlist(tab[["char_max_length"]])
rm(tab_method1,tab_method2,tab_method3)

# grouping and computing statistical results of each method

global_indic <- tab %>% group_by(meth) %>% summarize(CP1_prec_mean = mean(Prec_2up_M1),
                                   CP1_prec_med = median(Prec_2up_M1),
                                   CP2_prec_mean_up = mean(Prec_up2up_M2),
                                   CP2_prec_med_up = median(Prec_up2up_M2),
                                   CP2_prec_mean_low = mean(Prec_low2low_M2),
                                   CP2_prec_med_low = median(Prec_low2low_M2),
                                   CP2_prec_mean_up_rel = mean(Prec_up2up_M2_rel),
                                   CP2_prec_med_up_rel = median(Prec_up2up_M2_rel),
                                   CP2_prec_mean_low_rel = mean(Prec_low2low_M2_rel),
                                   CP2_prec_med_low_rel = median(Prec_low2low_M2_rel),
                                   ideal_order1 = sum(Pts_order_M1 == "Up_emp->CP1->Low_emp"),
                                   ideal_order2 = sum(Pts_order_M2 == "Up_CP2->Up_emp->Low_emp->Low_CP2"),
                                   Mean_emp_over_est = mean(Overlap_emp_over_est_M2),
                                   Mean_est_over_emp = mean(Overlap_est_over_emp_M2),
                                   Med_emp_over_est = median(Overlap_emp_over_est_M2),
                                   Med_est_over_emp = median(Overlap_est_over_emp_M2))

# preparing the data needed for the clustering

load("fullsection_w_constraint.Rdata") # file produced by the fullsection_w_constraint.R script
L_key <- which(unique(RES$rivers) %in% comp_points$Nom.riv)
river_names_cal <- unique(RES$rivers)[L_key]
river_names_all <- unique(RES$rivers)

# creating a table with elevation points of each river located upstream of the highest confluence point 

ntab <- table_points %>% filter(name == river_names_all[1],Z >= conf_points$Z[conf_points$Name == river_names_all[1]])
for(i in river_names_all[-1])
{
  if(i %in% conf_points$Name)
  {
    temp <- table_points %>% filter(name == i,Z >= conf_points$Z[conf_points$Name == i])
    ntab <- rbind(ntab,temp)
  }else
  {
    temp <- table_points %>% filter(name == i)
    ntab <- rbind(ntab,temp)
  }
}
rm(temp)

# smoothing step of the elevation profiles 

SMR <- as.list(rep(NA,length(river_names_all))) #L_key
t_t <- rep(NA,length(river_names_all))
for(i in 1:length(SMR))
{
    pos <- which(ntab$name == river_names_all[i])
    SMR[[i]] <- smooth.spline(x = seq(0,length(pos)-1,length.out = length(pos)),y = ntab$Z[pos],df = 58)
    t_t[i] <- max(SMR[[i]]$y)/max(ntab$Z[pos])
    # SMR[[i]] <- smooth.spline(x = seq(0,1,length.out = length(pos)),y = (ntab$Z[pos]-min(ntab$Z[pos]))/(max(ntab$Z[pos])-min(ntab$Z[pos])),df = 58)
    SMR[[i]]$y <- (SMR[[i]]$y-min(SMR[[i]]$y))/(max(SMR[[i]]$y)-min(SMR[[i]]$y))
}

# preparing the elvation profiles that we are going to cluster and
# storing them into a matrix

names(SMR) <- river_names_all
M_traj <- matrix(NA,nrow = length(river_names_all),ncol = 1000) 
for(i in 1:length(river_names_all)) 
{
  pos <- which(ntab$name == river_names_all[i])
  M_traj[i,] <- SMR[[which(names(SMR) == river_names_all[i])]]$y[round(seq(1,length(pos),length.out = 1000))]
}

# performing the clustering and formatting the results

res_clust <- fdahclust(x = seq(0,1,length.out = 1000),
                       y = M_traj,
                       n_clusters = 5,
                       warping_class = "affine",
                       cluster_on_phase = TRUE,
                       linkage_criterion = "ward.D2")
summary(as.factor(res_clust$memberships))
names(res_clust$memberships) <- river_names_all
summary(as.factor(res_clust$memberships[names(res_clust$memberships) %in% river_names_cal]))

# preparing data to plot the clustering results

orig_curve <- as.data.frame(apply(res_clust$original_curves, 1, c))
orig_grid <- as.data.frame(apply(res_clust$original_grids, 1, c))
align_grid <- as.data.frame(apply(res_clust$aligned_grids, 1, c))
colnames(orig_curve) <- river_names_all
colnames(orig_grid) <- river_names_all
colnames(align_grid) <- river_names_all
orig_curve <- orig_curve %>% 
  pivot_longer(cols = 1:ncol(orig_curve)) %>% 
  arrange(name)
orig_grid <- orig_grid %>% 
  pivot_longer(cols = 1:ncol(orig_grid)) %>% 
  arrange(name)
align_grid <- align_grid %>% 
  pivot_longer(cols = 1:ncol(align_grid)) %>% 
  arrange(name)
colnames(orig_curve)[2] <- "Elevation"
colnames(orig_grid)[2] <- "Point"
colnames(align_grid)[2] <- "Align_point"
orig_curve <- orig_curve %>% mutate(point = orig_grid$Point,
                      aligned_point = align_grid$Align_point)
orig_curve$clust_member <- rep(res_clust$memberships,each = 1000)

# plot of the clustering results

G1 <- ggplot(data = orig_curve)+
  geom_line(aes(x = point,
                y = Elevation,
                group = name,
                col = as.factor(clust_member)))+
  ylab("Elevation")+xlab("Original distance grid")+
  scale_color_discrete("Cluster")+theme(legend.position = "bottom")

G2 <- ggplot(data = orig_curve)+
  geom_line(aes(x = aligned_point,
                y = Elevation,
                group = name,
                col = as.factor(clust_member)))+
  ylab("Elevation")+xlab("Transformed distance grid")+
  scale_color_discrete("Cluster")+theme(legend.position = "bottom")

G <- (G1/G2)+
  plot_layout(guides = "collect")&theme(legend.position = "bottom",text = element_text(size = 16))

ggsave(filename = "clustering.pdf",plot = G,device = "pdf",width = 19,height = 12)

# preparing the table with the indicators statistics to compare the methods across the clusters

col_clust <- res_clust$memberships[names(res_clust$memberships) %in% river_names_cal]
tab$clust_shape <- rep(col_clust,3)
h_tab <- tab %>% filter(Name %in% river_names_cal & meth == 3)
h_tab <- h_tab %>% select(Name,Prec_up2up_M2,Prec_low2low_M2,Pts_order_M2,Pts_order_M2,char_max_length,clust_shape)
h_tab$Prec_up2up_M2_rel <- round(h_tab$Prec_up2up_M2/unlist(h_tab$char_max_length)*100,2)
h_tab$Prec_low2low_M2_rel <- round(h_tab$Prec_low2low_M2/unlist(h_tab$char_max_length)*100,2)
h_tab$Prez_up <- paste(as.character(h_tab$Prec_up2up_M2)," (",as.character(h_tab$Prec_up2up_M2_rel)," %)", sep = "")
h_tab$Prez_down <- paste(as.character(h_tab$Prec_low2low_M2)," (",as.character(h_tab$Prec_low2low_M2_rel)," %)", sep = "")
h_tab <- h_tab %>% select(Name,Prez_up,Prez_down,Pts_order_M2,clust_shape)
h_tab <- h_tab[order(h_tab$clust_shape),]

# xtable(h_tab)

shape_indic <- tab %>% group_by(meth,clust_shape) %>% summarize(#CP1_prec_mean = mean(Prec_2up_M1),
                                                     #CP1_prec_med = median(Prec_2up_M1),
                                                     CP2_prec_mean_up = mean(Prec_up2up_M2),
                                                     CP2_prec_med_up = median(Prec_up2up_M2),
                                                     CP2_prec_mean_low = mean(Prec_low2low_M2),
                                                     CP2_prec_med_low = median(Prec_low2low_M2),
                                                     CP2_prec_mean_up_rel = mean(Prec_up2up_M2_rel),
                                                     CP2_prec_med_up_rel = median(Prec_up2up_M2_rel),
                                                     CP2_prec_mean_low_rel = mean(Prec_low2low_M2_rel),
                                                     CP2_prec_med_low_rel = median(Prec_low2low_M2_rel),
                                                     # ideal_order1 = sum(Pts_order_M1 == "Up_emp->CP1->Low_emp"),
                                                     # ideal_order2 = sum(Pts_order_M2 == "Up_CP2->Up_emp->Low_emp->Low_CP2"),
                                                     Mean_emp_over_est = mean(Overlap_emp_over_est_M2),
                                                     Mean_est_over_emp = mean(Overlap_est_over_emp_M2),
                                                     Med_emp_over_est = median(Overlap_emp_over_est_M2),
                                                     Med_est_over_emp = median(Overlap_est_over_emp_M2))


# plotting the indicators statistics across the clusters

G1 <- ggplot(data = tab,aes(x = as.factor(clust_shape),y = Prec_up2up_M2))+geom_boxplot(aes(fill = as.factor(meth)))+
  xlab("Shape cluster")+
  ylab("Precision upper point")+
  scale_fill_discrete("Method")

G2 <- ggplot(data = tab,aes(x = as.factor(clust_shape),y = Prec_low2low_M2))+geom_boxplot(aes(fill = as.factor(meth)))+
  xlab("Shape cluster")+
  ylab("Precision lower point")+
  scale_fill_discrete("Method")

g1 <- ggplot(data = shape_indic)+geom_line(aes(x = clust_shape,y = CP2_prec_mean_up,group = meth,col = as.factor(meth)))+
  xlab("Cluster")+ylab("Mean upper precision (m)")+scale_color_discrete("Method")+
  ggtitle("(a)")
g2 <- ggplot(data = shape_indic)+geom_line(aes(x = clust_shape,y = CP2_prec_mean_low,group = meth,col = as.factor(meth)))+
  xlab("Cluster")+ylab("Mean lower precision (m)")+scale_color_discrete("Method")+
  ggtitle("(b)")
g1/g2

g1_ter <- ggplot(data = shape_indic)+geom_line(aes(x = clust_shape,y = Mean_est_over_emp,group = meth,col = as.factor(meth)))+
  xlab("Cluster")+ylab("Empirical coverage (%)")+
  scale_color_discrete("Method")+ylim(c(5,100))+
  ggtitle("(c)")
g2_ter <- ggplot(data = shape_indic)+geom_line(aes(x = clust_shape,y = Mean_emp_over_est,group = meth,col = as.factor(meth)))+
  xlab("Cluster")+ylab("Estimated coverage (%)")+
  scale_color_discrete("Method")+ylim(c(5,100))+
  ggtitle("(d)")

g1_quad <- ggplot(data = shape_indic)+geom_line(aes(x = clust_shape,y = CP2_prec_mean_up_rel*100,group = meth,col = as.factor(meth)))+
  xlab("Cluster")+ylab("Mean relative\nupper precision (%)")+scale_color_discrete("Method")+ylim(c(0,50))+
  ggtitle("(e)")
g2_quad <- ggplot(data = shape_indic)+geom_line(aes(x = clust_shape,y = CP2_prec_mean_low_rel*100,group = meth,col = as.factor(meth)))+
  xlab("Cluster")+ylab("Mean relative\nlower precision (%)")+scale_color_discrete("Method")+ylim(c(0,50))+
  ggtitle("(f)")
G <- ((g1|g2)/(g1_quad|g2_quad)/(g1_ter|g2_ter))+
  plot_layout(guides = "collect")&theme(legend.position = "bottom",text = element_text(size = 8))

ggsave(filename = "clustering_results.pdf",plot = G,device = "pdf",width = 19,height = 12,units = "cm")

# table comparing the length of the navigable section estimated by the model 
# to the length of the empirical navigable section

load("partsection_wo_constraint.Rdata")
tab_dist <- table_points
tab_dist <- st_as_sf(x = tab_dist,coords = c(1,2))
st_crs(tab_dist) <- 4326
tab_dist <- st_transform(tab_dist,3035) 
load("fullsection_w_constraint.Rdata")
RES <- RES %>% filter(stat == "Est",model_type == "M2",id == "Lower_point")
RES <- st_zm(x = RES,drop = TRUE)
comp_points <- comp_points %>% filter(Type == "Low")
comp_points <- st_transform(comp_points,3035)
r_l <- data.frame(matrix(NA,ncol = 3,nrow = 18))
colnames(r_l) <- c("Name","l_emp","l_est")
r_l$Name <- river_names_cal
for(i in r_l$Name)
{
  pos <- which(tab_dist$name == i)
  p_e <- comp_points[comp_points$Nom.riv == i,]
  p_est <- RES[RES$rivers == i,]
  l_emp <- which.min(st_distance(x = tab_dist$geometry[pos],y = p_e$geometry))
  l_est <- which.min(st_distance(x = tab_dist$geometry[pos],y = p_est$geometry))
  r_l$l_emp[which(r_l$Name == i)] <- (length(pos)-l_emp)*100
  r_l$l_est[which(r_l$Name == i)] <- (length(pos)-l_est)*100
}
write.csv2(x = r_l,file = "navi_sec_dist.csv")




#################################################################
######################## DRAFT CODE #############################  
#################################################################

# sil <- rep(NA,8)
# for(i in 2:8)
# {
#   res_clust <- fdahclust(x = seq(0,1,length.out = 1000),
#                          y = M_traj,
#                          n_clusters = i,
#                          warping_class = "affine",
#                          cluster_on_phase = TRUE,
#                          linkage_criterion = "ward.D2")
#   sil[i] <- mean(res_clust$silhouettes)
# }


# g1_bis <- ggplot(data = shape_indic)+geom_line(aes(x = clust_shape,y = CP2_prec_med_up,group = meth,col = as.factor(meth)))+
#   xlab("Cluster")+ylab("Median upper precision (m)")+scale_color_discrete("Method")+ylim(c(27500,200000))+
#   ggtitle("(c)")
# g2_bis <- ggplot(data = shape_indic)+geom_line(aes(x = clust_shape,y = CP2_prec_med_low,group = meth,col = as.factor(meth)))+
#   xlab("Cluster")+ylab("Median lower precision (m)")+scale_color_discrete("Method")+ylim(c(27500,200000))+
#   ggtitle("(d)")
# g1_bis/g2_bis
# /(g1_bis|g2_bis)


# ntab_clust_fig <- data.table(matrix(NA,nrow = length(river_names)*1000,ncol = 3))
# ntab_clust_fig$V1 <- rep(river_names,each = 1000)
# ntab_clust_fig$V2 <- rep(1:1000,length(river_names))
# for(i in 1:nrow(M_traj))
# {
#   ntab_clust_fig$V3[(1:1000)+(i-1)*1000] <- M_traj[i,] 
# }
# colnames(ntab_clust_fig) <- c("name","x","y")
# ggplot(data = ntab_clust_fig[ntab_clust_fig$name %in% names(res_clust$memberships[res_clust$memberships==5]),])+
#   geom_line(aes(x = x,y = y,group = name,color = name),linewidth = 0.5)

# g1_ter <- ggplot(shape_indic[shape_indic$meth!=1,])+
#   geom_point(aes(x = Mean_emp_over_est,
#                  y = Mean_est_over_emp,
#                  col = as.factor(clust_shape)))+
#   scale_color_discrete("Shape cluster")+xlab("Emp over Est")+ylab("Est over Emp")+theme(legend.position = "bottom")
# g2_ter <- ggplot(shape_indic[shape_indic$meth!=1,])+
#   geom_point(aes(x = Mean_emp_over_est,
#                  y = Mean_est_over_emp,
#                  col = as.factor(meth)))+
#   scale_color_discrete("Method")+xlab("Emp over Est")+ylab("Est over Emp")+theme(legend.position = "bottom")
# g1_ter|g2_ter

# ggplot(ntab[name %in% river_names[which(res_clust$memberships == 5)],])+geom_line(aes(x = D,y = Z, group = name,col = name))









# ECDF <- as.list(rep(NA,length(river_names)))#L_key
# for(i in 1:length(ECDF))
# {
#   # pos <- which(table_points$name == river_names[i])
#   ECDF[[i]] <- ecdf(SMR[[i]]$y)
# # ecdf()
# }
# n_per_river_conf <- unlist((ntab %>% group_by(name) %>% summarize(n()))[,2])
# n_per_river <- unlist((table_points[(name%in%river_names)&(!name%in%conf_points$Name),] %>% group_by(name) %>% summarize(n()))[,2])
# n_min <- min(c(n_per_river_conf,n_per_river))
# n_mean <- mean(c(n_per_river_conf,n_per_river))


# comp_fromr <- tab[tab$meth == 1,c("Name","Prec_low2low_M2","Prec_up2up_M2")]
# p <- high(comp_fromr$Prec_low2low_M2) * high(comp_fromr$Prec_up2up_M2)
# comp_fromr <- psel(comp_fromr, p, top = nrow(comp_fromr))
# comp_fromr$.level <- as.factor(comp_fromr$.level)
# ggplot(data = comp_fromr, aes(x= Prec_low2low_M2, y= Prec_up2up_M2, colour=.level, label=Name))+
#   geom_point()+
#   geom_text(aes(label = Name))
#   
# 
# 
# 
# 
# 
# # load("fullsection_wo_constraint.Rdata")
# # load("partsection_wo_constraint.Rdata")
# 
# tab_indic$Name[which(res_clust$cluster == 2)]
# 
# 
# 
# 
# 
# 
# ggplot(data = tab_indic)+geom_point(aes(x = unlist(char_max_length),
#                                         y = unlist(char_max_height),
#                                         col = as.factor(clust)))+
#   scale_color_discrete("")
# # different kmeans
# K1 <- kmeans(x = tab_indic[,c("char_max_length","char_max_height")],centers = 4)
# K1 <- K1$cluster
# K2 <- kmeans(x = tab_indic[,c("char_max_length","char_max_height","char_pos_conf")],centers = 4)
# K2 <- K2$cluster
# K3 <- kmeans(x = tab_indic[,c("char_max_length","char_max_height","char_pos_conf","char_AUC")],centers = 4)
# K3 <- K3$cluster  
# K4 <- kmeans(x = tab_indic[,c("char_max_length","char_max_height","char_pos_conf","char_AUC","Linearity_NCP")],centers = 4)
# K4 <- K4$cluster
# K5 <- kmeans(x = tab_indic[,c("char_max_length","char_max_height","char_pos_conf","char_AUC","Linearity_NCP","Linearity_res")],centers = 4)
# K5 <- K5$cluster
# tab_indic$Name[K3 == 1]
# 
# R1 <- kmeans(x = tab_indic[,c("Prec_2up_M1","Prec_2low_M1")],centers = 6)
# R1 <- R1$cluster
# R2 <- kmeans(x = tab_indic[,c("Prec_up2up_M2","Prec_low2low_M2")],centers = 6)
# R2 <- R2$cluster
# R3 <- kmeans(x = tab_indic[,c("Prec_2up_M1","Prec_2low_M1","Prec_up2up_M2","Prec_low2low_M2")],centers = 6)
# R3 <- R3$cluster
# R4 <- kmeans(x = tab_indic[,c("Prec_2up_M1","Prec_2low_M1","Prec_up2up_M2","Prec_low2low_M2","Overlap_est_over_emp_M2","Overlap_emp_over_est_M2")],centers = 6)
# R4 <- R4$cluster
# 
# 
# wasserstein_d <- function(x,y)
# {
#   cdf1 <- ecdf(x)
#   cdf2 <- ecdf(y)
#   abs_diff <- function(t)
#   {
#     abs(cdf1(t)-cdf2(t))
#   }
#   res <- integrate(f = abs_diff,lower = 0,upper = 1,subdivisions = 1000000)
#   return(res$value)
# }
# 
# ggplot(data = new_points,aes(x = V1,y = V2,label = Name))+
#   geom_point(aes(col = unlist(tab_indic$Pts_order_M2)))+
#   scale_color_discrete("Prec")+geom_text()
# 
# 
# 
# 
# 
# 
# p <- high(tab_indic$Prec_low2low_M2) * high(tab_indic$Prec_up2up_M2)
# D1 <- psel(tab_indic, p, top = nrow(tab_indic))
# D1$.level <- as.factor(D1$.level)
# 
# 
# 
# ggplot(data = D1,aes(x = Prec_low2low_M2,y = Prec_up2up_M2,label = Name))+
#   geom_point(aes(col = .level))+geom_text_repel()
# 
# 
# res_clust <- epmeans(elist = ECDF,k = 7)
# res_clust <- as.factor(res_clust$cluster)
# summary(res_clust)
# 
# M_dist <- matrix(NA,nrow = length(L_key),ncol = length(L_key))
# diag(x = M_dist) <- 0
# row.names(M_dist) <- river_names
# colnames(M_dist) <- river_names
# for(i in 1:(nrow(M_dist)-1))
# {
#   for(j in (i+1):nrow(M_dist))
#   {
#     pos_x <- which(table_points$name == row.names(M_dist)[i])
#     pos_y <- which(table_points$name == colnames(M_dist)[j])
#     ind_x <- round(seq(1,length(pos_x),length.out = n_min))
#     ind_y <- round(seq(1,length(pos_y),length.out = n_min))
#     F_x <- table_points$Z[pos_x[ind_x]]
#     F_y <- table_points$Z[pos_y[ind_y]]
#     # ind_x <- ind_x/max(ind_x)
#     # ind_y <- ind_y/max(ind_y)
#     M_dist[i,j] <- #wasserstein1d(a = F_x,b = F_y,p = 1)
#       dtw(x = F_x,y = F_y,distance.only = TRUE)$distance
#     # M_dist[i,j] <- distFrechet(Px = ind_x,Py = F_x,Qx = ind_y,Qy = F_y,timeScale = 0.01,FrechetSumOrMax = "max")
#     M_dist[j,i] <- M_dist[i,j]
#   }
# }
# new_points <- cmdscale(d = M_dist,k = 2)
# new_points <- as.data.frame(new_points)
# new_points$Name <- row.names(new_points)
# plot(new_points$V1,new_points$V2)
# # res_clust <- kmeans(x = new_points[,c("V1","V2")],centers = 4)
# # res_clust <- as.factor(res_clust$cluster)
# # summary(res_clust)
# res_clust <- hclust(d = M_dist,method = "average")
# plot(res_clust)
# summary(as.factor(cutree(res_clust, k = 6)))
# cutree(res_clust, k = 6)
# 
# try <- table_points[which(table_points$name %in% river_names),]
# try$clust <- tab$clust_shape[match(try$name,river_names)]
# try <- st_as_sf(x = try,coords = c("X","Y"))
# st_crs(try) <- 4326
# try <- st_transform(try,crs = 2154)
# ggplot(data = try)+geom_point(aes(x = st_coordinates(try)[,1],
#                                  y = st_coordinates(try)[,2],group = L1,col = as.factor(clust)))

# essai <- compare_caps(x = seq(0,1,length.out = 1000),
# y = M_traj,
# n_clusters = 1:7,
# metric = c("l2"),
# clustering_method = c("hclust-average"),
# warping_class = c("affine"),
# centroid_type = c("mean"),
# cluster_on_phase = TRUE)


# river_test <- c()
# river_dist <- c()
# for(i in river_names_all[which(!river_names_all %in% river_names_cal)])
# {
#   d <- rep(NA,5)
#   for(j in 1:5)
#   {
#     p <- length(SMR[[which(names(SMR) == i)]]$y)
#     x_c <- res_clust$center_grids[j,]
#     # pos <- which(x_c >= 0 & x_c <= 1)
#     # x_c <- x_c[pos]
#     y_c <- res_clust$center_curves[j,,]
#     # [pos]
#     # pos_x_t <- round(p*x_c)
#     x_t <- seq(0,1,length.out = 1000)
#     # (pos_x_t/p)
#     y_t <- SMR[[which(names(SMR) == i)]]$y[round(seq(1,p,length.out = 1000))]
#     # pos_x_t
#     d[j] <- fdadist(x = rbind(x_c,x_t),
#                   y = rbind(y_c,y_t),
#                   warping_class = "affine",cluster_on_phase = TRUE)
#   }
#   river_dist <- rbind(river_dist,d)
#   cluster <- which.min(d)
#   river_test <- c(river_test,cluster)
# }
# summary(as.factor(river_test))
# names(river_test) <- river_names_all[which(!river_names_all %in% river_names_cal)]
# clust_1 <- names(which(river_test == 3))
# C1 <- c()
# for(i in clust_1)
# {
#   p <- length(SMR[[which(names(SMR) == i)]]$y)
#   C1 <- rbind(C1,SMR[[which(names(SMR) == i)]]$y[round(seq(1,p,length.out = 1000))])
# }
# plot(res_clust$center_grids[1,],res_clust$center_curves[1,,],type = 'l',col = "red")
# lines(res_clust$center_grids[2,],res_clust$center_curves[2,,],type = 'l',col = "blue")
# lines(res_clust$center_grids[3,],res_clust$center_curves[3,,],type = 'l',col = "violet")
# lines(res_clust$center_grids[4,],res_clust$center_curves[4,,],type = 'l',col = "lightgreen")
# lines(res_clust$center_grids[5,],res_clust$center_curves[5,,],type = 'l',col = "lightblue")
# lines(seq(0,1,length.out = 1000),C1[1,])
# lines(seq(0,1,length.out = 1000),C1[2,])
# lines(seq(0,1,length.out = 1000),C1[3,])
# lines(seq(0,1,length.out = 1000),C1[4,])
# lines(seq(0,1,length.out = 1000),C1[5,])
# lines(seq(0,1,length.out = 1000),C1[6,])
# lines(seq(0,1,length.out = 1000),C1[7,])
# lines(seq(0,1,length.out = 1000),C1[8,])
# 
# dtwDist(mx = matrix(data = res_clust$center_curves[1,,],nrow = 1),my = matrix(SMR$Adda$y,nrow = 1),)
# dtwDist(mx = matrix(data = res_clust$center_curves[2,,],nrow = 1),my = matrix(SMR$Adda$y,nrow = 1))
# dtwDist(mx = matrix(data = res_clust$center_curves[3,,],nrow = 1),my = matrix(SMR$Adda$y,nrow = 1))
# dtwDist(mx = matrix(data = res_clust$center_curves[4,,],nrow = 1),my = matrix(SMR$Adda$y,nrow = 1))
# dtwDist(mx = matrix(data = res_clust$center_curves[5,,],nrow = 1),my = matrix(SMR$Adda$y,nrow = 1))

