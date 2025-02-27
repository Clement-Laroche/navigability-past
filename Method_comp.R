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

load(file = "results/Indic_partsection_wo_constraint.Rdata") # file created by the partsection_wo_constraint.R script
tab_method1 <- tab_indic
load(file = "results/Indic_fullsection_wo_constraint.Rdata") # file created by the fullsection_wo_constraint.R script
tab_method2 <- tab_indic
load(file = "results/Indic_fullsection_w_constraint.Rdata") # file created by the fullsection_w_constraint.R script
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

load("results/fullsection_w_constraint.Rdata") # file produced by the fullsection_w_constraint.R script
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

# plot of the clustering results stored in Figure 11 of the paper

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

ggsave(filename = "Fig_11.pdf",plot = G,device = "pdf",width = 19,height = 12)

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

load("results/partsection_wo_constraint.Rdata")
tab_dist <- table_points
tab_dist <- st_as_sf(x = tab_dist,coords = c(1,2))
st_crs(tab_dist) <- 4326
tab_dist <- st_transform(tab_dist,3035) 
load("results/fullsection_w_constraint.Rdata")
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
write.csv2(x = r_l,file = "results/navi_sec_dist.csv")

