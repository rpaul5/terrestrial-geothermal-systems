# R script used for generating heatmaps

# library 
library(ggplot2)
library(tidyverse)
library(scales)
library(readxl)
library(dplyr)
library(cluster)
library(ggdendro)
library(dendextend)

# making the heatmap for EC# 3
signalp_norm_site_no_ds <- as.dendrogram(
  hclust(d= dist(cor(new_sig_ec_3[,-1], method = "spearman")))) # cluster the sites together
plot(signalp_norm_site_no_ds) 

new_signalp_norm_site_order <- order.dendrogram(signalp_norm_site_no_ds) # reorder df to match dendrogram

new_signalp_norm_cor_site <- data.frame(new_sig_ec_3)[,c(1,(new_signalp_norm_site_order+length(1)))] # create df with reordered values

new_signalp_norm_cor_site_2 <- new_signalp_norm_cor_site %>% column_to_rownames(., "EC") # copy the df
new_signalp_norm_cor_site_3  <- as.data.frame(t(new_signalp_norm_cor_site_2)) # transpose the copied df

new_signalp_norm_hc <- as.dendrogram(
  hclust(d= dist(cor(new_signalp_norm_cor_site_3, method = "spearman")))) # create hierarchical dendrogram for enzymes
plot(new_signalp_norm_hc)

new_signalp_norm_order <- order.dendrogram(new_signalp_norm_hc) # reorder the df to match second dendrogram


new_signalp_norm_new_long <- gather(new_signalp_norm_cor_site, site, coverage, 2:101, factor_key = TRUE) # gather the reordered df


new_signalp_norm_new_long$`EC` <- factor(x = new_signalp_norm_new_long$`EC` ,
                                         levels = new_signalp_norm_new_long$`EC`[new_signalp_norm_order], # match both ordered df together
                                         ordered = TRUE)

new_signalp_norm_plot_test <- ggplot(data=new_signalp_norm_new_long, aes(x=site, y=`EC`, fill=log(coverage+0.5))) + # generate the heatmap
  geom_tile(colour=NA) +
  coord_equal(-1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
  scale_x_discrete(position = "top") +
  theme(text = element_text(size = 8)) +
  scale_fill_gradientn(colours = c("black", "yellow", "red")) +
  theme(legend.position = "left") +
  scale_y_discrete(position = "right")

print(new_signalp_norm_plot_test) 


# coloring the dendrogram 
# import classifying files for province as assembly_province


df1 <- as.data.frame(t(assembly_province))                          # transpose df 
df2 <- tibble::rownames_to_column(df1, "EC#")                 # make col 1  into row names

new_order_site <- data.frame(df2)[,c(1,(signalp_norm_cazy_site_order+length(1)))]     # order the sites based off dendrogram
new_df <- as.data.frame(t(new_order_site))                                         # transpose again
new_df_2 <- janitor::row_to_names(new_df, 1, remove_rows_above = TRUE)            # make row 1 into col names      


color = rep(NA, length=length(new_df_2$Province))                                # create a vector and assign colors
color[which(new_df_2$Province=="Argentina backarc")] = "purple"
color[which(new_df_2$Province=="Active volcanic arc")] = "dodgerblue"
color[which(new_df_2$Province=="Active volcanic backarc")] = "blue4"
color[which(new_df_2$Province=="Costa Rica outer forearc")] = "darkturquoise"
color[which(new_df_2$Province=="Cordillera Talamanca")] = "chartreuse3"
color[which(new_df_2$Province=="Panama slab window")] = "goldenrod"
color[which(new_df_2$Province=="Spreading center hot spot")] = "red"

signalp_norm_cazy_site <- as.dendrogram(
  hclust(d= dist(cor(cazy_spread_75[,-1], method = "spearman"))))
# run dend again

labels_colors(signalp_norm_cazy_site) <- color                          # add colors

plot(signalp_norm_cazy_site) 
# plot the dend


# attempting to add the shapes
new_df_2$Type = as.numeric(as.factor(new_df_2$Type))                # create the column into numeric
new_df_2$Type[new_df_2$Type == 3] <- 1                              # changed the shapes i dont like the pch 3 shape

pch=c(1:max(new_df_2$Type))                                        # sets pch to be that column
nodes=pch[new_df_2$Type]                                           # creates the vector 

#nodePar = list(lab.cex = 0.6, pch = c(21,22),cex = 0.7) #node parameters

signalp_norm_cazy_site_2 = signalp_norm_cazy_site %>% set("leaves_pch", c(nodes))              # set leaves based on the vector

#par(mar=c(3,4,1,15))                      # not sure what this does
plot(signalp_norm_cazy_site_2, cex = 12)                       # plot


# this is the CAZyme family heatmap dataframes 
# the same exact steps used for the EC3 heatmap were implemented

signalp_norm_cazy <- as.dendrogram(
  hclust(d= dist(cor(cazy_spread_75[,-1], method = "spearman"))))
plot(signalp_norm_cazy, cex = 10) 

signalp_norm_cazy_site_order <- order.dendrogram(signalp_norm_cazy)

new_signalp_cazy_cor_site <- data.frame(cazy_spread_75)[,c(1,(signalp_norm_cazy_site_order+length(1)))]

new_signalp_cazy_cor_site_2 <- new_signalp_cazy_cor_site %>% column_to_rownames(., "DIAMOND")
new_signalp_cazy_cor_site_3  <- as.data.frame(t(new_signalp_cazy_cor_site_2))

new_signalp_cazy_hc <- as.dendrogram(
  hclust(d= dist(cor(new_signalp_cazy_cor_site_3, method = "spearman"))))
plot(new_signalp_cazy_hc, cex = 10)

new_signalp_cazy_order <- order.dendrogram(new_signalp_cazy_hc)

new_signalp_cazy_new_long <- gather(new_signalp_cazy_cor_site, site, coverage, 2:101, factor_key = TRUE)

new_signalp_cazy_new_long$DIAMOND <- factor(x = new_signalp_cazy_new_long$DIAMOND ,
                                            levels = new_signalp_cazy_new_long$DIAMOND[new_signalp_cazy_order],
                                            ordered = TRUE)
new_signalp_cazy_new_plot_test <- ggplot(data=new_signalp_cazy_new_long, aes(x=site, y=DIAMOND, fill=log(coverage+0.5))) +
  geom_tile(colour=NA) +
  coord_equal(-1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
  scale_x_discrete(position = "top") +
  theme(text = element_text(size = 12)) +
  scale_fill_gradientn(colours = c("black", "yellow", "red")) +
  theme(legend.position = "left") +
  scale_y_discrete(position = "right")

print(new_signalp_cazy_new_plot_test)


# heatmap for the MEROPs annotations
# the same exact steps utilized for the EC3 heatmap are used for this heatmap


new_signalp_merop_site <- as.dendrogram(
  hclust(d= dist(cor(new_merops_75[,-1], method = "spearman"))))
plot(new_signalp_merop_site) 

new_signalp_merops_site_order <- order.dendrogram(new_signalp_merop_site)

new_signalp_merops_cor_site <- data.frame(new_merops_75)[,c(1,(new_signalp_merops_site_order+length(1)))]

new_signalp_merops_cor_site_2 <- new_signalp_merops_cor_site %>% column_to_rownames(., "peptidase_family")
new_signalp_merops_cor_site_3  <- as.data.frame(t(new_signalp_merops_cor_site_2))

new_signalp_merops_hc <- as.dendrogram(
  hclust(d= dist(cor(new_signalp_merops_cor_site_3, method = "spearman"))))
plot(new_signalp_merops_hc)

new_signalp_merops_order <- order.dendrogram(new_signalp_merops_hc)


new_signalp_merops_new_long <- gather(new_signalp_merops_cor_site, site, coverage, 2:101, factor_key = TRUE)


new_signalp_merops_new_long$peptidase_family <- factor(x = new_signalp_merops_new_long$peptidase_family ,
                                                       levels = new_signalp_merops_new_long$peptidase_family[new_signalp_merops_order],
                                                       ordered = TRUE)
new_signalp_merops_plot_test <- ggplot(data=new_signalp_merops_new_long, aes(x=site, y=peptidase_family, fill=log(coverage+0.5))) +
  geom_tile(colour=NA) +
  coord_equal(-1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
  scale_x_discrete(position = "top") +
  theme(text = element_text(size = 8)) +
  scale_fill_gradientn(colours = c("black", "yellow", "red")) +
  theme(legend.position = "left") +
  scale_y_discrete(position = "right")

print(new_signalp_merops_plot_test)


new_signalp_merops_cor_site_2 <- new_signalp_merops_cor_site %>% column_to_rownames(., "peptidase_family")
new_signalp_merops_cor_site_3  <- as.data.frame(t(new_signalp_merops_cor_site_2))

new_signalp_merops_hc <- as.dendrogram(
  hclust(d= dist(cor(new_signalp_merops_cor_site_3, method = "pearson"))))
plot(new_signalp_merops_hc)

new_signalp_merops_order <- order.dendrogram(new_signalp_merops_hc)


new_signalp_merops_new_long <- gather(new_signalp_merops_cor_site, site, coverage, 2:101, factor_key = TRUE)


new_signalp_merops_new_long$peptidase_family <- factor(x = new_signalp_merops_new_long$peptidase_family ,
                                                       levels = new_signalp_merops_new_long$peptidase_family[new_signalp_merops_order],
                                                       ordered = TRUE)
new_signalp_merops_plot_test <- ggplot(data=new_signalp_merops_new_long, aes(x=site, y=peptidase_family, fill=log(coverage+0.5))) +
  geom_tile(colour=NA) +
  coord_equal(-1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
  scale_x_discrete(position = "top") +
  theme(text = element_text(size = 8)) +
  scale_fill_gradientn(colours = c("black", "yellow", "red")) +
  theme(legend.position = "left") +
  scale_y_discrete(position = "right")

print(new_signalp_merops_plot_test)

