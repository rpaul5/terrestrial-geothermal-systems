# this is the new pca plot script with the updated temperatures
library(tidyverse)
library(vegan)
library(devtools)
library(ggbiplot)
library(basetheme)
library(plotly)
library(ggfortify)


# Generating PCA plot for MEROPs annotations

merops_75_pca<- read.csv("merops_75.csv", # importing the dataframe to be used
                     row.names = 1)


merop.pca <- prcomp(merops_75_pca, center = TRUE, scale. = TRUE) # Performing the PCA analysis
p_merop <- autoplot(merop.pca, data = merops_75_pca, colour = color_mer , shape = FALSE, label.size = 5,                      #plotting the PCA analysis
                    loadings = TRUE, loadings.label = TRUE, loadings.colour = 'black', loadings.label.colour ="black",
                    loadings.label.size = 9, label.repel=TRUE, loadings.size = 1.5) + theme_classic()

plot(p_merop)



color_mer = rep(NA, length=length(merops_75_pca$Province))                                # create a vector and assign colors
color_mer[which(merops_75_pca$Province=="3")] = "purple"
color_mer[which(merops_75_pca$Province=="1")] = "dodgerblue"
color_mer[which(merops_75_pca$Province=="2")] = "blue4"
color_mer[which(merops_75_pca$Province=="5")] = "darkturquoise"
color_mer[which(merops_75_pca$Province=="4")] = "chartreuse3"
color_mer[which(merops_75_pca$Province=="6")] = "goldenrod"
color_mer[which(merops_75_pca$Province=="7")] = "red"



# cazy family pca plot 

cazy_75_pca<- read.csv("cazy_75_pca.csv",
                   row.names = 1)


color_caz = rep(NA, length=length(cazy_75_pca$Province))                                # create a vector and assign colors
color_caz[which(cazy_75_pca$Province=="3")] = "purple"
color_caz[which(cazy_75_pca$Province=="1")] = "dodgerblue"
color_caz[which(cazy_75_pca$Province=="2")] = "blue4"
color_caz[which(cazy_75_pca$Province=="5")] = "darkturquoise"
color_caz[which(cazy_75_pca$Province=="4")] = "chartreuse3"
color_caz[which(cazy_75_pca$Province=="6")] = "goldenrod"
color_caz[which(cazy_75_pca$Province=="7")] = "red"

cazy.pca <- prcomp(cazy_75_pca, center = TRUE,scale. = TRUE,)
p_cazy <- autoplot(cazy.pca, data = cazy_75_pca, colour = color_caz , shape = FALSE, label.size = 5, 
                   loadings = TRUE, loadings.label = TRUE, loadings.colour = 'black', loadings.label.colour ="black",
                   loadings.label.size = 5, label.repel=TRUE, loadings.size = 1.5) + theme_classic()

plot(p_cazy)


# enzyme commission 3 pca plot 

ec_75_pca <- read.csv("enzyme_com_3.csv", row.names = 1)

color_ec = rep(NA, length=length(ec_75_pca$Province))                                # create a vector and assign colors
color_ec[which(ec_75_pca$Province=="3")] = "purple"
color_ec[which(ec_75_pca$Province=="1")] = "dodgerblue"
color_ec[which(ec_75_pca$Province=="2")] = "blue4"
color_ec[which(ec_75_pca$Province=="5")] = "darkturquoise"
color_ec[which(ec_75_pca$Province=="4")] = "chartreuse3"
color_ec[which(ec_75_pca$Province=="6")] = "goldenrod"
color_ec[which(ec_75_pca$Province=="7")] = "red"


ec.pca <- prcomp(ec_75_pca, center = TRUE,scale. = TRUE,)
p_ec <- autoplot(ec.pca, data = ec_75_pca, colour = color_ec , shape = FALSE, label.size = 5, 
                 loadings = TRUE, loadings.label = TRUE, loadings.colour = 'black', loadings.label.colour ="black",
                 loadings.label.size = 5, label.repel=TRUE, loadings.size = 1.5) + theme_classic()

plot(p_ec)



# just chitin and peptidoglycan 

chitin_peptid <- read.csv("chitin_peptidoglycan_pca.csv",
                          row.names = 1)

color_pep = rep(NA, length=length(chitin_peptid$Province))                                # create a vector and assign colors
color_pep[which(chitin_peptid$Province=="3")] = "purple"
color_pep[which(chitin_peptid$Province=="1")] = "dodgerblue"
color_pep[which(chitin_peptid$Province=="2")] = "blue4"
color_pep[which(chitin_peptid$Province=="5")] = "darkturquoise"
color_pep[which(chitin_peptid$Province=="4")] = "chartreuse3"
color_pep[which(chitin_peptid$Province=="6")] = "goldenrod"
color_pep[which(chitin_peptid$Province=="7")] = "red"

chi_pep.pca <- prcomp(chitin_peptid, center = TRUE,scale. = TRUE,)
p_chi_pep <- autoplot(chi_pep.pca, data = chitin_peptid, colour = color_pep , shape = FALSE, label.size = 5, 
                      loadings = TRUE, loadings.label = TRUE, loadings.colour = 'black', loadings.label.colour ="black",
                      loadings.label.size = 5, label.repel=TRUE, loadings.size = 1.5) + theme_classic()

plot(p_chi_pep)

# pca plot for xylan associated enzymes
xylan <- read.csv("xylan_pca.csv",
                      row.names = 1)

color_xylan = rep(NA, length=length(xylan$Province))                                # create a vector and assign colors
color_xylan[which(xylan$Province=="3")] = "purple"
color_xylan[which(xylan$Province=="1")] = "dodgerblue"
color_xylan[which(xylan$Province=="2")] = "blue4"
color_xylan[which(xylan$Province=="5")] = "darkturquoise"
color_xylan[which(xylan$Province=="4")] = "chartreuse3"
color_xylan[which(xylan$Province=="6")] = "goldenrod"
color_xylan[which(xylan$Province=="7")] = "red"

xylan_2.pca <- prcomp(xylan, center = TRUE,scale. = TRUE)

p_xylan_2 <- autoplot(xylan_2.pca, data = xylan, colour = color_xylan , shape = FALSE, label.size = 5, 
                      loadings = TRUE, loadings.label = TRUE, loadings.colour = 'black', loadings.label.colour ="black",
                      loadings.label.size = 5, label.repel=TRUE, loadings.size = 1.5) + theme_classic()

plot(p_xylan_2)

