# R script of analyses reported in
# Sexual selection and speciation: a meta-analysis of comparative studies
# by Tim Janicke*, Tamra C. Mendelson, Michael G. Ritchie, Lucas Marie-Orleach, Jeanne Tonnabel
# *Corresponding author (tim.janicke@cefe.cnrs.fr)

rm(list = ls())
print(.libPaths())
print(sessionInfo())
print(version)
citation()

library(data.table)
library(metafor)
library(ape)
library(ggplot2)
library(cowplot)
library(orchaRd)
library(outliers)
library(dplyr)
library(tidyr)
library(rphylopic)
library(ggtree)


#setwd()

# 1. Import data and check for outliers ####

Data <- read.csv("Meta_SexSelSpec_Data_v01.csv", header=TRUE, sep=",", na.strings="", dec=".", strip.white=TRUE)
Tree <- read.tree("Meta_SexSelSpec_Tree_v01.nwk")

Grubbs_test_r <- grubbs.test(Data$r, opposite = FALSE)
Grubbs_test_r

# Generate filtered dataset with more conservative inclusion criteria
Data_conservative <-  subset(Data, Obs_ID  != 16) # only identified outlier
Data_conservative <-  subset(Data_conservative, N_r  > 9)
Data_conservative <-  subset(Data_conservative, Study_ID != 13) # Exclude Misof 2002 because of high Cook's Distance


## ___ ####
# 2. Prune Tree and convert into correlation matrix ####

Tree_List <- sort(Tree$tip.label)
Species_List_ALL <- unique(Data$Taxon_Phylo)
Tree_pruned_ALL <-drop.tip(Tree, Tree$tip.label[-na.omit(match(Species_List_ALL, Tree$tip.label))])
Tree_List <- sort(Tree_pruned_ALL$tip.label)
Data_list <- sort(unique(Data$Taxon_Phylo))
PhyloTree <- as.matrix(forceSymmetric(vcv(Tree_pruned_ALL, corr=TRUE)))


## ___ ####
# 3. Meta-analysis  — Null model r ####

## * Traditional (not correcting for phylogeny), full dataset ####

REML_Model_All_Null_trad_r <- rma.mv(r, r_Var,
                                     mod =~1,
                                     random=list(~ 1 | Taxon, # non-phylo effect
                                                 ~ 1 | Study_ID, 
                                                 ~ 1 | Obs_ID), 
                                     method="REML",
                                     test="t", # using t dist rather than z 
                                     data=Data)
REML_Model_All_Null_trad_r

# extract model statistics for Table 1 
REML_Model_All_Null_phylo_trad_r_PredictionIntervals <- predict(REML_Model_All_Null_trad_r, digits=3)
W <- diag(1/REML_Model_All_Null_trad_r$vi)
X <- model.matrix(REML_Model_All_Null_trad_r)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
Model_All_Null_trad_r_I2_Total <- 100 * sum(REML_Model_All_Null_trad_r$sigma2) / (sum(REML_Model_All_Null_trad_r$sigma2) + (REML_Model_All_Null_trad_r$k-REML_Model_All_Null_trad_r$p)/sum(diag(P)))
Model_All_Null_trad_r_Sum_I2 <- sum(REML_Model_All_Null_trad_r$sigma2)
I2_Model_All_Null_trad_r_Taxon_Taxon    <- (REML_Model_All_Null_trad_r$sigma2[1]/Model_All_Null_trad_r_Sum_I2)*Model_All_Null_trad_r_I2_Total
I2_Model_All_Null_trad_r_Taxon_Study_ID <- (REML_Model_All_Null_trad_r$sigma2[2]/Model_All_Null_trad_r_Sum_I2)*Model_All_Null_trad_r_I2_Total
I2_Model_All_Null_trad_r_Taxon_Obs_ID   <- (REML_Model_All_Null_trad_r$sigma2[3]/Model_All_Null_trad_r_Sum_I2)*Model_All_Null_trad_r_I2_Total

Model_All_Null_trad_r_Table <- data.frame(model= "Traditional (full dataset)",
                                          n=length(unique(Data$Study_ID)),
                                          k=nrow(Data),
                                          r=round(REML_Model_All_Null_phylo_trad_r_PredictionIntervals[[1]],3),
                                          CI=paste0("(",round(REML_Model_All_Null_phylo_trad_r_PredictionIntervals[[3]],3),",",round(REML_Model_All_Null_phylo_trad_r_PredictionIntervals[[4]],3),")"),
                                          PI=paste0("(",round(REML_Model_All_Null_phylo_trad_r_PredictionIntervals[[5]],3),",",round(REML_Model_All_Null_phylo_trad_r_PredictionIntervals[[6]],3),")"),
                                          tval=round((coef(summary(REML_Model_All_Null_trad_r)))$tval,3),
                                          pval=round((coef(summary(REML_Model_All_Null_trad_r)))$pval,3),
                                          I2_ObsID=round(I2_Model_All_Null_trad_r_Taxon_Obs_ID,2),
                                          I2_StudyID=round(I2_Model_All_Null_trad_r_Taxon_Study_ID,2),
                                          I2_Nonphylo=round(I2_Model_All_Null_trad_r_Taxon_Taxon,2), 
                                          I2_Phylo="-",
                                          I2_Total=round(Model_All_Null_trad_r_I2_Total,2)
)
Model_All_Null_trad_r_Table


##  * Phylogenetic (correcting for phylogeny), full dataset ####

REML_Model_All_Null_phylo_r <- rma.mv(r, r_Var,
                                      mod =~1,
                                      random=list(~ 1 | Taxon_Phylo, # phylo effect
                                                  ~ 1 | Taxon, # non-phylo effect
                                                  ~ 1 | Study_ID,
                                                  ~ 1 | Obs_ID), 
                                      R = list(Taxon_Phylo = PhyloTree), #phylogenetic relatedness
                                      method="REML",
                                      test="t", # using t dist rather than z 
                                      data=Data)
REML_Model_All_Null_phylo_r

# extract model statistics for Table 1 
REML_Model_All_Null_phylo_PredictionIntervals <- predict(REML_Model_All_Null_phylo_r, digits=3)
W <- diag(1/REML_Model_All_Null_phylo_r$vi)
X <- model.matrix(REML_Model_All_Null_phylo_r)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
REML_Model_All_Null_phylo_I2 <- i2_ml(REML_Model_All_Null_phylo_r, method = "matrix")

REML_Model_All_Null_phylo_r_Table <- data.frame(model= "Phylogenetic (full dataset)",
                                              n=length(unique(Data$Study_ID)),
                                              k=nrow(Data),
                                              r=round(REML_Model_All_Null_phylo_PredictionIntervals[[1]],3),
                                              CI=paste0("(",round(REML_Model_All_Null_phylo_PredictionIntervals[[3]],3),",",round(REML_Model_All_Null_phylo_PredictionIntervals[[4]],3),")"),
                                              PI=paste0("(",round(REML_Model_All_Null_phylo_PredictionIntervals[[5]],3),",",round(REML_Model_All_Null_phylo_PredictionIntervals[[6]],3),")"),
                                              tval=round((coef(summary(REML_Model_All_Null_phylo_r)))$tval,3),
                                              pval=round((coef(summary(REML_Model_All_Null_phylo_r)))$pval,3),
                                              I2_ObsID=round(REML_Model_All_Null_phylo_I2[["I2_Obs_ID"]],2),
                                              I2_StudyID=round(REML_Model_All_Null_phylo_I2[["I2_Study_ID"]],2),
                                              I2_Nonphylo=round(REML_Model_All_Null_phylo_I2[["I2_Taxon"]],2), 
                                              I2_Phylo=round(REML_Model_All_Null_phylo_I2[["I2_Taxon_Phylo"]],2),
                                              I2_Total=round(REML_Model_All_Null_phylo_I2[["I2_Total"]],2)
)
REML_Model_All_Null_phylo_r_Table


##  ** Test for publication bias and time-lag bias  ####
Data$z_sei <- sqrt(Data$z_Var)
REML_Model_All_Null_phylo_PubBias_z <- rma.mv(z, z_Var,
                                              mod =~z_sei + Year + PhyloCorrection,
                                              random=list(~ 1 | Taxon_Phylo, # phylo effect
                                                          ~ 1 | Taxon, # non-phylo effect
                                                          ~ 1 | Study_ID, 
                                                          ~ 1 | Obs_ID), 
                                              R = list(Taxon_Phylo = PhyloTree), #phylogenetic relatedness
                                              method="REML",
                                              test="t", # using t dist rather than z 
                                              data=Data)
REML_Model_All_Null_phylo_PubBias_z


##  * Phylogenetic (correcting for phylogeny), filtered dataset ####

REML_Model_All_Null_phylo_r_cons <- rma.mv(r, r_Var,
                                           mod =~1,
                                           random=list(~ 1 | Taxon_Phylo, # phylo effect
                                                       ~ 1 | Taxon, # non-phylo effect
                                                       ~ 1 | Study_ID,
                                                       ~ 1 | Obs_ID), 
                                           R = list(Taxon_Phylo = PhyloTree), #phylogenetic relatedness
                                           method="REML",
                                           test="t", # using t dist rather than z 
                                           data=Data_conservative)
REML_Model_All_Null_phylo_r_cons

# extract model statistics for Table 1 
REML_Model_All_Null_phylo_r_cons_PredictionIntervals <- predict(REML_Model_All_Null_phylo_r_cons, digits=3)
REML_Model_All_Null_phylo_r_cons_I2 <- i2_ml(REML_Model_All_Null_phylo_r_cons, method = "matrix")
REML_Model_All_Null_phylo_r_cons_Table <- data.frame(model= "Phylogenetic (filtered dataset)",
                                                     n=length(unique(Data_conservative$Study_ID)),
                                                     k=nrow(Data_conservative),
                                                     r=round(REML_Model_All_Null_phylo_r_cons_PredictionIntervals[[1]],3),
                                                     CI=paste0("(",round(REML_Model_All_Null_phylo_r_cons_PredictionIntervals[[3]],3),",",round(REML_Model_All_Null_phylo_r_cons_PredictionIntervals[[4]],3),")"),
                                                     PI=paste0("(",round(REML_Model_All_Null_phylo_r_cons_PredictionIntervals[[5]],3),",",round(REML_Model_All_Null_phylo_r_cons_PredictionIntervals[[6]],3),")"),
                                                     tval=round((coef(summary(REML_Model_All_Null_phylo_r_cons)))$tval,3),
                                                     pval=round((coef(summary(REML_Model_All_Null_phylo_r_cons)))$pval,3),
                                                     I2_ObsID=round(REML_Model_All_Null_phylo_r_cons_I2[["I2_Obs_ID"]],2),
                                                     I2_StudyID=round(REML_Model_All_Null_phylo_r_cons_I2[["I2_Study_ID"]],2),
                                                     I2_Nonphylo=round(REML_Model_All_Null_phylo_r_cons_I2[["I2_Taxon"]],2), 
                                                     I2_Phylo=round(REML_Model_All_Null_phylo_r_cons_I2[["I2_Taxon_Phylo"]],2),
                                                     I2_Total=round(REML_Model_All_Null_phylo_r_cons_I2[["I2_Total"]],2)
)
REML_Model_All_Null_phylo_r_cons_Table


## ___ ####
# 4. Meta-analysis  — Null model z ####
##  * Traditional (not correcting for phylogeny), full dataset  ####

REML_Model_All_Null_trad_z <- rma.mv(z, z_Var,
                                     mod =~1,
                                     random=list(~ 1 | Taxon, # non-phylo effect
                                                 ~ 1 | Study_ID, 
                                                 ~ 1 | Obs_ID), 
                                     method="REML",
                                     test="t", # using t dist rather than z 
                                     data=Data)
REML_Model_All_Null_trad_z

# extract model statistics for Table S4
REML_Model_All_Null_phylo_trad_z_PredictionIntervals <- predict(REML_Model_All_Null_trad_z, digits=3)
W <- diag(1/REML_Model_All_Null_trad_z$vi)
X <- model.matrix(REML_Model_All_Null_trad_z)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
Model_All_Null_trad_z_I2_Total <- 100 * sum(REML_Model_All_Null_trad_z$sigma2) / (sum(REML_Model_All_Null_trad_z$sigma2) + (REML_Model_All_Null_trad_z$k-REML_Model_All_Null_trad_z$p)/sum(diag(P)))
Model_All_Null_trad_z_Sum_I2 <- sum(REML_Model_All_Null_trad_z$sigma2)
I2_Model_All_Null_trad_z_Taxon_Taxon    <- (REML_Model_All_Null_trad_z$sigma2[1]/Model_All_Null_trad_z_Sum_I2)*Model_All_Null_trad_z_I2_Total
I2_Model_All_Null_trad_z_Taxon_Study_ID <- (REML_Model_All_Null_trad_z$sigma2[2]/Model_All_Null_trad_z_Sum_I2)*Model_All_Null_trad_z_I2_Total
I2_Model_All_Null_trad_z_Taxon_Obs_ID   <- (REML_Model_All_Null_trad_z$sigma2[3]/Model_All_Null_trad_z_Sum_I2)*Model_All_Null_trad_z_I2_Total

Model_All_Null_trad_z_Table <- data.frame(model= "Traditional (full dataset)",
                                          n=length(unique(Data$Study_ID)),
                                          k=nrow(Data),
                                          z=round(REML_Model_All_Null_phylo_trad_z_PredictionIntervals[[1]],3),
                                          CI=paste0("(",round(REML_Model_All_Null_phylo_trad_z_PredictionIntervals[[3]],3),",",round(REML_Model_All_Null_phylo_trad_z_PredictionIntervals[[4]],3),")"),
                                          PI=paste0("(",round(REML_Model_All_Null_phylo_trad_z_PredictionIntervals[[5]],3),",",round(REML_Model_All_Null_phylo_trad_z_PredictionIntervals[[6]],3),")"),
                                          tval=round((coef(summary(REML_Model_All_Null_trad_z)))$tval,3),
                                          pval=round((coef(summary(REML_Model_All_Null_trad_z)))$pval,3),
                                          I2_ObsID=round(I2_Model_All_Null_trad_z_Taxon_Obs_ID,2),
                                          I2_StudyID=round(I2_Model_All_Null_trad_z_Taxon_Study_ID,2),
                                          I2_Nonphylo=round(I2_Model_All_Null_trad_z_Taxon_Taxon,2), 
                                          I2_Phylo="-",
                                          I2_Total=round(Model_All_Null_trad_z_I2_Total,2)
)
Model_All_Null_trad_z_Table


##  * Phylogenetic (correcting for phylogeny), full dataset ####

REML_Model_All_Null_phylo_z <- rma.mv(z, z_Var,
                                      mod =~1,
                                      random=list(~ 1 | Taxon_Phylo, # phylo effect
                                                  ~ 1 | Taxon, # non-phylo effect
                                                  ~ 1 | Study_ID,
                                                  ~ 1 | Obs_ID), 
                                      R = list(Taxon_Phylo = PhyloTree), #phylogenetic relatedness
                                      method="REML",
                                      test="t", # using t dist rather than z 
                                      data=Data)
REML_Model_All_Null_phylo_z

# extract model statistics for Table S4
REML_Model_All_Null_phylo_PredictionIntervals <- predict(REML_Model_All_Null_phylo_z, digits=3)
W <- diag(1/REML_Model_All_Null_phylo_z$vi)
X <- model.matrix(REML_Model_All_Null_phylo_z)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
REML_Model_All_Null_phylo_I2 <- i2_ml(REML_Model_All_Null_phylo_z, method = "matrix")

REML_Model_All_Null_phylo_z_Table <- data.frame(model= "Phylogenetic (full dataset)",
                                              n=length(unique(Data$Study_ID)),
                                              k=nrow(Data),
                                              z=round(REML_Model_All_Null_phylo_PredictionIntervals[[1]],3),
                                              CI=paste0("(",round(REML_Model_All_Null_phylo_PredictionIntervals[[3]],3),",",round(REML_Model_All_Null_phylo_PredictionIntervals[[4]],3),")"),
                                              PI=paste0("(",round(REML_Model_All_Null_phylo_PredictionIntervals[[5]],3),",",round(REML_Model_All_Null_phylo_PredictionIntervals[[6]],3),")"),
                                              tval=round((coef(summary(REML_Model_All_Null_phylo_z)))$tval,3),
                                              pval=round((coef(summary(REML_Model_All_Null_phylo_z)))$pval,3),
                                              I2_ObsID=round(REML_Model_All_Null_phylo_I2[["I2_Obs_ID"]],2),
                                              I2_StudyID=round(REML_Model_All_Null_phylo_I2[["I2_Study_ID"]],2),
                                              I2_Nonphylo=round(REML_Model_All_Null_phylo_I2[["I2_Taxon"]],2), 
                                              I2_Phylo=round(REML_Model_All_Null_phylo_I2[["I2_Taxon_Phylo"]],2),
                                              I2_Total=round(REML_Model_All_Null_phylo_I2[["I2_Total"]],2)
)
REML_Model_All_Null_phylo_z_Table


##  * Phylogenetic (correcting for phylogeny), filtered dataset ####

REML_Model_All_Null_phylo_z_cons <- rma.mv(z, z_Var,
                                           mod =~1,
                                           random=list(~ 1 | Taxon_Phylo, # phylo effect
                                                       ~ 1 | Taxon, # non-phylo effect
                                                       ~ 1 | Study_ID,
                                                       ~ 1 | Obs_ID), 
                                           R = list(Taxon_Phylo = PhyloTree), #phylogenetic relatedness
                                           method="REML",
                                           test="t", # using t dist rather than z 
                                           data=Data_conservative)
REML_Model_All_Null_phylo_z_cons

# extract model statistics for Table S4
REML_Model_All_Null_phylo_z_cons_PredictionIntervals <- predict(REML_Model_All_Null_phylo_z_cons, digits=3)
REML_Model_All_Null_phylo_z_cons_I2 <- i2_ml(REML_Model_All_Null_phylo_z_cons, method = "matrix")
REML_Model_All_Null_phylo_z_cons_Table <- data.frame(model= "Phylogenetic (filtered dataset)",
                                                     n=length(unique(Data_conservative$Study_ID)),
                                                     k=nrow(Data_conservative),
                                                     z=round(REML_Model_All_Null_phylo_z_cons_PredictionIntervals[[1]],3),
                                                     CI=paste0("(",round(REML_Model_All_Null_phylo_z_cons_PredictionIntervals[[3]],3),",",round(REML_Model_All_Null_phylo_z_cons_PredictionIntervals[[4]],3),")"),
                                                     PI=paste0("(",round(REML_Model_All_Null_phylo_z_cons_PredictionIntervals[[5]],3),",",round(REML_Model_All_Null_phylo_z_cons_PredictionIntervals[[6]],3),")"),
                                                     tval=round((coef(summary(REML_Model_All_Null_phylo_z_cons)))$tval,3),
                                                     pval=round((coef(summary(REML_Model_All_Null_phylo_z_cons)))$pval,3),
                                                     I2_ObsID=round(REML_Model_All_Null_phylo_z_cons_I2[["I2_Obs_ID"]],2),
                                                     I2_StudyID=round(REML_Model_All_Null_phylo_z_cons_I2[["I2_Study_ID"]],2),
                                                     I2_Nonphylo=round(REML_Model_All_Null_phylo_z_cons_I2[["I2_Taxon"]],2), 
                                                     I2_Phylo=round(REML_Model_All_Null_phylo_z_cons_I2[["I2_Taxon_Phylo"]],2),
                                                     I2_Total=round(REML_Model_All_Null_phylo_z_cons_I2[["I2_Total"]],2)
)
REML_Model_All_Null_phylo_z_cons_Table


## ___ ####
# 5. Meta-analysis — Moderator models  r####

## * Taxonomic clade ####
Data_TaxonClass_Model  <-  subset(Data, Class  != 'Amphibia') # sample size k < 5
Data_TaxonClass_Model  <-  subset(Data_TaxonClass_Model, Class  != 'Arachnida') # sample size k < 5

REML_r_Mod_Class   = rma.mv(r ~ Class,     V = r_Var,   data = Data_TaxonClass_Model, random = c(~ 1 | Obs_ID, ~ 1 | Study_ID), method = "REML")
summary(REML_r_Mod_Class)

# extract model statistics for Table 2 
REML_r_Mod_Class_Name <- "Taxonomic clade"
REML_r_Mod_Class_k <- REML_r_Mod_Class$k
REML_r_Mod_Class_Qdf <- REML_r_Mod_Class$QMdf[1]
REML_r_Mod_Class_QM <- REML_r_Mod_Class$QM
REML_r_Mod_Class_Qpval <- REML_r_Mod_Class$QMp

Model_Results_REML_r_Mod_Class <- data.frame(Model    = REML_r_Mod_Class_Name,
                                             k        = REML_r_Mod_Class_k,
                                             num_dfs  = REML_r_Mod_Class_Qdf,
                                             QM       = round(REML_r_Mod_Class_QM,2),
                                             Pval     = round(REML_r_Mod_Class_Qpval,3))
Model_Results_REML_r_Mod_Class


##  * Speciation proxy ####

REML_r_Mod_Speciation_proxy   = rma.mv(r ~ Speciation_proxy,     V = r_Var,   data = Data, random = c(~ 1 | Obs_ID, ~ 1 | Study_ID,  ~ 1 | Taxon,~ 1 | Taxon_Phylo), R = list(Taxon_Phylo = PhyloTree), method = "REML")
summary(REML_r_Mod_Speciation_proxy)

# extract model statistics for Table 2 
REML_r_Mod_Speciation_proxy_Name <- "Speciation proxy"
REML_r_Mod_Speciation_proxy_k <- REML_r_Mod_Speciation_proxy$k
REML_r_Mod_Speciation_proxy_Qdf <- REML_r_Mod_Speciation_proxy$QMdf[1]
REML_r_Mod_Speciation_proxy_QM <- REML_r_Mod_Speciation_proxy$QM
REML_r_Mod_Speciation_proxy_Qpval <- REML_r_Mod_Speciation_proxy$QMp
Model_Results_REML_r_Mod_Speciation_proxy <- data.frame(Model    = REML_r_Mod_Speciation_proxy_Name,
                                                        k        = REML_r_Mod_Speciation_proxy_k,
                                                        num_dfs  = REML_r_Mod_Speciation_proxy_Qdf,
                                                        QM       = round(REML_r_Mod_Speciation_proxy_QM,2),
                                                        Pval     = round(REML_r_Mod_Speciation_proxy_Qpval,3))
Model_Results_REML_r_Mod_Speciation_proxy

Speciation_proxy_order <- c("Species richness",   "Diversification rate", "Speciation rate")
palette_Speciation_proxy=c("#A1DAB4", "#41B6C4","#253494")

OrchaRD_Plot_Speciation_proxy <- orchaRd::orchard_plot(REML_r_Mod_Speciation_proxy, 
                                                        mod = "Speciation_proxy", group = "Study_ID", 
                                                        xlab = expression(paste("Effect size (", italic(r),')')), angle = 0, alpha = 0.6,
                                                        g = TRUE,
                                                        k = TRUE,
                                                        transfm = "none", 
                                                        twig.size = 0.75, trunk.size = 0.7,  branch.size = 1.5,
                                                        tree.order = c("Species richness",   "Diversification rate", "Speciation rate")) +
  labs(title = "", x = "Speciation proxy", y = expression(paste("Effect size (", italic(r),')'))) +
  theme(axis.line = element_line(lineend = "square")) +
  theme(legend.position.inside = c(-0.02, 1.13), 
        legend.justification = c(0, 1), 
        legend.key.size = unit(1, "mm"),
        legend.direction = "horizontal", 
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 10),
        legend.key = element_blank()) +
  scale_fill_manual(values = palette_Speciation_proxy) + 
  scale_colour_manual(values = palette_Speciation_proxy) +
  scale_x_discrete(limits = Speciation_proxy_order) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "none",
        axis.line.x = element_line(colour = "black", linewidth = 0.7, lineend = "square"),
        axis.line.y = element_line(colour = "black", linewidth = 0.7, lineend = "square"),
        axis.text.x = element_text(face="plain", color="black", size=12, angle=0),
        axis.text.y = element_text(face="plain", color="black", size=12, angle=0),
        axis.title.x = element_text(size=14,face="plain", margin = margin(r=10,0,0,0), vjust = -1),
        axis.title.y = element_text(size=14,face="plain", margin = margin(r=10,0,0,0)),
        axis.ticks = element_line(linewidth = 0.7, colour = "black"),
        axis.ticks.length = unit(.3, "cm"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill= "transparent"),
        plot.margin = margin(t = 1,  # Top margin
                             r =  1,  # Right margin
                             b =  1,  # Bottom margin
                             l =  1)) # Left margin 
OrchaRD_Plot_Speciation_proxy


##  * Sexual selection proxy ####

REML_r_Mod_SexSel_proxy   = rma.mv(r ~ SexSel_proxy,     V = r_Var,   data = Data, random = c(~ 1 | Obs_ID, ~ 1 | Study_ID,  ~ 1 | Taxon,~ 1 | Taxon_Phylo), R = list(Taxon_Phylo = PhyloTree), method = "REML")
summary(REML_r_Mod_SexSel_proxy)

# extract model statistics for Table 2
REML_r_Mod_SexSel_proxy_Name <- "Sexual selection proxy"
REML_r_Mod_SexSel_proxy_k <- REML_r_Mod_SexSel_proxy$k
REML_r_Mod_SexSel_proxy_Qdf <- REML_r_Mod_SexSel_proxy$QMdf[1]
REML_r_Mod_SexSel_proxy_QM <- REML_r_Mod_SexSel_proxy$QM
REML_r_Mod_SexSel_proxy_Qpval <- REML_r_Mod_SexSel_proxy$QMp
Model_Results_REML_r_Mod_SexSel_proxy <- data.frame(Model    = REML_r_Mod_SexSel_proxy_Name,
                                                    k        = REML_r_Mod_SexSel_proxy_k,
                                                    num_dfs  = REML_r_Mod_SexSel_proxy_Qdf,
                                                    QM       = round(REML_r_Mod_SexSel_proxy_QM,2),
                                                    Pval     = round(REML_r_Mod_SexSel_proxy_Qpval,3))
Model_Results_REML_r_Mod_SexSel_proxy

palette_SexSelProxy=c("darkgrey", "#A1DAB4", "#41B6C4", "#2C7FB8", "#253494")
SexSelProxy_order <- c("Other",   "Trait", "Mating system", "Dichromatism",  "Sexual size dimorphism")

OrchaRD_Plot_SexSel_proxy <- orchaRd::orchard_plot(REML_r_Mod_SexSel_proxy, 
                                                    mod = "SexSel_proxy", group = "Study_ID", 
                                                    xlab = expression(paste("Effect size (", italic(r),')')), angle = 0, alpha = 0.6,
                                                    g = TRUE,
                                                    k = TRUE,
                                                    transfm = "none", 
                                                    twig.size = 0.75, trunk.size = 0.7,  branch.size = 1.5,
                                                    tree.order = c("Other",   "Trait", "Mating system", "Dichromatism",  "Sexual size dimorphism")) +
  labs(title = "", x = "Sexual selection proxy", y = expression(paste("Effect size (", italic(r),')'))) +
  theme(axis.line = element_line(lineend = "square")) +
  theme(legend.position.inside = c(0, 1.09), 
        legend.justification = c(0, 1), 
        legend.key.size = unit(1, "mm"),
        legend.direction = "horizontal", 
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 12),
        legend.key = element_blank()) +
  scale_fill_manual(values = palette_SexSelProxy) + 
  scale_colour_manual(values = palette_SexSelProxy) +
  scale_x_discrete(limits = SexSelProxy_order) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line.x = element_line(colour = "black", linewidth = 0.7, lineend = "square"),
        axis.line.y = element_line(colour = "black", linewidth = 0.7, lineend = "square"),
        axis.text.x = element_text(face="plain", color="black", size=12, angle=0),
        axis.text.y = element_text(face="plain", color="black", size=12, angle=0),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=14,face="plain", margin = margin(r=10,0,0,0)),
        axis.ticks = element_line(linewidth = 0.7, colour = "black"),
        axis.ticks.length = unit(.3, "cm"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill= "transparent"),
        plot.margin = margin(t =  1,  # Top margin
                             r =  1,  # Right margin
                             b =  1,  # Bottom margin
                             l =  1)) # Left margin 
OrchaRD_Plot_SexSel_proxy


##  * Sex-specific sexual selection ####

REML_r_Mod_SexSel_proxy_sex   = rma.mv(r ~ SexSel_proxy_sex,     V = r_Var,   data = Data, random = c(~ 1 | Obs_ID, ~ 1 | Study_ID,  ~ 1 | Taxon,~ 1 | Taxon_Phylo), R = list(Taxon_Phylo = PhyloTree), method = "REML")
summary(REML_r_Mod_SexSel_proxy_sex)

# extract model statistics for Table 2
REML_r_Mod_SexSel_proxy_sex_Name <- "Sex-specific sexual selection"
REML_r_Mod_SexSel_proxy_sex_k <- REML_r_Mod_SexSel_proxy_sex$k
REML_r_Mod_SexSel_proxy_sex_Qdf <- REML_r_Mod_SexSel_proxy_sex$QMdf[1]
REML_r_Mod_SexSel_proxy_sex_QM <- REML_r_Mod_SexSel_proxy_sex$QM
REML_r_Mod_SexSel_proxy_sex_Qpval <- REML_r_Mod_SexSel_proxy_sex$QMp

Model_Results_REML_r_Mod_SexSel_proxy_sex <- data.frame(Model    = REML_r_Mod_SexSel_proxy_sex_Name,
                                                        k        = REML_r_Mod_SexSel_proxy_sex_k,
                                                        num_dfs  = REML_r_Mod_SexSel_proxy_sex_Qdf,
                                                        QM       = round(REML_r_Mod_SexSel_proxy_sex_QM,2),
                                                        Pval     = round(REML_r_Mod_SexSel_proxy_sex_Qpval,3))
Model_Results_REML_r_Mod_SexSel_proxy_sex


palette_SexSel_proxy_sex=c("#A1DAB4", "#41B6C4","#253494")
SexSel_proxy_sex_order <- c("Both",   "Female", "Male")

OrchaRD_Plot_SexSel_proxy_sex <- orchaRd::orchard_plot(REML_r_Mod_SexSel_proxy_sex, 
                                                       mod = "SexSel_proxy_sex", group = "Study_ID", 
                                                       xlab = expression(paste("Effect size (", italic(r),')')), angle = 0, alpha = 0.6,
                                                       g = TRUE,
                                                       k = TRUE,
                                                       transfm = "none", 
                                                       twig.size = 0.75, trunk.size = 0.7,  branch.size = 1.5,
                                                       tree.order = c("Both",   "Female", "Male")) +
  labs(title = "", x = "Sex", y = expression(paste("Effect size (", italic(r),')'))) +
  theme(axis.line = element_line(lineend = "square")) +
  theme(legend.position.inside = c(0.00, 1.13), 
        #legend.position.inside = c(-0.02, 1.13), 
        legend.justification = c(0, 1), 
        legend.key.size = unit(1, "mm"),
        legend.direction = "horizontal", 
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 12),
        legend.key = element_blank()) +
  scale_fill_manual(values = palette_SexSel_proxy_sex) + 
  scale_colour_manual(values = palette_SexSel_proxy_sex) +
  scale_x_discrete(limits = SexSel_proxy_sex_order) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "none",
        axis.line.x = element_line(colour = "black", linewidth = 0.7, lineend = "square"),
        axis.line.y = element_line(colour = "black", linewidth = 0.7, lineend = "square"),
        axis.text.x = element_text(face="plain", color="black", size=12, angle=0),
        axis.text.y = element_text(face="plain", color="black", size=12, angle=0),
        axis.title.x = element_text(size=14,face="plain", margin = margin(r=10,0,0,0), vjust = -1),
        axis.title.y = element_text(size=14,face="plain", margin = margin(r=10,0,0,0)),
        axis.ticks = element_line(linewidth = 0.7, colour = "black"),
        axis.ticks.length = unit(.3, "cm"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill= "transparent"),
        plot.margin = margin(t =  1,  # Top margin
                             r =  1,  # Right margin
                             b =  1,  # Bottom margin
                             l =  1)) # Left margin 
OrchaRD_Plot_SexSel_proxy_sex


##  * Sexual selection mechanism ####

REML_r_Mod_SexSel_proxy_mechanism   = rma.mv(r ~ SexSel_proxy_mechanism,     V = r_Var,   data = Data, random = c(~ 1 | Obs_ID, ~ 1 | Study_ID,  ~ 1 | Taxon,~ 1 | Taxon_Phylo), R = list(Taxon_Phylo = PhyloTree), method = "REML")
summary(REML_r_Mod_SexSel_proxy_mechanism)

# extract model statistics for Table 2
REML_r_Mod_SexSel_proxy_mechanism_Name  <- "Sexual selection mechanism"
REML_r_Mod_SexSel_proxy_mechanism_k     <- REML_r_Mod_SexSel_proxy_mechanism$k
REML_r_Mod_SexSel_proxy_mechanism_Qdf   <- REML_r_Mod_SexSel_proxy_mechanism$QMdf[1]
REML_r_Mod_SexSel_proxy_mechanism_QM    <- REML_r_Mod_SexSel_proxy_mechanism$QM
REML_r_Mod_SexSel_proxy_mechanism_Qpval <- REML_r_Mod_SexSel_proxy_mechanism$QMp

Model_Results_REML_r_Mod_SexSel_proxy_mechanism <- data.frame(Model    = REML_r_Mod_SexSel_proxy_mechanism_Name,
                                                              k         = REML_r_Mod_SexSel_proxy_mechanism_k,
                                                              num_dfs   = REML_r_Mod_SexSel_proxy_mechanism_Qdf,
                                                              QM        = round(REML_r_Mod_SexSel_proxy_mechanism_QM,2),
                                                              Pval      = round(REML_r_Mod_SexSel_proxy_mechanism_Qpval,3))
Model_Results_REML_r_Mod_SexSel_proxy_mechanism

palette_SexSel_proxy_mechanism=c("#A1DAB4", "#41B6C4","#253494")
SexSel_proxy_mechanism_order <- c("Both",   "Competition", "Choice")

OrchaRD_Plot_SexSel_proxy_mechanism <- orchaRd::orchard_plot(REML_r_Mod_SexSel_proxy_mechanism, 
                                                             mod = "SexSel_proxy_mechanism", group = "Study_ID", 
                                                             xlab = expression(paste("Effect size (", italic(r),')')), angle = 0, alpha = 0.6,
                                                             g = TRUE,
                                                             k = TRUE,
                                                             transfm = "none", 
                                                             twig.size = 0.75, trunk.size = 0.7,  branch.size = 1.5,
                                                             tree.order = c("Both", "Competition", "Choice")) +
  labs(title = "", x = "Sexual selection mechanism", y = expression(paste("Effect size (", italic(r),')'))) +
  theme(axis.line = element_line(lineend = "square")) +
  theme(legend.position.inside = c(0.00, 1.09), 
        legend.justification = c(0, 1), 
        legend.key.size = unit(1, "mm"),
        legend.direction = "horizontal", 
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 12),
        legend.key = element_blank()) +
  scale_fill_manual(values = palette_SexSel_proxy_mechanism) + 
  scale_colour_manual(values = palette_SexSel_proxy_mechanism) +
  scale_x_discrete(limits = SexSel_proxy_mechanism_order) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line.x = element_line(colour = "black", linewidth = 0.7, lineend = "square"),
        axis.line.y = element_line(colour = "black", linewidth = 0.7, lineend = "square"),
        axis.text.x = element_text(face="plain", color="black", size=12, angle=0),
        axis.text.y = element_text(face="plain", color="black", size=12, angle=0),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=14,face="plain", margin = margin(r=10,0,0,0)),
        axis.ticks = element_line(linewidth = 0.7, colour = "black"),
        axis.ticks.length = unit(.3, "cm"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill= "transparent"),
        plot.margin = margin(t =  1,  # Top margin
                             r =  1,  # Right margin
                             b =  1,  # Bottom margin
                             l =  1)) # Left margin 
OrchaRD_Plot_SexSel_proxy_mechanism


##  * Mating stage ####
REML_r_Mod_SexSel_proxy_stage   = rma.mv(r ~ SexSel_proxy_stage,     V = r_Var,   data = Data, random = c(~ 1 | Obs_ID, ~ 1 | Study_ID,  ~ 1 | Taxon,~ 1 | Taxon_Phylo), R = list(Taxon_Phylo = PhyloTree), method = "REML")
summary(REML_r_Mod_SexSel_proxy_stage)

# extract model statistics for Table 2
REML_r_Mod_SexSel_proxy_stage_Name  <- "Mating stage"
REML_r_Mod_SexSel_proxy_stage_k     <- REML_r_Mod_SexSel_proxy_stage$k
REML_r_Mod_SexSel_proxy_stage_Qdf   <- REML_r_Mod_SexSel_proxy_stage$QMdf[1]
REML_r_Mod_SexSel_proxy_stage_QM    <- REML_r_Mod_SexSel_proxy_stage$QM
REML_r_Mod_SexSel_proxy_stage_Qpval <- REML_r_Mod_SexSel_proxy_stage$QMp

Model_Results_REML_r_Mod_SexSel_proxy_stage <- data.frame(Model    = REML_r_Mod_SexSel_proxy_stage_Name,
                                                          k         = REML_r_Mod_SexSel_proxy_stage_k,
                                                          num_dfs   = REML_r_Mod_SexSel_proxy_mechanism_Qdf,
                                                          QM        = round(REML_r_Mod_SexSel_proxy_stage_QM,2),
                                                          Pval      = round(REML_r_Mod_SexSel_proxy_stage_Qpval,3))
Model_Results_REML_r_Mod_SexSel_proxy_stage

palette_SexSel_proxy_stage=c("#A1DAB4", "#41B6C4","#253494")
SexSel_proxy_stage_order <- c("Both",   "Post-mating", "Pre-mating")

OrchaRD_Plot_SexSel_proxy_stage <- orchaRd::orchard_plot(REML_r_Mod_SexSel_proxy_stage, 
                                                         mod = "SexSel_proxy_stage", group = "Study_ID", 
                                                         xlab = expression(paste("Effect size (", italic(r),')')), angle = 0, alpha = 0.6,
                                                         g = TRUE,
                                                         k = TRUE,
                                                         transfm = "none", 
                                                         twig.size = 0.75, trunk.size = 0.7,  branch.size = 1.5,
                                                         tree.order = c("Both",   "Post-mating", "Pre-mating")) +
  labs(title = "", x = "Mating stage", y = expression(paste("Effect size (", italic(r),')'))) +
  theme(axis.line = element_line(lineend = "square")) +
  theme(legend.position.inside = c(0.00, 1.13), 
        #legend.position.inside = c(-0.02, 1.13), 
        legend.justification = c(0, 1), 
        legend.key.size = unit(1, "mm"),
        legend.direction = "horizontal", 
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 12),
        legend.key = element_blank()) +
  scale_fill_manual(values = palette_SexSel_proxy_stage) + 
  scale_colour_manual(values = palette_SexSel_proxy_stage) +
  scale_x_discrete(limits = SexSel_proxy_stage_order) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "none",
        axis.line.x = element_line(colour = "black", linewidth = 0.7, lineend = "square"),
        axis.line.y = element_line(colour = "black", linewidth = 0.7, lineend = "square"),
        axis.text.x = element_text(face="plain", color="black", size=12, angle=0),
        axis.text.y = element_text(face="plain", color="black", size=12, angle=0),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=14,face="plain", margin = margin(r=10,0,0,0)),
        axis.ticks = element_line(linewidth = 0.7, colour = "black"),
        axis.ticks.length = unit(.3, "cm"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill= "transparent"),
        plot.margin = margin(t =  1,  # Top margin
                             r =  1,  # Right margin
                             b =  1,  # Bottom margin
                             l =  1)) # Left margin 
OrchaRD_Plot_SexSel_proxy_stage


##  * Phylogenetic correction ####

REML_r_Mod_PhyloCorrection   = rma.mv(r ~ PhyloCorrection,     V = r_Var,   data = Data, random = c(~ 1 | Obs_ID, ~ 1 | Study_ID, ~ 1 | Taxon, ~ 1 | Taxon_Phylo), R = list(Taxon_Phylo = PhyloTree), method = "REML")
summary(REML_r_Mod_PhyloCorrection)

# extract model statistics for Table 2
REML_r_Mod_PhyloCorrection_Name <- "Phylogenetic correction"
REML_r_Mod_PhyloCorrection_k <- REML_r_Mod_PhyloCorrection$k
REML_r_Mod_PhyloCorrection_Qdf <- REML_r_Mod_PhyloCorrection$QMdf[1]
REML_r_Mod_PhyloCorrection_QM <- REML_r_Mod_PhyloCorrection$QM
REML_r_Mod_PhyloCorrection_Qpval <- REML_r_Mod_PhyloCorrection$QMp

Model_Results_REML_r_Mod_PhyloCorrection <- data.frame(Model    = REML_r_Mod_PhyloCorrection_Name,
                                                       k        = REML_r_Mod_PhyloCorrection_k,
                                                       num_dfs  = REML_r_Mod_PhyloCorrection_Qdf,
                                                       QM       = round(REML_r_Mod_PhyloCorrection_QM,2),
                                                       Pval     = round(REML_r_Mod_PhyloCorrection_Qpval,3))
Model_Results_REML_r_Mod_PhyloCorrection

palette_PhyloCorrection=c("#41B6C4", "#253494")
PhyloCorrection_order <- c("No",   "Yes")

OrchaRD_Plot_PhyloCorrection <- orchaRd::orchard_plot(REML_r_Mod_PhyloCorrection, 
                                                      mod = "PhyloCorrection", group = "Study_ID", 
                                                      xlab = expression(paste("Effect size (", italic(r),')')), angle = 0, alpha = 0.6,
                                                      g = TRUE,
                                                      k = TRUE,
                                                      transfm = "none", 
                                                      twig.size = 0.75, trunk.size = 0.7,  branch.size = 1.5,
                                                      tree.order = c("No",   "Yes")) +
  labs(title = "", x = "Phylogenetic correction", y = expression(paste("Effect size (", italic(r),')'))) +
  theme(axis.line = element_line(lineend = "square")) +
  theme(legend.position.inside = c(-0.0005, 1.15), 
        legend.justification = c(0, 1), 
        legend.key.size = unit(1, "mm"),
        legend.direction = "horizontal", 
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 12),
        legend.key = element_blank()) +
  scale_fill_manual(values = palette_PhyloCorrection) + 
  scale_colour_manual(values = palette_PhyloCorrection) +
  scale_x_discrete(limits = PhyloCorrection_order) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "none",
        axis.line.x = element_line(colour = "black", linewidth = 0.7, lineend = "square"),
        axis.line.y = element_line(colour = "black", linewidth = 0.7, lineend = "square"),
        axis.text.x = element_text(face="plain", color="black", size=12, angle=0),
        axis.text.y = element_text(face="plain", color="black", size=12, angle=0),
        axis.title.x = element_text(size=14,face="plain", margin = margin(r=10,0,0,0), vjust = -1),
        axis.title.y = element_text(size=14,face="plain", margin = margin(r=10,0,0,0)),
        axis.ticks = element_line(linewidth = 0.7, colour = "black"),
        axis.ticks.length = unit(.3, "cm"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill= "transparent"),
        plot.margin = margin(t =  1,  # Top margin
                             r =  1,  # Right margin
                             b =  1,  # Bottom margin
                             l =  1)) # Left margin 
OrchaRD_Plot_PhyloCorrection


## ___ ####
# 6. Meta-analysis  — Moderator models  z ####

##  * Taxonomic clade ####
REML_z_Mod_Class   = rma.mv(z ~ Class,     V = z_Var,   data = Data_TaxonClass_Model, random = c(~ 1 | Obs_ID, ~ 1 | Study_ID), method = "REML")
summary(REML_z_Mod_Class)

# extract model statistics for Table S5 
REML_z_Mod_Class_Name <- "Taxonomic clade"
REML_z_Mod_Class_k <- REML_z_Mod_Class$k
REML_z_Mod_Class_Qdf <- REML_z_Mod_Class$QMdf[1]
REML_z_Mod_Class_QM <- REML_z_Mod_Class$QM
REML_z_Mod_Class_Qpval <- REML_z_Mod_Class$QMp

Model_Results_REML_z_Mod_Class <- data.frame(Model    = REML_z_Mod_Class_Name,
                                             k        = REML_z_Mod_Class_k,
                                             num_dfs  = REML_z_Mod_Class_Qdf,
                                             QM       = round(REML_z_Mod_Class_QM,2),
                                             Pval     = round(REML_z_Mod_Class_Qpval,3))
Model_Results_REML_z_Mod_Class


##  * Speciation proxy ####

REML_z_Mod_Speciation_proxy   = rma.mv(z ~ Speciation_proxy,     V = z_Var,   data = Data, random = c(~ 1 | Obs_ID, ~ 1 | Study_ID,  ~ 1 | Taxon,~ 1 | Taxon_Phylo), R = list(Taxon_Phylo = PhyloTree), method = "REML")
summary(REML_z_Mod_Speciation_proxy)

# extract model statistics for Table S5
REML_z_Mod_Speciation_proxy_Name <- "Speciation proxy"
REML_z_Mod_Speciation_proxy_k <- REML_z_Mod_Speciation_proxy$k
REML_z_Mod_Speciation_proxy_Qdf <- REML_z_Mod_Speciation_proxy$QMdf[1]
REML_z_Mod_Speciation_proxy_QM <- REML_z_Mod_Speciation_proxy$QM
REML_z_Mod_Speciation_proxy_Qpval <- REML_z_Mod_Speciation_proxy$QMp
Model_Results_REML_z_Mod_Speciation_proxy <- data.frame(Model    = REML_z_Mod_Speciation_proxy_Name,
                                                        k        = REML_z_Mod_Speciation_proxy_k,
                                                        num_dfs  = REML_z_Mod_Speciation_proxy_Qdf,
                                                        QM       = round(REML_z_Mod_Speciation_proxy_QM,2),
                                                        Pval     = round(REML_z_Mod_Speciation_proxy_Qpval,3))
Model_Results_REML_z_Mod_Speciation_proxy


## * Sexual selection proxy ####

REML_z_Mod_SexSel_proxy   = rma.mv(z ~ SexSel_proxy,     V = z_Var,   data = Data, random = c(~ 1 | Obs_ID, ~ 1 | Study_ID,  ~ 1 | Taxon,~ 1 | Taxon_Phylo), R = list(Taxon_Phylo = PhyloTree), method = "REML")
summary(REML_z_Mod_SexSel_proxy)

# extract model statistics for Table S5
REML_z_Mod_SexSel_proxy_Name <- "Sexual selection proxy"
REML_z_Mod_SexSel_proxy_k <- REML_z_Mod_SexSel_proxy$k
REML_z_Mod_SexSel_proxy_Qdf <- REML_z_Mod_SexSel_proxy$QMdf[1]
REML_z_Mod_SexSel_proxy_QM <- REML_z_Mod_SexSel_proxy$QM
REML_z_Mod_SexSel_proxy_Qpval <- REML_z_Mod_SexSel_proxy$QMp
Model_Results_REML_z_Mod_SexSel_proxy <- data.frame(Model    = REML_z_Mod_SexSel_proxy_Name,
                                                    k        = REML_z_Mod_SexSel_proxy_k,
                                                    num_dfs  = REML_z_Mod_SexSel_proxy_Qdf,
                                                    QM       = round(REML_z_Mod_SexSel_proxy_QM,2),
                                                    Pval     = round(REML_z_Mod_SexSel_proxy_Qpval,3))
Model_Results_REML_z_Mod_SexSel_proxy


##  * Sex-specific sexual selection ####

REML_z_Mod_SexSel_proxy_sex   = rma.mv(z ~ SexSel_proxy_sex,     V = z_Var,   data = Data, random = c(~ 1 | Obs_ID, ~ 1 | Study_ID,  ~ 1 | Taxon,~ 1 | Taxon_Phylo), R = list(Taxon_Phylo = PhyloTree), method = "REML")
summary(REML_z_Mod_SexSel_proxy_sex)

# extract model statistics for Table S5
REML_z_Mod_SexSel_proxy_sex_Name <- "Sex-specific sexual selection"
REML_z_Mod_SexSel_proxy_sex_k <- REML_z_Mod_SexSel_proxy_sex$k
REML_z_Mod_SexSel_proxy_sex_Qdf <- REML_z_Mod_SexSel_proxy_sex$QMdf[1]
REML_z_Mod_SexSel_proxy_sex_QM <- REML_z_Mod_SexSel_proxy_sex$QM
REML_z_Mod_SexSel_proxy_sex_Qpval <- REML_z_Mod_SexSel_proxy_sex$QMp

Model_Results_REML_z_Mod_SexSel_proxy_sex <- data.frame(Model    = REML_z_Mod_SexSel_proxy_sex_Name,
                                                        k        = REML_z_Mod_SexSel_proxy_sex_k,
                                                        num_dfs  = REML_z_Mod_SexSel_proxy_sex_Qdf,
                                                        QM       = round(REML_z_Mod_SexSel_proxy_sex_QM,2),
                                                        Pval     = round(REML_z_Mod_SexSel_proxy_sex_Qpval,3))
Model_Results_REML_z_Mod_SexSel_proxy_sex


##  * Sexual selection mechanism ####

REML_z_Mod_SexSel_proxy_mechanism   = rma.mv(z ~ SexSel_proxy_mechanism,     V = z_Var,   data = Data, random = c(~ 1 | Obs_ID, ~ 1 | Study_ID,  ~ 1 | Taxon,~ 1 | Taxon_Phylo), R = list(Taxon_Phylo = PhyloTree), method = "REML")
summary(REML_z_Mod_SexSel_proxy_mechanism)

# extract model statistics for Table S5
REML_z_Mod_SexSel_proxy_mechanism_Name  <- "Sexual selection mechanism"
REML_z_Mod_SexSel_proxy_mechanism_k     <- REML_z_Mod_SexSel_proxy_mechanism$k
REML_z_Mod_SexSel_proxy_mechanism_Qdf   <- REML_z_Mod_SexSel_proxy_mechanism$QMdf[1]
REML_z_Mod_SexSel_proxy_mechanism_QM    <- REML_z_Mod_SexSel_proxy_mechanism$QM
REML_z_Mod_SexSel_proxy_mechanism_Qpval <- REML_z_Mod_SexSel_proxy_mechanism$QMp

Model_Results_REML_z_Mod_SexSel_proxy_mechanism <- data.frame(Model    = REML_z_Mod_SexSel_proxy_mechanism_Name,
                                                              k         = REML_z_Mod_SexSel_proxy_mechanism_k,
                                                              num_dfs   = REML_z_Mod_SexSel_proxy_mechanism_Qdf,
                                                              QM        = round(REML_z_Mod_SexSel_proxy_mechanism_QM,2),
                                                              Pval      = round(REML_z_Mod_SexSel_proxy_mechanism_Qpval,3))
Model_Results_REML_z_Mod_SexSel_proxy_mechanism


##  * Mating stage ####

REML_z_Mod_SexSel_proxy_stage   = rma.mv(z ~ SexSel_proxy_stage,     V = z_Var,   data = Data, random = c(~ 1 | Obs_ID, ~ 1 | Study_ID,  ~ 1 | Taxon,~ 1 | Taxon_Phylo), R = list(Taxon_Phylo = PhyloTree), method = "REML")
summary(REML_z_Mod_SexSel_proxy_stage)

# extract model statistics for Table S5
REML_z_Mod_SexSel_proxy_stage_Name  <- "Mating stage"
REML_z_Mod_SexSel_proxy_stage_k     <- REML_z_Mod_SexSel_proxy_stage$k
REML_z_Mod_SexSel_proxy_stage_Qdf   <- REML_z_Mod_SexSel_proxy_stage$QMdf[1]
REML_z_Mod_SexSel_proxy_stage_QM    <- REML_z_Mod_SexSel_proxy_stage$QM
REML_z_Mod_SexSel_proxy_stage_Qpval <- REML_z_Mod_SexSel_proxy_stage$QMp

Model_Results_REML_z_Mod_SexSel_proxy_stage <- data.frame(Model    = REML_z_Mod_SexSel_proxy_stage_Name,
                                                          k         = REML_z_Mod_SexSel_proxy_stage_k,
                                                          num_dfs   = REML_z_Mod_SexSel_proxy_mechanism_Qdf,
                                                          QM        = round(REML_z_Mod_SexSel_proxy_stage_QM,2),
                                                          Pval      = round(REML_z_Mod_SexSel_proxy_stage_Qpval,3))
Model_Results_REML_z_Mod_SexSel_proxy_stage


## * Phylogenetic correction ####

REML_z_Mod_PhyloCorrection   = rma.mv(z ~ PhyloCorrection,     V = z_Var,   data = Data, random = c(~ 1 | Obs_ID, ~ 1 | Study_ID, ~ 1 | Taxon, ~ 1 | Taxon_Phylo), R = list(Taxon_Phylo = PhyloTree), method = "REML")
summary(REML_z_Mod_PhyloCorrection)

# extract model statistics for Table S5
REML_z_Mod_PhyloCorrection_Name <- "Phylogenetic correction"
REML_z_Mod_PhyloCorrection_k <- REML_z_Mod_PhyloCorrection$k
REML_z_Mod_PhyloCorrection_Qdf <- REML_z_Mod_PhyloCorrection$QMdf[1]
REML_z_Mod_PhyloCorrection_QM <- REML_z_Mod_PhyloCorrection$QM
REML_z_Mod_PhyloCorrection_Qpval <- REML_z_Mod_PhyloCorrection$QMp

Model_Results_REML_z_Mod_PhyloCorrection <- data.frame(Model    = REML_z_Mod_PhyloCorrection_Name,
                                                       k        = REML_z_Mod_PhyloCorrection_k,
                                                       num_dfs  = REML_z_Mod_PhyloCorrection_Qdf,
                                                       QM       = round(REML_z_Mod_PhyloCorrection_QM,2),
                                                       Pval     = round(REML_z_Mod_PhyloCorrection_Qpval,3))
Model_Results_REML_z_Mod_PhyloCorrection


#______________________________________________________________####
# 7. Figures ####

##  * Figure 1 ####

sort_by_year <- function(data, year_column, ascending_year = TRUE) {
  if (!is.data.frame(data)) {
    stop("The input must be a data.frame")
  }
  if (!year_column %in% colnames(data)) {
    stop("The specified columns must exist in the data.frame")
  }
  year_order <- ifelse(ascending_year, 1, -1)
  sorted_data <- data[order(year_order * data[[year_column]]), ]
  return(sorted_data)
}
sorted_data <- sort_by_year(Data,"Year", ascending_year = TRUE)
sorted_data$Year_order <- cumsum(c(TRUE,diff(sorted_data$Year) !=0))
n_yearspublished <- length(unique(sorted_data$Year_order))

Table_Cum_REML_Model_All_Null_by_Year_ID <- matrix(NA,n_yearspublished,8) 

for (i in 2:(n_yearspublished)) {
  tempData <- sorted_data[-(1:length(sorted_data$Obs_ID)),]
  for (g in 1:i) {
    subdata <- subset(sorted_data, sorted_data$Year_order == g ) ## do that for all possible values and then merge (another loop that adds the dataset to a bigger dataset)
    tempData <- rbind(tempData,subdata)
  }
  try(REML_Model_All_Null_phylo <- rma.mv(r, r_Var,
                                          mod =~1,
                                          random=list(~ 1 | Taxon_Phylo, # phylo effect
                                                      ~ 1 | Taxon, # non-phylo effect
                                                      ~ 1 | Study_ID, 
                                                      ~ 1 | Obs_ID), 
                                          R = list(Taxon_Phylo = PhyloTree), #phylogenetic relatedness
                                          method="REML",
                                          test="t", # using t dist rather than z
                                          control=list(optimizer="nlminb", rel.tol=1e-8),
                                          data=tempData))
  Table_Cum_REML_Model_All_Null_by_Year_ID[i,1] <- coef(summary(REML_Model_All_Null_phylo))$estimate
  Table_Cum_REML_Model_All_Null_by_Year_ID[i,2] <- coef(summary(REML_Model_All_Null_phylo))$se
  Table_Cum_REML_Model_All_Null_by_Year_ID[i,3] <- coef(summary(REML_Model_All_Null_phylo))$tval
  Table_Cum_REML_Model_All_Null_by_Year_ID[i,4] <- coef(summary(REML_Model_All_Null_phylo))$df
  Table_Cum_REML_Model_All_Null_by_Year_ID[i,5] <- coef(summary(REML_Model_All_Null_phylo))$pval
  Table_Cum_REML_Model_All_Null_by_Year_ID[i,6] <- coef(summary(REML_Model_All_Null_phylo))$ci.lb
  Table_Cum_REML_Model_All_Null_by_Year_ID[i,7] <- coef(summary(REML_Model_All_Null_phylo))$ci.ub
  Table_Cum_REML_Model_All_Null_by_Year_ID[i,8] <- g
}

colnames(Table_Cum_REML_Model_All_Null_by_Year_ID)<- c("cum_ES","cum_ES_SE","cum_ES_t", "cum_ES_dfs","cum_ES_Pvalue", "cum_ES_l95CI", "cum_ES_u95CI", "Sort_ID")
Table_Cum_REML_Model_All_Null_by_Year_ID[1,1] <- sorted_data[1,17 ]
Table_Cum_REML_Model_All_Null_by_Year_ID[1,2] <- sqrt(sorted_data[1,18])
Table_Cum_REML_Model_All_Null_by_Year_ID[1,3] <- sorted_data[1,17 ]/sqrt(sorted_data[1,18])
Table_Cum_REML_Model_All_Null_by_Year_ID[1,6] <- sorted_data[1,17 ] - 1.96*(sqrt(sorted_data[1,18]))
Table_Cum_REML_Model_All_Null_by_Year_ID[1,7] <- sorted_data[1,17 ] + 1.96*(sqrt(sorted_data[1,18]))
Table_Cum_REML_Model_All_Null_by_Year_ID[1,8] <- 1

Table_SortedData_ES_by_Year<- as.data.frame(aggregate(r ~ sorted_data$Year, data = sorted_data, FUN = length))
colnames(Table_SortedData_ES_by_Year)<- c("Year", "N_ES")
Table_SortedData_Studies_by_Year <- as.data.frame(aggregate(r ~ sorted_data$Year + sorted_data$Study_ID, data = sorted_data, FUN = length))
colnames(Table_SortedData_Studies_by_Year)<- c("Year","Study_ID",  "N_ES")
Table_SortedData_Studies_by_Year <- as.data.frame(aggregate(N_ES ~ Table_SortedData_Studies_by_Year$Year, data = Table_SortedData_Studies_by_Year, FUN = length))
colnames(Table_SortedData_Studies_by_Year)<- c("Year", "N_ES")
Table_ES_and_Studies_by_Year <-merge(Table_SortedData_ES_by_Year,Table_SortedData_Studies_by_Year, by = "Year")
colnames(Table_ES_and_Studies_by_Year)<- c("Year", "N_ES", "N_Studies")
Table_ES_and_Studies_by_Year$Sort_ID <- seq.int(nrow(Table_ES_and_Studies_by_Year))

Cumulative_Data_by_Year <- merge(Table_ES_and_Studies_by_Year, Table_Cum_REML_Model_All_Null_by_Year_ID, by = "Sort_ID")
Cumulative_Data_by_Year <- within(Cumulative_Data_by_Year, N_ES_cum <- cumsum(N_ES))
Cumulative_Data_by_Year <- within(Cumulative_Data_by_Year, N_Studies_cum <- cumsum(N_Studies))
Cumulative_Data_by_Year <- as.data.frame(Cumulative_Data_by_Year)
full_years <- data.frame(Year = seq(min(Cumulative_Data_by_Year$Year), max(Cumulative_Data_by_Year$Year), by = 1))

Data_cumulative_by_Year <- full_years %>%
  left_join(Cumulative_Data_by_Year, by = "Year")

Data_cumulative_by_Year <- Data_cumulative_by_Year %>%
  fill(cum_ES, cum_ES_SE,cum_ES_t, cum_ES_dfs, cum_ES_Pvalue, cum_ES_l95CI, cum_ES_u95CI, N_ES_cum, N_Studies_cum,
       .direction = "down")
Data_cumulative_by_Year

Plot_SampleSize_by_Year <- ggplot() +
  geom_rect(aes(xmin=1994, xmax=2008.5, ymin=0, ymax=51), fill="grey92", alpha=1, size = 0, inherit.aes = FALSE) +
  geom_bar(data = Data_cumulative_by_Year, aes(x = Year, y = N_Studies_cum), stat = "identity", fill ="#3771c8ff", alpha = 0.7) +
  geom_bar(data = Data_cumulative_by_Year, aes(x = Year, y = N_Studies), stat = "identity", fill = "#d5e5ffff", alpha = 0.7) +
  geom_line(data = Data_cumulative_by_Year, aes(x = Year, y = N_ES_cum/3), linewidth = 0.6, color="black")+
  labs(title = "", x = "Year", y = "Number of publications") +
  scale_x_continuous(limits = c(1994, 2025), breaks=seq(1995, 2025, by=5), expand=c(0.02,0.02)) +
  scale_y_continuous(limits = c(0, 51),
                     sec.axis = sec_axis(transform=~.*3, name = "Number of effect sizes")) +
  theme(panel.border = element_blank(), 
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = c(0.11, 0.95),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", size = 14),
        axis.line.x = element_line(colour = "black", size = 1),
        axis.line.y = element_line(colour = "black", size = 1),
        axis.text.x = element_text(face = "plain", colour = "black", size = 16, angle = 0),
        axis.text.y = element_text(face = "plain", colour = "black", size = 16, angle = 0),
        axis.title.x = element_text(face = "plain", size = 16, margin = margin(r = 0, 10, 0, 0)),
        axis.title.y = element_text(face = "plain", size = 16, margin = margin(r = 10, 0, 0, 0)),
        axis.ticks = element_line(size = 1, colour = "black"),
        axis.ticks.length = unit(.3, "cm")) +
  annotate("text", x = 2001, y = 47, label = "previous meta-analysis", size = 4.5, angle=0, lineheight = 1, hjust = 0.5)
Plot_SampleSize_by_Year  

Estimate_observed <- coef(summary(REML_Model_All_Null_phylo))$estimate

Plot_ES_by_Year <- ggplot(Data_cumulative_by_Year, aes(x = Year, y = cum_ES)) +
  geom_rect(aes(xmin=1994, xmax=2008.5, ymin=-0.6, ymax=1.1), fill="grey92", alpha=1, size = 0, inherit.aes = FALSE) +
  geom_hline(yintercept=0, color='black', linetype='longdash', alpha=1.0, linewidth=0.6) +
  geom_hline(yintercept=Estimate_observed, color='black', linetype='solid', alpha=1.0, linewidth=0.6) +
  geom_errorbar(aes(ymin = cum_ES_l95CI, ymax = cum_ES_u95CI), size = 4, width = 0.0, color = "#3771c8ff", alpha = 0.7) +
  geom_point(size = 2.5, shape=21, alpha=1, fill = 'white', colour='black') +  
  labs( x = "Year", y = "Effect size (r)") +
  scale_x_continuous(limits = c(1994, 2025), breaks=seq(1995, 2025, by=5), expand=c(0.02,0.02)) +
  theme(panel.border = element_blank(), 
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = c(0.11, 0.95),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", size = 14),
        axis.line.x = element_line(colour = "black", size = 1),
        axis.line.y = element_line(colour = "black", size = 1),
        axis.text.x = element_text(face = "plain", colour = "black", size = 16, angle = 0),
        axis.text.y = element_text(face = "plain", colour = "black", size = 16, angle = 0),
        axis.title.x = element_text(face = "plain", size = 16, margin = margin(r = 0, 10, 0, 0)),
        axis.title.y = element_text(face = "plain", size = 16, margin = margin(r = 10, 0, 0, 0)),
        axis.ticks = element_line(size = 1, colour = "black"),
        axis.ticks.length = unit(.3, "cm"))
Plot_ES_by_Year

Figure_Meta_by_Year <- plot_grid(Plot_SampleSize_by_Year, Plot_ES_by_Year,
                                 labels = "AUTO",
                                 label_size = 15,
                                 hjust = 0, 
                                 vjust = 1, 
                                 align = "hv",
                                 ncol = 1,
                                 nrow = 2)
Figure_Meta_by_Year
ggsave("Figure_1.pdf", width = 7, height = 10, useDingbats=FALSE)


##  * Figure 2 ####

icon_insect <- get_phylopic(uuid = "267f1e53-333a-4dbd-943c-1031f18b569e")
icon_bird <- get_phylopic(uuid = "7af0b4cf-c6bf-4b54-879b-36d957d795df")
icon_fish <- get_phylopic(uuid = "23a7d09d-4a4d-4ad5-ad07-49a6b59a7fba")
icon_amphibian <- get_phylopic(uuid = "31763631-b109-46e5-9588-7f5de7fac7c0")
icon_reptile <- get_phylopic(uuid = "264fa655-afd7-451c-8f27-e0a9557376e6")
icon_bird <- get_phylopic(uuid = "7af0b4cf-c6bf-4b54-879b-36d957d795df")
icon_mammal <- get_phylopic(uuid = "b37d2e2f-4753-4414-b275-882bbc860e24")
icon_mammal <-flip_phylopic(icon_mammal, horizontal = TRUE, vertical = FALSE)

palette_1=c("#8C510A", "#BF812D", "#DFC27D", "#F6E8C3", "#C7EAE5", "#80CDC1", "#35978F", "#737373") # BrBG
ForestPlot_Data <- Data[, c("Obs_ID", "Study_ID", "Authors_and_Year", "Taxon",
                            "Class_Order", "PhyloCorrection", "N_r",
                            "r", "r_Var", "r_lCI", "r_uCI")]

ForestPlot_Data$Class_Order=as.factor(ForestPlot_Data$Class_Order)
ForestPlot_Data <- ForestPlot_Data[order(ForestPlot_Data$r, decreasing = TRUE),]
ForestPlot_Data <- ForestPlot_Data[order(ForestPlot_Data$Class_Order, decreasing = TRUE),]
ForestPlot_Data$Observation <- 1:nrow(ForestPlot_Data) 

REML_Model_All_Null_phylo_r_estimate     <- (coef(summary(REML_Model_All_Null_phylo_r)))$estimate
REML_Model_All_Null_phylo_r_estimate_lCI <- (coef(summary(REML_Model_All_Null_phylo_r)))$ci.lb
REML_Model_All_Null_phylo_r_estimate_uCI <- (coef(summary(REML_Model_All_Null_phylo_r)))$ci.ub
Table_REML_Model_All_Null_phylo_r_Predict <- as.data.frame(predict(REML_Model_All_Null_phylo_r))
REML_Model_All_Null_phylo_r_estimate_lPI <- Table_REML_Model_All_Null_phylo_r_Predict$pi.lb
REML_Model_All_Null_phylo_r_estimate_uPI <- Table_REML_Model_All_Null_phylo_r_Predict$pi.ub

Forest_r <- ggplot(data=ForestPlot_Data, aes(y=Observation, x=r, xmin=r_lCI, xmax=r_uCI,color=Class_Order), inherit.aes = FALSE) +
  geom_rect(aes(xmin=REML_Model_All_Null_phylo_r_estimate_lPI, xmax=REML_Model_All_Null_phylo_r_estimate_uPI, ymin=0, ymax=147), fill="gray95", alpha=0.5, linewidth = 0, inherit.aes = FALSE) +
  geom_rect(aes(xmin=REML_Model_All_Null_phylo_r_estimate_lCI, xmax=REML_Model_All_Null_phylo_r_estimate_uCI, ymin=0, ymax=147), fill="grey85", alpha=0.5, linewidth = 0, inherit.aes = FALSE) +
  geom_vline(xintercept=REML_Model_All_Null_phylo_r_estimate, color='black', linetype='solid', alpha=1.0, linewidth=0.4) +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=1.0, linewidth=0.6) +
  geom_errorbarh(height=0, size = 1.5, alpha=0.9)+
  geom_point(size = 1.4, shape=21, alpha=1, fill = 'white', colour='black')+
  scale_x_continuous(limits = c(-1.3, 1.4), expand=c(0,0)) +
  scale_y_continuous(limits = c(0, 148), breaks=seq(0, 147, by=20), expand=c(0.02,0.02)) +
  labs(title='', x=expression(paste("Effect size ", italic(r),' (95% CI)')), y = 'Study ID') +
  theme_classic()+
  theme(axis.line = element_line(lineend = "square")) +
  theme(axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14),
        axis.title=element_text(size=16),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill= "transparent"),
        plot.margin = margin(t = 0,  # Top margin
                             r = 30,  # Right margin
                             b = 5,  # Bottom margin
                             l = 5)) + # Left margin) +
  add_phylopic(img = icon_insect,    x = 1.4, y = 142, height = 6,  alpha = 0.8, fill = "black") +
  add_phylopic(img = icon_fish,      x = 1.4, y = 125, height = 4,  alpha = 0.8, fill = "black") +
  add_phylopic(img = icon_amphibian, x = 1.4, y = 108, height = 8,  alpha = 0.8, fill = "black") +
  add_phylopic(img = icon_reptile,   x = 1.4, y = 95,  height = 12, alpha = 0.8, fill = "black") +  
  add_phylopic(img = icon_bird,      x = 1.4, y = 67,  height = 8, alpha = 0.8, fill = "black") +
  add_phylopic(img = icon_mammal,    x = 1.4, y = 28,  height = 8, alpha = 0.8, fill = "black") +
  scale_colour_manual(values=palette_1, labels=NULL,guide = "none") 
Forest_r

ggsave("Figure_2.pdf", width=6.5, height=8.5, useDingbats=FALSE)


##  * Figure 3 ####

Figure_major_moderators <- plot_grid(OrchaRD_Plot_SexSel_proxy,
                                     OrchaRD_Plot_Speciation_proxy,
                                     labels = "AUTO",
                                     label_size = 15,
                                     hjust = 0, 
                                     vjust = 1, 
                                     align = "hv",
                                     ncol = 1,
                                     nrow = 2)
Figure_major_moderators
ggsave("Figure_3.pdf", width=6, height=9, useDingbats=FALSE)


##  * Figure S2 ####

Data_N <- as.data.frame(aggregate(r ~ Class + Study_ID + Obs_ID, data = Data, FUN = mean, na.rm = TRUE))
Data_N_by_Study_ID <- as.data.frame(aggregate(r ~ Data_N$Class + Data_N$Study_ID, data = Data_N, FUN = length))
colnames(Data_N_by_Study_ID) <- c("Class", "Study_ID", "N_ES")
Data_N_by_Class_ES <- as.data.frame(aggregate(N_ES ~ Class, data = Data_N_by_Study_ID, FUN = sum))
Data_N_by_Class_Studies <- as.data.frame(aggregate(N_ES ~ Class, data = Data_N_by_Study_ID, FUN = length))

Data_N_by_Class_ES$N_Studies <- Data_N_by_Class_Studies$N_ES
Data_N_by_Class_ES$order      <- c(3, 4, 8, 1, 6, 2, 7, 5)
rownames(Data_N_by_Class_ES) <- 1:nrow(Data_N_by_Class_ES)
Data_N_by_Class_ES

palette_1=c("#8C510A", "#BF812D", "#DFC27D", "#F6E8C3", "#C7EAE5", "#80CDC1", "#35978F", "#F0F0F0") # BrBG

ES_Sum <- sum(Data_N_by_Class_ES$N_ES)

DChart_N_ES <- ggplot(data = Data_N_by_Class_ES, mapping = aes(x = 1, y = N_ES, fill = as.factor(order))) +
  scale_fill_manual(values = palette_1) +
  geom_bar(stat = "identity") +
  xlim(0, 1.5) +
  coord_polar(theta = "y", direction = -1) +
  labs(x = NULL, y = NULL) +
  labs(fill = "") +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank()) +
  annotate("text", x = 0, y = 1, label = ES_Sum, size = 5, fontface = "bold", angle = 0) 
DChart_N_ES

Studies_Sum <- length(summary(as.factor(Data$Authors_and_Year))) # Number of studies

DChart_N_Studies <- ggplot(data = Data_N_by_Class_ES, mapping = aes(x = 1, y = N_Studies, fill = as.factor(order))) +
  scale_fill_manual(values = palette_1) +
  geom_bar(stat = "identity") +
  xlim(0, 1.5) +
  coord_polar(theta = "y", direction = -1) +   # polar coordinate system
  labs(x = NULL, y = NULL) +
  labs(fill = "") +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank()) +
  annotate("text", x = 0, y = 1, label = Studies_Sum, size = 5, fontface = "bold", angle = 0) 
DChart_N_Studies

Figure_Donut_ALL <- plot_grid(DChart_N_ES, DChart_N_Studies,
                              #labels = c("a","d","b","e","c","f"),
                              #                         labels = "AUTO",
                              #label_size = 15,
                              #hjust = -5, 
                              #vjust = -5, 
                              #align = "hv",
                              ncol = 2)

p_nodes <- ggtree(Tree_pruned_ALL) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab()  + xlim(-100, 550)
Phylo_Tree <- rotate(p_nodes, 55) %>% rotate(83) %>% rotate(85)
Phylo_Tree

Phylo <- ggtree(Tree_pruned_ALL, size = 0.7) +
  #ylim(-2, 51) +
  xlim(-20, 650) +
  annotate("rect", xmin = 425, xmax = 560, ymin = 54.5, ymax = 53.5, fill = palette_1[8], alpha = .8) + # Animalia
  annotate("rect", xmin = 425, xmax = 560, ymin = 53.5, ymax = 52.5, fill = palette_1[1], alpha = .8) + # Araneae
  annotate("rect", xmin = 425, xmax = 560, ymin = 52.5, ymax = 45.5, fill = palette_1[2], alpha = .8) + # Insecta
  annotate("rect", xmin = 425, xmax = 560, ymin = 45.5, ymax = 44.5, fill = palette_1[8], alpha = .8) + # Vertebrata
  annotate("rect", xmin = 425, xmax = 560, ymin = 44.5, ymax = 22.5, fill = palette_1[3], alpha = .8) + # Fishes
  annotate("rect", xmin = 425, xmax = 560, ymin = 22.5, ymax = 19.5, fill = palette_1[4], alpha = .8) + # Amphibia
  annotate("rect", xmin = 425, xmax = 560, ymin = 19.5, ymax = 13.5, fill = palette_1[5], alpha = .8) + # Reptilia
  annotate("rect", xmin = 425, xmax = 560, ymin = 13.5, ymax =  2.5, fill = palette_1[6], alpha = .8) + # Aves
  annotate("rect", xmin = 425, xmax = 560, ymin =  2.5, ymax =  0.4, fill = palette_1[7], alpha = .8) + # Mammalia
  geom_tiplab(size = 3.4, color = "black", hjust = 0.0) +
  add_phylopic(img = icon_insect,    x = 620, y = 49, height = 2.5,  alpha = 0.8, fill = "black") +
  add_phylopic(img = icon_fish,      x = 620, y = 34, height = 1.8,  alpha = 0.8, fill = "black") +
  add_phylopic(img = icon_amphibian, x = 620, y = 22, height = 3.2,  alpha = 0.8, fill = "black") +
  add_phylopic(img = icon_reptile,   x = 620, y = 16, height = 5,    alpha = 0.8, fill = "black") +  
  add_phylopic(img = icon_bird,      x = 620, y = 9,  height = 3.5,  alpha = 0.8, fill = "black") +
  add_phylopic(img = icon_mammal,    x = 620, y = 2,  height = 3.5,    alpha = 0.8, fill = "black") +
  theme_tree()

Phylo <- rotate(Phylo, 55) %>% rotate(83) %>% rotate(85)
Phylo

ggdraw(Phylo) + 
  draw_plot(DChart_N_ES,      x = -0.4, y = -0.36, scale = 0.23) +
  draw_plot(DChart_N_Studies, x = -0.2, y = -0.36, scale = 0.23) +
  annotate("text", x = 0.095, y = 0.05, label = "Number of\neffect sizes", size = 5, angle=0, lineheight = 1, hjust = 0.5) +
  annotate("text", x = 0.295, y = 0.05, label = "Number of\nstudies",       size = 5, angle=0, lineheight = 1, hjust = 0.5)

ggsave("Figure_S2.pdf", width=8, height=10, useDingbats=FALSE)


##  * Figure S3 ####

Dataset <-  data.table(Data)
years <- min(Dataset$Year):max(Dataset$Year)  # Ensure all years from 1995 to 2025 are included
classes <- unique(Dataset$Class)  # Explicitly define all classes
all_combinations <- CJ(Year = years, Class = classes, sorted = TRUE)

Dataset_by_Year_and_Class <- Dataset[, .(
  NumEffectSizes = .N, 
  NumStudies = uniqueN(Study_ID)
), by = .(Year, Class)]

sum(Dataset_by_Year_and_Class$NumEffectSizes)
sum(Dataset_by_Year_and_Class$NumStudies)

Dataset_by_Year_and_Class <- merge(all_combinations, Dataset_by_Year_and_Class, by = c("Year", "Class"), all.x = TRUE)
Dataset_by_Year_and_Class[is.na(NumEffectSizes), NumEffectSizes := 0]
Dataset_by_Year_and_Class[is.na(NumStudies), NumStudies := 0]
Dataset_by_Year_and_Class[, CumulativeEffectSizes := cumsum(NumEffectSizes), by = Class]
Dataset_by_Year_and_Class[, CumulativeStudies := cumsum(NumStudies), by = Class]
Dataset_by_Year_and_Class[, TotalCumulativeEffectSizes := sum(CumulativeEffectSizes), by = Year]
Dataset_by_Year_and_Class[, TotalCumulativeStudies := sum(CumulativeStudies), by = Year]
Dataset_by_Year_and_Class[, PropCumulativeEffectSizes := CumulativeEffectSizes / TotalCumulativeEffectSizes, by = Year]
Dataset_by_Year_and_Class[, PropCumulativeStudies := CumulativeStudies / TotalCumulativeStudies, by = Year]

Class_order <- c("Arachnida", "Insecta", "Actinopterygii", "Amphibia", "Squamata", "Aves", "Mammalia","Animalia")
Dataset_by_Year_and_Class$Class <- factor(Dataset_by_Year_and_Class$Class, levels = Class_order)
palette <- c("#8C510A", "#BF812D", "#DFC27D", "#F6E8C3", "#C7EAE5", "#80CDC1", "#35978F", "#737373")

scale_x_discrete(expand=c(0.05, 0.05))

Plot_ES_Prop_by_Year_and_Class <- ggplot(Dataset_by_Year_and_Class, aes(x = Year, y = PropCumulativeEffectSizes, fill = Class)) +
  geom_bar(stat = "identity") +
  labs( x = "Year", y = "Proportion of effect sizes") +
  scale_fill_manual(values = palette) +
  scale_x_continuous(limits = c(1994, 2025), breaks=seq(1995, 2025, by=5), expand=c(0.02,0.02)) +
  theme(panel.border = element_blank(), 
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "right",
        legend.margin = margin(5, 5, 5, 5),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", size = 11),
        axis.line.x = element_line(colour = "black", size = 1),
        axis.line.y = element_line(colour = "black", size = 1),
        axis.text.x = element_text(face = "plain", colour = "black", size = 16, angle = 0),
        axis.text.y = element_text(face = "plain", colour = "black", size = 16, angle = 0),
        axis.title.x = element_text(face = "plain", size = 16, margin = margin(r = 0, 10, 0, 0)),
        axis.title.y = element_text(face = "plain", size = 16, margin = margin(r = 10, 0, 0, 0)),
        axis.ticks = element_line(size = 1, colour = "black"),
        axis.ticks.length = unit(.3, "cm"))
Plot_ES_Prop_by_Year_and_Class

Plot_Studies_Prop_by_Year_and_Class <- ggplot(Dataset_by_Year_and_Class, aes(x = Year, y = PropCumulativeStudies, fill = Class)) +
  geom_bar(stat = "identity") +
  labs( x = "Year", y = "Proportion of studies") +
  scale_fill_manual(values = palette) +
  scale_x_continuous(limits = c(1994, 2025), breaks=seq(1995, 2025, by=5), expand=c(0.02,0.02)) +
  theme(panel.border = element_blank(), 
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "right",
        legend.margin = margin(5, 5, 5, 5),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", size = 11),
        axis.line.x = element_line(colour = "black", size = 1),
        axis.line.y = element_line(colour = "black", size = 1),
        axis.text.x = element_text(face = "plain", colour = "black", size = 16, angle = 0),
        axis.text.y = element_text(face = "plain", colour = "black", size = 16, angle = 0),
        axis.title.x = element_text(face = "plain", size = 16, margin = margin(r = 0, 10, 0, 0)),
        axis.title.y = element_text(face = "plain", size = 16, margin = margin(r = 10, 0, 0, 0)),
        axis.ticks = element_line(size = 1, colour = "black"),
        axis.ticks.length = unit(.3, "cm"))
Plot_Studies_Prop_by_Year_and_Class

Figure_Proportions_by_Year_and_Class <- plot_grid(Plot_ES_Prop_by_Year_and_Class, Plot_Studies_Prop_by_Year_and_Class,
                                                  labels = "AUTO",
                                                  label_size = 15,
                                                  hjust = 0, 
                                                  vjust = 1, 
                                                  align = "hv",
                                                  ncol = 1,
                                                  nrow = 2)
Figure_Proportions_by_Year_and_Class
ggsave("Figure_S3.pdf", width = 8, height = 10, useDingbats=FALSE)


##  * Figure S4 ####

palette_ClassPlot=c("#737373",     "#35978F",   "#80CDC1",   "#C7EAE5",   "#DFC27D",   "#BF812D") # BrBG
Class_order <- c("Animalia", "Mammalia",  "Aves", "Squamata",  "Actinopterygii", "Insecta")

OrchaRD_Plot_Class <- orchaRd::orchard_plot(REML_r_Mod_Class, 
                                            mod = "Class", group = "Study_ID", 
                                            #xlab = "Effect size (r)", angle = 0, alpha = 0.6, 
                                            xlab = expression(paste("Effect size (", italic(r),')')), angle = 0, alpha = 0.8,
                                            g = TRUE,
                                            k = TRUE,
                                            transfm = "none", 
                                            twig.size = 0.75, trunk.size = 0.7,  branch.size = 1.5,
                                            tree.order = c("Animalia", "Mammalia",  "Aves", "Squamata",  "Actinopterygii", "Insecta")) +
  theme(axis.line = element_line(lineend = "square")) +
  theme(legend.position.inside = c(0, 1.05), 
        legend.justification = c(0, 1), 
        legend.key.size = unit(1, "mm"),
        legend.direction = "horizontal", 
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 12),
        legend.key = element_blank()) +
  scale_fill_manual(values = palette_ClassPlot) + 
  scale_colour_manual(values = palette_ClassPlot) +
  scale_x_discrete(limits = Class_order) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line.x = element_line(colour = "black", linewidth = 0.7, lineend = "square"),
        axis.line.y = element_line(colour = "black", linewidth = 0.7, lineend = "square"),
        axis.text.x = element_text(face="plain", color="black", size=12, angle=0),
        axis.text.y = element_text(face="plain", color="black", size=14, angle=0),
        axis.title.x = element_text(size=14,face="plain", margin = margin(r=10,0,0,0), vjust = -1),
        axis.title.y = element_text(size=14,face="plain", margin = margin(r=10,0,0,0)),
        axis.ticks = element_line(linewidth = 0.7, colour = "black"),
        axis.ticks.length = unit(.3, "cm"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill= "transparent"),
        plot.margin = margin(t =  25,  # Top margin
                             r =  1,  # Right margin
                             b =  1,  # Bottom margin
                             l =  1)) # Left margin 
OrchaRD_Plot_Class
ggsave("Figure_S4.pdf", width=5.5, height=7.5, useDingbats=FALSE)


##  * Figure S5 ####

Figure_minor_moderators <- plot_grid(OrchaRD_Plot_SexSel_proxy_mechanism,
                                     OrchaRD_Plot_SexSel_proxy_stage,
                                     OrchaRD_Plot_SexSel_proxy_sex,
                                     OrchaRD_Plot_PhyloCorrection,
                                     labels = "AUTO",
                                     label_size = 15,
                                     hjust = 0, 
                                     vjust = 1, 
                                     align = "hv",
                                     ncol = 2,
                                     nrow = 2)
Figure_minor_moderators
ggsave("Figure_S5.pdf", width=9, height=8, useDingbats=FALSE)


##  * Figure S6 ####

REML_Model_All_Full_phylo_Moderators_Except_SE_and_Year <- rma.mv(z, z_Var,
                                                                  mod =~factor(PhyloCorrection),
                                                                  random=list(~ 1 | Taxon_Phylo, # phylo effect
                                                                              ~ 1 | Taxon, # non-phylo effect
                                                                              ~ 1 | Study_ID, 
                                                                              ~ 1 | Obs_ID), 
                                                                  R = list(Taxon_Phylo = PhyloTree), #phylogenetic relatedness
                                                                  method="REML",
                                                                  test="t", # using t dist rather than z 
                                                                  data=Data)
REML_Model_All_Full_phylo_Moderators_Except_SE_and_Year

Data$Res_REML_Model_All_Full_phylo_Moderators_Except_SE_and_Year <- resid(REML_Model_All_Full_phylo_Moderators_Except_SE_and_Year)
Data$z_SE <- sqrt(Data$z_Var)

Scatter_z_versus_z_SE <- ggplot(data = Data,
                                mapping = aes(x = z_SE, y = Res_REML_Model_All_Full_phylo_Moderators_Except_SE_and_Year)) +
  geom_hline(yintercept = 0, linetype = "longdash", color = "black", size = 0.5) +
  labs(title = "", x = "Standard Error", y = "Residual Fisher's z") +
  scale_y_continuous(limits = c(-1.3, 0.6), breaks = c(-1.2, -0.6, 0, 0.6)) + 
  scale_x_continuous(limits = c(0, 0.75), breaks = c(0, 0.2, 0.4, 0.6, 0.8)) + 
  geom_smooth(method = lm, se = TRUE, linetype = "longdash", color = "#08519C", fill = "#9ECAE1", size = 0.8) +
  geom_point(shape = 21, size = 5, fill = "grey", alpha = .6) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "none",
        axis.line.x = element_line(colour = "black", linewidth = 1, lineend = "square"),
        axis.line.y = element_line(colour = "black", linewidth = 1, lineend = "square"),
        axis.text.x = element_text(face="plain", color="black", size=16, angle=0),
        axis.text.y = element_text(face="plain", color="black", size=16, angle=0),
        axis.title.x = element_text(size=18,face="plain", margin = margin(r=0,10,0,0)),
        axis.title.y = element_text(size=18,face="plain", margin = margin(r=5,0,0,0)),
        axis.ticks = element_line(linewidth = 1, colour = "black"),
        axis.ticks.length = unit(.3, "cm"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill= "transparent"),
        plot.margin = margin(t = 1,  # Top margin
                             r = 1,  # Right margin
                             b = 1,  # Bottom margin
                             l = 1)) # Left margin
Scatter_z_versus_z_SE

Scatter_z_versus_year <- ggplot(data = Data,
                                mapping = aes(x = Year, y = Res_REML_Model_All_Full_phylo_Moderators_Except_SE_and_Year)) +
  geom_hline(yintercept = 0, linetype = "longdash", color = "black", size = 0.5) +
  labs(title = "", x = "Year of publication", y = "Residual Fisher's z") +
  scale_y_continuous(limits = c(-1.4, 0.6), breaks = c(-1.2, -0.6, 0, 0.6)) + 
  geom_smooth(method = lm, se = TRUE, linetype = "longdash", color = "#08519C", fill = "#9ECAE1", size = 0.8) +
  geom_point(shape = 21, size = 5, fill = "grey", alpha = .6) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "none",
        axis.line.x = element_line(colour = "black", linewidth = 1, lineend = "square"),
        axis.line.y = element_line(colour = "black", linewidth = 1, lineend = "square"),
        axis.text.x = element_text(face="plain", color="black", size=16, angle=0),
        axis.text.y = element_text(face="plain", color="black", size=16, angle=0),
        axis.title.x = element_text(size=18,face="plain", margin = margin(r=0,10,0,0)),
        axis.title.y = element_text(size=18,face="plain", margin = margin(r=5,0,0,0)),
        axis.ticks = element_line(linewidth = 1, colour = "black"),
        axis.ticks.length = unit(.3, "cm"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill= "transparent"),
        plot.margin = margin(t = 1,  # Top margin
                             r = 1,  # Right margin
                             b = 1,  # Bottom margin
                             l = 1)) # Left margin
Scatter_z_versus_year

Figure_PubBias_vertical <- plot_grid(Scatter_z_versus_z_SE, Scatter_z_versus_year,
                                     labels = "AUTO",
                                     label_size = 15,
                                     hjust = 0, 
                                     vjust = 1, 
                                     align = "hv",
                                     ncol = 1,
                                     nrow = 2)
Figure_PubBias_vertical
ggsave("Figure_S6.pdf", width=6, height=9, useDingbats=FALSE)


#______________________________________________________________####
# 8. Tables ####

##  * Table 1 ####

Table_Null_Models_r <- rbind(Model_All_Null_trad_r_Table,
                             REML_Model_All_Null_phylo_r_Table,
                             REML_Model_All_Null_phylo_r_cons_Table)
Table_Null_Models_r <- as.data.table(Table_Null_Models_r)
Table_Null_Models_r_transposed <- t(Table_Null_Models_r)
Table_Null_Models_r_transposed <- as.data.frame(Table_Null_Models_r_transposed)
Table_Null_Models_r_transposed <- cbind(RowNames = rownames(Table_Null_Models_r_transposed), Table_Null_Models_r_transposed)
rownames(Table_Null_Models_r_transposed) <- NULL  # Remove row names
colnames(Table_Null_Models_r_transposed) <- Table_Null_Models_r_transposed[1, ]  # Set first row as column names
Table_Null_Models_r_transposed <- Table_Null_Models_r_transposed[-1, ]  # Remove the first row since it's now the header
print(Table_Null_Models_r_transposed)
write.csv(Table_Null_Models_r_transposed, "Table_1.csv", row.names = FALSE)


##  * Table 2 ####

Table_Mod_Models_r <- rbind(Model_Results_REML_r_Mod_Class,
                            Model_Results_REML_r_Mod_Speciation_proxy,
                            Model_Results_REML_r_Mod_SexSel_proxy,
                            Model_Results_REML_r_Mod_SexSel_proxy_sex,
                            Model_Results_REML_r_Mod_SexSel_proxy_mechanism,
                            Model_Results_REML_r_Mod_SexSel_proxy_stage,
                            Model_Results_REML_r_Mod_PhyloCorrection)
Table_Mod_Models_r
write.csv(Table_Mod_Models_r, "Table_2.csv", row.names = FALSE)


##  * Table S1 ####

factors <- c("Class", 
             "Speciation_proxy",
             "SexSel_proxy", 
             "SexSel_proxy_mechanism",
             "SexSel_proxy_stage",
             "PhyloCorrection",
             "SexSel_proxy_sex")

results_list <- lapply(factors, function(factor) {
  if (!(factor %in% colnames(Data))) return(NULL)
  df <- Data %>%
    group_by(.data[[factor]]) %>%
    summarise(Number_of_studies = n_distinct(Study_ID),
              Number_of_observations = n(), .groups = "drop") %>%
    mutate(Moderator = factor, Level = as.character(.data[[factor]]))
  
  df <- df %>%
    mutate(Proportion_of_observations = round(Number_of_observations / sum(Number_of_observations), 2)) %>%
    arrange(factor(Level, levels = c(setdiff(unique(Level), c("Other", "Both")), "Other", "Both"))) %>%
    select(Moderator, Level, Number_of_studies, Number_of_observations, Proportion_of_observations)
  
  return(df)
})

final_results <- bind_rows(results_list)
print(final_results)
write.csv(final_results, "Table_S1.csv", row.names = FALSE)


##  * Table S2 ####

get_unique_values <- function(x) {
  unique_vals <- unique(na.omit(x))
  paste(unique_vals, collapse = ", ")
}

SummaryTable_SexSelProxy <- Data %>%
  group_by(SexSel_measure) %>%
  summarise(
    Total_Cases = n(),
    SexSel_proxy_mechanism = get_unique_values(SexSel_proxy_mechanism),
    SexSel_proxy_stage = get_unique_values(SexSel_proxy_stage),
    SexSel_proxy = get_unique_values(SexSel_proxy),
    SexSel_proxy_sex = get_unique_values(SexSel_proxy_sex),
    .groups = 'drop'
  ) %>%
  arrange(SexSel_proxy) %>%  # Sort by category
  select(
    SexSel_measure,
    Total_Cases,
    SexSel_proxy,  # Move this right after proxy
    SexSel_proxy_mechanism,
    SexSel_proxy_stage,
    SexSel_proxy_sex
  ) %>%
  rename(
    "Sexual selection measure"= SexSel_measure,
    "Sexual selection proxy" = SexSel_proxy,
    "Cases" = Total_Cases,
    "sexual selection mechanism" = SexSel_proxy_mechanism,
    "Mating stage"=  SexSel_proxy_stage,
    "Sex"= SexSel_proxy_sex
  )

print(SummaryTable_SexSelProxy)
write.csv(SummaryTable_SexSelProxy, "Table_S2.csv", row.names = FALSE)


##  * Table S4 ####

Table_Null_Models_z <- rbind(Model_All_Null_trad_z_Table,
                             REML_Model_All_Null_phylo_z_Table,
                             REML_Model_All_Null_phylo_z_cons_Table)
Table_Null_Models_z <- as.data.table(Table_Null_Models_z)
Table_Null_Models_z_transposed <- t(Table_Null_Models_z)
Table_Null_Models_z_transposed <- as.data.frame(Table_Null_Models_z_transposed)
Table_Null_Models_z_transposed <- cbind(RowNames = rownames(Table_Null_Models_z_transposed), Table_Null_Models_z_transposed)
rownames(Table_Null_Models_z_transposed) <- NULL  # Remove row names
colnames(Table_Null_Models_z_transposed) <- Table_Null_Models_z_transposed[1, ]  # Set first row as column names
Table_Null_Models_z_transposed <- Table_Null_Models_z_transposed[-1, ]  # Remove the first row since it's now the header
print(Table_Null_Models_z_transposed)
write.csv(Table_Null_Models_z_transposed, "Table_S4.csv", row.names = FALSE)


##  * Table S5 ####

Table_Mod_Models_z <- rbind(Model_Results_REML_z_Mod_Class,
                            Model_Results_REML_z_Mod_Speciation_proxy,
                            Model_Results_REML_z_Mod_SexSel_proxy,
                            Model_Results_REML_z_Mod_SexSel_proxy_sex,
                            Model_Results_REML_z_Mod_SexSel_proxy_mechanism,
                            Model_Results_REML_z_Mod_SexSel_proxy_stage,
                            Model_Results_REML_z_Mod_PhyloCorrection)
Table_Mod_Models_z
write.csv(Table_Mod_Models_z, "Table_S5.csv", row.names = FALSE)
