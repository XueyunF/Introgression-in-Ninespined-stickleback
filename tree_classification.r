
# Script used for classification of individual ML tree
########################
#definition of individual (14 ingroup), lineage assignment, e.g. WL, EL, admixed (test individuals)
# ind1 - WL
# ind2 - EL
# ind3 - EL
# ind4 - EL
# ind5 - least admixed WL individual
# ind 5 - 14 Test
########################
library(ape)
library(phangorn)
library(phytools)

#read trees in
setwd("path_to_trees")
t=list.files(path = "path_to_trees",pattern = "RAxML_bestTree.*.phy")
ntrees<-length(t)

treesraw<- vector("list",ntrees)
for (i in 1:ntrees){
  trees2<-read.tree(paste0("path_to_trees",t[i]))
  treesraw[[i]] <-list(trees2)
  }
alltreesraw <-do.call(rbind,treesraw)
class(alltreesraw) <- "multiPhylo"

#align trees
ntreesraw<-lapply(alltreesraw, compute.brlen, keep.multi = TRUE)
class(ntreesraw) <- "multiPhylo"

# drop tips that not relevant to the test individual to reduce complexty
# for example test for ind6 (adxmixed)
Ti <- ind6
ntreesraw<-lapply(ntreesraw,drop.tip,tip=c("ind5","ind7","ind8","ind9","ind10","ind11","ind12","ind13","ind14"))
class(ntreesraw) <- "multiPhylo"

#function to categorize trees
category <- lapply(ntreesraw, function(y){
  if (is.monophyletic(y, c(Ti,ind1)) ==TRUE & is.monophyletic(y, c(ind2,ind3,ind4)) ==TRUE){
    print("WL_branchoff") #trees that show monophyly of Ti and WL(ind1) and monophyly within EL (ind2,ind3,ind4) are claasified as WL monophyly(or WL branch-off)
  } 
  else if (!is.monophyletic(y, c(Ti,ind1)) ==TRUE & is.monophyletic(y, c(Ti,ind2,ind3,ind4)) ==TRUE){
    print("EL_branchoff") #trees that *NOT* show monophyly of Ti and WL(ind1) but monophyly between Ti (test individual) and EL (ind2,ind3,ind4) are claasified as EL monophyly(or EL branch-off)
  }
    else print("unknown")
})

# test if the ind1 is admixed with EL
Ti <- ind1
WL <- ind5 #least EL admixed WL sample

category_ind1 <- lapply(ntreesraw, function(y){
  if (is.monophyletic(y, c(Ti,WL)) ==TRUE & is.monophyletic(y, c(ind2,ind3,ind4)) ==TRUE){
    print("WL_branchoff") #trees that show monophyly of Ti and WL(ind1) and monophyly within EL (ind2,ind3,ind4) are claasified as WL monophyly(or WL branch-off)
  } 
  else if (!is.monophyletic(y, c(Ti,WL)) ==TRUE & is.monophyletic(y, c(Ti,ind2,ind3,ind4)) ==TRUE){
    print("EL_branchoff") #trees that *NOT* show monophyly of Ti and WL(ind1) but monophyly between Ti (test individual) and EL (ind2,ind3,ind4) are claasified as EL monophyly(or EL branch-off)
  }
    else print("unknown")
})

#check if EL individuals are admixed
# ind2 for example
Ti <- ind2
WL <- ind1
EL <- c(ind3,ind4) 

category_ind2 <- lapply(ntreesraw, function(y){
  if (is.monophyletic(y, c(Ti,WL)) ==TRUE & is.monophyletic(y, c(EL)) ==TRUE){
    print("WL_branchoff") #trees that show monophyly of Ti and WL(ind1) and monophyly within EL (ind2,ind3,ind4) are claasified as WL monophyly(or WL branch-off)
  } 
  else if (!is.monophyletic(y, c(Ti,WL)) ==TRUE & is.monophyletic(y, c(Ti,EL)) ==TRUE){
    print("EL_branchoff") #trees that *NOT* show monophyly of Ti and WL(ind1) but monophyly between Ti (test individual) and EL (ind2,ind3,ind4) are claasified as EL monophyly(or EL branch-off)
  }
    else print("unknown")
})
