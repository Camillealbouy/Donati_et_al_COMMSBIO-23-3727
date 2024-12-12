################################################################################
###     Continuity in morphological disparity in tropical reef fishes        ###
###                       across evolutionary scales                         ###
### Dr. Giulia Donati, giulia.donati@usys.ethz.ch &                          ###
### Dr Camille Albouy calbouy@ethz.ch                                        ###
###                                                                          ###
# Assembling morphological datasets                                          ###
# ##############################################################################

###############################################################################
# Part 1) Computation of the morphological variability across the WIO and the CA. 
###############################################################################
library(dplyr)
set.seed(150)

# load data -----------------------------------------------------------------
BigSp_list <- readRDS("results/intra_sp_morpho_traits.RDS") #the WIO,CA metadata

# median of morphological traits across locations of the WIO
res_sp_regions <- BigSp_list$data %>% group_by(genus_species, Site2) %>%  summarise_all(median)
res_2 <- split(res_sp_regions, res_sp_regions$genus_species)#split by species

# WIO intra-specific variability ----------------------------------------

REs <-  do.call(rbind,lapply(res_2, function (i) {
  i <- data.frame(i)
  Combo <- combn(1:nrow(i),2)
  
  sapply(3:ncol(i), function(Trait_nm){
    Num <- mean(sapply(1:ncol(Combo),function(z){abs(i[Combo[1,z],Trait_nm]-i[Combo[2,z],Trait_nm])}))
    Denum <- mean(i[,3],na.rm=T)
    Num/Denum
  }) 
})) #mean difference pairwise combination of locations/ mean all median values per location
colnames(REs) <- colnames(res_2[[1]])[3:ncol(res_2[[1]])]
REs <- REs[complete.cases(REs), ]

individuals_sp_n <- subset(as.data.frame(table(BigSp_list$data$genus_species)), Freq > 0)
colnames(individuals_sp_n) <- c("genus_species", "n")
rownames(individuals_sp_n) <- individuals_sp_n$genus_species

#adding sample size per species and location
dataframe_var <- merge(REs,individuals_sp_n,by="row.names")
rownames(dataframe_var) <- dataframe_var[,1]
observed_variability  <- dataframe_var[,-c(1,29)]

###############################################################################
# Part 1.1) Computation of the morphological variance and standard deviation across the WIO 
###############################################################################

# Variance and Stdev
res_sp_regions_variance <- BigSp_list$data %>% group_by(genus_species) %>%  summarise_all(var) 
res_sp_regions_sdev <- BigSp_list$data %>% group_by(genus_species) %>%  summarise_all(sd) 

###############################################################################
#
# Randomisation of the population within species - for the variability calcualation
# 
###############################################################################

list_999 <- lapply(1:999,function(j){BigSp_list$data})

# treat each list element: I split by genus species, sample the site and re-assemble the data
for (a in 1:length(list_999)){
  split_1 <- split(list_999[[a]], list_999[[a]]$genus_species) #split by species
  test <- list()
  for(i in 1:length(split_1)){
    test[[i]] <- split_1[[i]]
    test[[i]]$Site2 <- sample(test[[i]]$Site2)
  }
  list_999[[a]] <- do.call(rbind,test)
} # end of a

res_2 <- list()
dataframe_var_random <- list()

for(x in 1:999){
  cat("x=",x,"\n")
  # median of morphological traits across locations of the WIO
  res_sp_regions <- list_999[[x]] %>% group_by(genus_species, Site2) %>%  summarise_all(median)
  res_2 <- split(res_sp_regions, res_sp_regions$genus_species)#split by species
 
  # WIO intra-specific variability ----------------------------------------
  REs<-  do.call(rbind,lapply(res_2,function(i){
    Combo <- combn(1:nrow(i),2) # a faire varier
    
    y <- as.data.frame(i)
    sapply(3:ncol(y), function(Trait_nm){
      Num <- mean(sapply(1:ncol(Combo),function(z){abs(y[Combo[1,z],Trait_nm]-y[Combo[2,z],Trait_nm])}))
      Denum <- mean(y[,3],na.rm=T)
      Num/Denum
    }) 
  })) #mean difference pairwise combination of locations/ mean all median values per location
  colnames(REs) <- colnames(res_2[[1]])[3:ncol(res_2[[1]])]
  REs <- REs[complete.cases(REs), ]
  
  individuals_sp_n <- data.frame(table(list_999[[x]]$genus_species))
  individuals_sp_n <- individuals_sp_n[which(individuals_sp_n$Freq>0),]
  
  rownames(individuals_sp_n) <- individuals_sp_n$Var1
  colnames(individuals_sp_n) <- c("genus_species","n")
  
  #adding sample size per species and location
  df <- merge(REs,individuals_sp_n,by="row.names")
  rownames(df) <- df[,1]
  dataframe_var_random[[x]] <- df[,-c(1,29)]#removing excess columns
}

### join all in a list: micro_data_3 used in next scripts

Intra_specific <- list(landmark=BigSp_list$landmark,metadata=BigSp_list$data, 
                   intra_variability_pop=observed_variability,
                   intra_variance_species=res_sp_regions_variance,
                   intra_sd_species=res_sp_regions_sdev,
                   intra_variability_pop_random=dataframe_var_random)

saveRDS(Intra_specific,"results/Intra_specific_dataset.RDS") #micro_data_3
rm(list=ls())


######################## end of script #########################################
################################################################################ 