################################################################################
###     Continuity in morphological disparity in tropical reef fishes        ###
###                       across evolutionary scales                         ###
### Dr. Giulia Donati, giulia.donati@usys.ethz.ch &                          ###
### Dr Camille Albouy calbouy@ethz.ch                                        ###
###                                                                          ###
################################################################################
###                  Compile the Inter-specific dataset                      ###
################################################################################

rm(list=ls())

### Loading of the Inter data 
Inter_data  <- readRDS("results/TC_dataset_landmarks_metadata.RDS")$data
Inter_data <-Inter_data[,-c(3,4)]
Inter_data <- Inter_data[-which(Inter_data$Family=="Pinguipedidae"),]

### Inter-specific variability
res_sp_family <- Inter_data %>% group_by(Family, Genus) %>% summarise_all(median) #median morphological trait value per genus
res_2 <- split(res_sp_family, res_sp_family$Family)#split by species

### Caclualtion of the variability for each family
Var_Res <- do.call(rbind,lapply(res_2, function (i) {
  Combo <- combn(1:nrow(i),2)
  y <- as.data.frame(i)
  
  sapply(3:ncol(y), function(Trait_nm){
    Num <- mean(sapply(1:ncol(Combo),function(z){abs(y[Combo[1,z],Trait_nm]-y[Combo[2,z],Trait_nm])}))
    Denum <- mean(y[,3],na.rm=T)
    Num/Denum
  }) 
}))
colnames(Var_Res) <- colnames(res_2[[1]])[3:ncol(res_2[[1]])]

fam_sp_n <- data.frame(table(Inter_data$Family))
rownames(fam_sp_n) <- fam_sp_n$Var1
colnames(fam_sp_n) <- c("Family","nb_species")

yop <- data.frame(do.call(rbind, lapply(res_2,dim)))
colnames(yop) <- c("nb_genus","col")

fam_sp_n <- merge(fam_sp_n,yop, by="row.names")
fam_sp_n <- fam_sp_n[,-c(1,5)]
rownames(fam_sp_n) <- fam_sp_n$Family

### Adding sample size per species and location
df <- merge(Var_Res,fam_sp_n,by="row.names")
rownames(df) <- df[,1]
df_observed <- df[,-c(1,29)]#removing excess columns

################################################################################
# Null model: randomisation of the genus
################################################################################

### Caluclation of the variance and standard deviation for each family
res_family_variance <- Inter_data %>% group_by(Family) %>%  summarise_all(var) 
res_family_sdev <- Inter_data %>% group_by(Family) %>%  summarise_all(sd) 

### Preparation of the randomisation
list_999 <- lapply(1:999,function(j){Inter_data})
# treat each list element: I split by genus species, sample the site and re-assemble the data
for (a in 1:length(list_999)){
  
  split_1 <- split(list_999[[a]], list_999[[a]]$Family)#split by species
  test <- list()
  for(i in 1:length(split_1)){
    test[[i]]<-split_1[[i]]
    test[[i]]$Genus <- sample(test[[i]]$Genus)
  }
  
  list_999[[a]] <- do.call(rbind,test)
}

### Random inter-specific variability
set.seed(150)
dataframe_var_random <- list()

for(x in 1:999){
  cat("x=",x,"\n")
  
  # median of morphological traits across locations of the WIO
  res_sp_genus <- list_999[[x]] %>% group_by(Family, Genus) %>%  summarise_all(median)
  res_2 <- split(res_sp_genus, res_sp_genus$Family)#split by species
  
  # inter-specific variability ----------------------------------------
  
  REs <- do.call(rbind,lapply(res_2, function (i) {
    Combo <- combn(1:nrow(i),2)
    y <- as.data.frame(i)
    sapply(3:ncol(y), function(Trait_nm){
      Num <- mean(sapply(1:ncol(Combo),function(z){abs(y[Combo[1,z],Trait_nm]-y[Combo[2,z],Trait_nm])}))
      Denum <- mean(y[,3],na.rm=T)
      Num/Denum
    }) 
  })) #mean difference pairwise combination of locations/ mean all median values per location
  colnames(REs) <- colnames(res_2[[1]])[3:ncol(res_2[[1]])]
  REs <- REs[complete.cases(REs), ]
  
  fam_sp_n <- data.frame(table(Inter_data$Family))
  rownames(fam_sp_n) <- fam_sp_n$Var1
  colnames(fam_sp_n) <- c("Family","nb_species")
  
  #adding sample size per species and location
  df2 <- merge(REs,fam_sp_n,by="row.names")
  rownames(df2) <- df2[,1]
  dataframe_var_random[[x]] <- df2[,-c(1,29)]#removing excess columns
}

################################################################################
# join all in a list: Inter_data_3 used in next scripts
################################################################################

Inter_data_set <- list(landmark= Inter_data$landmark,
                       metadata=Inter_data$data,
                       inter_variability_obs=df_observed,
                       inter_variance_family=res_family_variance,
                       inter_sd_family=res_family_sdev,
                       inter_variability_random=dataframe_var_random)

#Save to add the random site list 
saveRDS(Inter_data_set,"results/Inter_specific_dataset.RDS")

################################################################################ 