################################################################################
###     Continuity in morphological disparity in tropical reef fishes        ###
###                       across evolutionary scales                         ###
### Dr. Giulia Donati, giulia.donati@usys.ethz.ch &                          ###
### Dr Camille Albouy calbouy@ethz.ch                                        ###
###  Assembly of the intra-specific data set                                 ###  
### Part 1) assembly of the landmark data for the WIO and the CA and         ###
### matching to the metadata containing information on species id, site, etc.###
### Part 2) creation of the morphological traits data including length       ###
### (euclidean distance between points) and ratios                           ###
### Part 3) adding the morphological traits to the metadata.                 ###
################################################################################

###libraries loading
lib_vect <- c("geomorph","rgl","ape","shapes","abind","FD","rredlist",
              "dplyr","reshape2","tibble","retistruct","scales")
sapply(lib_vect,library,character.only=TRUE)

### Function loading
source("scripts/functions/morphological_functions.R")

################################################################################
###                      Part 1. Data of WIO assembly                       ####
################################################################################ 
### selection of target species retained for anaylsis.
database <- read.csv2("data/intra_specific_database.csv", h=T,sep=",")
database$Sample_ID <- gsub("Morpho","MORPHO",as.character(database$Sample_ID))
database[database$Sample_ID=="17_MF_038","genus_species"] <- "Dascyllus_trimaculatus"#species id error

### corrections Monotaxis grandoculis
heterodon <- read.table("data/Correction_pop_map_Mon_h_data.txt",h=T)
heterodon <- gsub("Mon_g_","", as.character(heterodon$Sample_ID))
heterodon <- heterodon[-grep("TR",heterodon)]
levels(database$genus_species) <- c(levels(database$genus_species), "Monotaxis_heterodon")
database$genus_species[which(database$Sample_ID%in%heterodon)] <- "Monotaxis_heterodon"

### Corrections for other species
samples_to_remove <- database[which(database$genus_species%in%"Canthigaster_valentini"),"Sample_ID"]
samples_to_remove_l_bengalensis <-c("17_SC_058","17_SC_059","17_SC_060","17_SC_061","17_SC_062","17_SC_063",
                                    "17_SC_064","17_SC_065","17_SC_066","17_SC_067","17_SC_068","17_SC_069") 
#of Lutjanus - bengalensis removal and
samples_to_remove <- c(samples_to_remove,"16_MY_178", "16_MY_761", "17_MV_332","17_MF_205",samples_to_remove_l_bengalensis) #205, 761 is juv Naso; 332 bad angle Naso (large adult); 178 M.violacea outlier
database <- database[-which(database$Sample_ID%in%samples_to_remove),]

# upload of tps files - selection of landmarks (matching Claverie and removal of semi-landmarks)
tps_load <- list.files(path="data/TPS_intraspecific_dataset/")
All_land_list <- lapply(1:length(tps_load),function(i) {
  readland.tps(paste0("data/TPS_intraspecific_dataset/",tps_load[i]), specID="IMAGE",readcurves=T,warnmsg=T)})
All_tps <- abind::abind(All_land_list[1:length(All_land_list)],along=3) #Combine datasets of shape coordinates
All_tps <- All_tps[-c(4,21,22,23,24,25,26,27,28,29,30,31,32,33,34,18,19),,]#removing semi landmarks for now
All_tps <- All_tps[,,-which(dimnames(All_tps)[[3]]%in%samples_to_remove)] # removal of samples
All_tps <- All_tps[-10,,] # removal of samples

# the database and the tps data need to have the same order of sample ID:
dimnames(All_tps)[[3]] <- gsub("Morpho","MORPHO", as.character(dimnames(All_tps)[[3]]))# correction of dimnames
database <- database[which(database$Sample_ID%in% dimnames(All_tps)[[3]]),]

# check whether TPS are in the same order as the sample ID of the database.
database$Sample_ID <- as.character(database$Sample_ID)
database <- database[match(dimnames(All_tps)[[3]],database$Sample_ID),] 

### gpagen and visualisation
All_tps_WIO <- gpagen(All_tps[,,],ProcD=T) # Procrustes analysis corrects for scale, rotation and translation
#plotAllSpecimens(All_tps_WIO$coords)# gpa-aligned landmark data

All_tps <- All_tps_WIO$coords # replacement of non gpa aligned coordinates with superimposed ones
#plotOutliers(All_tps) # outliers are exclusively Zanclidae so I kept them !

### plotTangentspace - micro with 14 landmarks matching Claverie et al. 2014
All_tps14 <- All_tps[-c(6,16),,]

# assigning a color vector to each species
dimnames(All_tps14)[[3]] <- database$genus_species

# Generate a color palette from cyan to blue
cool <- substr(colorRampPalette(c("cyan", "blue"))(25), 1, 7)

# Assign colors to groups based on sample IDs (from the 3rd dimension of All_tps14)
groups <- setNames(cool[as.factor(dimnames(All_tps14)[[3]])], dimnames(All_tps14)[[3]])

#plotTangentSpace_tot(All_tps14,gp=dimnames(All_tps14)[[3]],col_sp =unique(groups),groups=groups,
#                     legend=unique(dimnames(All_tps14)[[3]]),label=dimnames(All_tps14)[[3]],txt.cex=0.2)

#####################################################################################
### Part 2. Creation of the morphological traits data including length and ratios ###
#####################################################################################
#  to visualise the traits see script Traits_visualization on fish
# creation of additional points (intersection points) that will serve to create different length measures

back <-do.call(rbind,lapply(1:dim(All_tps)[3],function(i){
  line.line.intersection(P1=c(max(All_tps[,1,i]), max(All_tps[,2,i])), 
    P2=c(max(All_tps[,1,i]),min(All_tps[,2,i])),P3=c(min(All_tps[,1,i]),min(All_tps[,2,i])), 
    P4=c(min(All_tps[,1,i])+3,min(All_tps[,2,i])),interior.only=F)})) 

front <- do.call(rbind,lapply(1:dim(All_tps)[3],function(i){
  line.line.intersection(P1=c(min(All_tps[,1,i]),min(All_tps[,2,i])), 
    P2=c(min(All_tps[,1,i]),max(All_tps[,2,i])),P3=c(back[i,1],back[i,2]), 
    P4=c(back[i,1]-3,back[i,2]),interior.only=F)})) 

top <- do.call(rbind,lapply(1:dim(All_tps)[3],function(i){ 
  line.line.intersection(P1=c(min(All_tps[,1,i]),max(All_tps[,2,i])), 
    P2=c(max(All_tps[,1,i]),max(All_tps[,2,i])), P3=c(back[i,1],back[i,2]), 
    P4=c(back[i,1],back[i,2]+5),interior.only=F)}))

hl <-do.call(rbind,lapply(1:dim(All_tps)[3],function(i){
  line.line.intersection(P1=c(front[i,1],front[i,2]),P2=c(back[i,1],back[i,2]), 
    P3=c(All_tps[12,1,i],All_tps[12,2,i]),P4=c(All_tps[12,1,i],All_tps[12,2,i]-3),
    interior.only=F)}))

pOD <- do.call(rbind,lapply(1:dim(All_tps)[3],function(i){
   line.line.intersection(P1=c(front[i,1],front[i,2]),P2=c(back[i,1],back[i,2]), 
        P3=c(All_tps[14,1,i],All_tps[14,2,i]),P4=c(All_tps[14,1,i],All_tps[14,2,i]-3),
        interior.only=F)}))

pAL <- do.call(rbind,lapply(1:dim(All_tps)[3],function(i){
  line.line.intersection(P1=c(front[i,1],front[i,2]),P2=c(back[i,1],back[i,2]), 
      P3=c(All_tps[9,1,i],All_tps[9,2,i]),P4=c(All_tps[9,1,i],All_tps[9,2,i]-3),
      interior.only=F)}))

prePBL <- do.call(rbind,lapply(1:dim(All_tps)[3],function(i){
  line.line.intersection(P1=c(front[i,1],front[i,2]),P2=c(back[i,1],back[i,2]), 
      P3=c(All_tps[11,1,i],All_tps[11,2,i]),P4=c(All_tps[11,1,i],All_tps[11,2,i]-5),
      interior.only=F)}))

dl_start <- do.call(rbind,lapply(1:dim(All_tps)[3],function(i){
  line.line.intersection(P1=c(front[i,1],front[i,2]),P2=c(back[i,1],back[i,2]), 
     P3=c(All_tps[2,1,i],All_tps[2,2,i]),P4=c(All_tps[2,1,i],All_tps[2,2,i]-5),
     interior.only=F)}))

dl_end <- do.call(rbind,lapply(1:dim(All_tps)[3],function(i){
  line.line.intersection(P1=c(front[i,1],front[i,2]),P2=c(back[i,1],back[i,2]), 
    P3=c(All_tps[4,1,i],All_tps[4,2,i]),P4=c(All_tps[4,1,i],All_tps[4,2,i]-5),
    interior.only=F)}))

af_end<- do.call(rbind,lapply(1:dim(All_tps)[3],function(i){ 
  line.line.intersection(P1=c(front[i,1],front[i,2]),P2=c(back[i,1],back[i,2]), 
    P3=c(All_tps[8,1,i],All_tps[8,2,i]),P4=c(All_tps[8,1,i],All_tps[8,2,i]-5),
    interior.only=F)}))

ped_bottom <- do.call(rbind,lapply(1:dim(All_tps)[3],function(i){ 
  line.line.intersection(P1=c(top[i,1],top[i,2]),P2=c(back[i,1],back[i,2]), 
    P3=c(All_tps[8,1,i],All_tps[8,2,i]),P4=c(All_tps[8,1,i]+5,All_tps[8,2,i]),
    interior.only=F)}))

ped_top <- do.call(rbind,lapply(1:dim(All_tps)[3],function(i){
  line.line.intersection(P1=c(top[i,1],top[i,2]),P2=c(back[i,1],back[i,2]), 
   P3=c(All_tps[4,1,i],All_tps[4,2,i]),P4=c(All_tps[4,1,i]+5,All_tps[4,2,i]), 
   interior.only = FALSE)}))

# euclidean distances between landmarks and/or intersection points

sp_HL <- sapply(1:dim(All_tps)[3],function(i){dist(rbind(front[i,],hl[i,]),method="euclidean")})
sp_ED <- sapply(1:dim(All_tps)[3],function(i){dist(rbind(All_tps[14,,i],All_tps[15,,i]),method="euclidean")})
sp_PreOD <- sapply(1:dim(All_tps)[3],function(i){dist(rbind(front[i,],pOD[i,]), method="euclidean")})
sp_ED_ratio <- sp_ED/sp_HL #eye diameter/head length
sp_PreOD_ratio <- sp_PreOD/sp_HL #pre_eye/HL 
sp_JL <-  sapply(1:dim(All_tps)[3],function(i){dist(rbind(All_tps[1,,i],All_tps[13,,i]),method="euclidean")})
sp_JL_ratio <- sp_JL/sp_HL #jaw length/head length
sp_SL <-  sapply(1:dim(All_tps)[3],function(i){dist(rbind(front[i,],back[i,]),method="euclidean")})
sp_HL_ratio <- sp_HL/sp_SL #head lenght/standard length
sp_SL <-  sapply(1:dim(All_tps)[3],function(i){dist(rbind(front[i,],back[i,]),method="euclidean")})
sp_BH <-  sapply(1:dim(All_tps)[3],function(i){dist(rbind(top[i,],back[i,]),method="euclidean")})
sp_ELO_ratio <- sp_SL/sp_BH

# tail
sp_PH <- sapply(1:dim(All_tps)[3],function(i){dist(rbind(ped_top[i,],ped_bottom[i,]), method = "euclidean")})
sp_PL <- sapply(1:dim(All_tps)[3],function(i){dist(rbind(dl_end[i,],back[i,]),method="euclidean")})
sp_PH_ratio <- sp_PH/sp_BH #peduncle heigh/ body height
sp_PELO_ratio <- sp_PL/sp_PH

# fins
sp_DL <- sapply(1:dim(All_tps)[3],function(i){dist(rbind(dl_start[i,],dl_end[i,]),method="euclidean")})
sp_PreDL <- sapply(1:dim(All_tps)[3],function(i){dist(rbind(front[i,],dl_start[i,]),method="euclidean")})
sp_DL_ratio <- sp_DL/sp_SL #dorsal fin length / SL
sp_PreDL_ratio <-sp_PreDL/sp_SL

sp_AL <- sapply(1:dim(All_tps)[3],function(i){dist(rbind(pAL[i,],af_end[i,]),method="euclidean")})
sp_Pre_AL <- sapply(1:dim(All_tps)[3],function(i){dist(rbind(front[i,],pAL[i,]), method="euclidean")})
sp_AL_ratio <- sp_AL/sp_SL
sp_Pre_AL_ratio <- sp_Pre_AL/sp_SL

sp_PBL <- sapply(1:dim(All_tps)[3],function(i){dist(rbind(All_tps[11,,i],All_tps[10,,i]),method="euclidean")})
sp_PrePBL <- sapply(1:dim(All_tps)[3],function(i){dist(rbind(front[i,],prePBL[i,]),method="euclidean")})
sp_PrePBL_ratio <- sp_PrePBL/sp_SL
sp_PBL_ratio <- sp_PBL/sp_SL

################################################################################
################  Adding the morphological traits to the metadata ##############
################################################################################ 

trait_db <- cbind(EyeD=sp_ED,JawL=sp_JL,HeadL=sp_HL,PreED=sp_PreOD,PeduncleH=sp_PH,
                  PeduncleL=sp_PL,DorsalL=sp_DL,AnalL=sp_AL,PectoralBL=sp_PBL,Pre_Anal=sp_Pre_AL,
                  SL=sp_SL,BH=sp_BH,Pre_PectoralBL=sp_PrePBL,Pre_DorsalL=sp_PreDL,ratio_ELO=sp_ELO_ratio
                  ,ratio_PeduncleELO=sp_PELO_ratio,ratio_EyeD=sp_ED_ratio,ratio_JawL=sp_JL_ratio,
                  ratio_HeadL=sp_HL_ratio,ratio_PeduncleH=sp_PH_ratio,ratio_DorsalL=sp_DL_ratio,
                  ratio_AnalL=sp_AL_ratio,ratio_PectoralBL=sp_PBL_ratio,ratio_PreED=sp_PreOD_ratio,
                  ratio_Pre_AL=sp_Pre_AL_ratio,ratio_Pre_PectoralBL=sp_PrePBL_ratio,ratio_Pre_DorsalL=sp_PreDL_ratio)

database <- data.frame(Site2=database[,"Site"],genus_species=database[,"genus_species"],trait_db)
database[,3:ncol(database)] <-  apply(database[,3:ncol(database)],2,as.numeric)


BigSp_list <- list("landmark"=All_tps,"data"=database)
saveRDS(BigSp_list,"results/intra_sp_morpho_traits.RDS") 
######################## end of script #########################################