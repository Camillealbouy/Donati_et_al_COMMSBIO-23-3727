###############################################################################################
# Giulia Donati, giulia.donati@usys.ethz.ch
# October 2019
# Project: a continuum of macro- to micro- evolution of morphological traits in tropical reef fishes
# Script aim: inter-specific morphological data collection landmarks and morpho traits:
# devtools::install_github("davidcsterratt/retistruct@v0.6.4", subdir="pkg/retistruct")
###############################################################################################
rm(list=ls())

## Library loading
lib_vect <- c("geomorph","rgl","shapes","abind","scales"
              ,"dplyr","reshape2","tibble","retistruct")
sapply(lib_vect,library,character.only=TRUE)

### Function loading
source("scripts/functions/morphological_functions.R")

################################################################################
###            Part 1. Morphological data assemblage                        ####
################################################################################ 

### Load morpho data - to extract species names and ID
targets <-read.tps(data="data/tps.trimmed.TPS")
subset2 <- unique(targets[,c(3:4)])

subset2$IMAGE <- gsub(".jpg","", as.character(subset2$IMAGE)) 
subset2$IMAGE <- gsub(" ","_", as.character(subset2$IMAGE)) 
colnames(subset2)<-c("genus_species","sample_ID")

# Load morphometric data - landmarks
targets1 <- readland.tps("data/tps.trimmed.TPS", specID = "ID")
targets1 <- targets1[-c(12:14),,]#removing non matching landmarks, i.e. lm 12,13 and 14 From Claverie et al. 2014

# checking samples ID and species names
dimnames(targets1)[[3]] <- subset2$genus_species

# extract only the genus as some specimen are identified only up to the genus
subset2$Genus <- do.call(rbind,strsplit(as.character(subset2$genus_species),'_'))[,1]

### now subsetting according to the micro data - that is species of Families sampled in the CA and WIO

micro_family <- c("Acanthuridae","Pomacentridae","Carangidae","Chaetodontidae",
                  "Pomacentridae","Labridae","Lutjanidae","Lethrinidae",
                  "Mullidae","Holocentridae","Monacanthidae","Pinguipedidae",
                  "Serranidae","Caesionidae","Zanclidae")

# Getting all the genuses of the species of the Families of the CA, WIO
gaspar <- read.csv("data/edited_clean_gaspar_gd.csv",h=T,row.names=1)[,c(1,8,9)]
colnames(gaspar) <- c("Genus_and_species", "Family","Genus")
gaspar <- gaspar[which(gaspar$Family%in%micro_family),]# 1923 species out of the 6315 are in the Families sampled in the CA and WIO species.

#reduction of subset2 to only the genus of interest
subset3 <- dplyr::semi_join(subset2,gaspar, by = "Genus") # All rows in subset2 that have a match in gaspar2.

#reduction of targets 1, which contains the corresponding landmark info of subset3
targets2 <- targets1[,,which(dimnames(targets1)[[3]]%in%subset3$genus_species)]

# replace the sample ID by the genus name
dimnames(targets2)[[3]] <- subset3$genus_species

# match the corresponding family names to the genus of gaspar2
gaspar2 <- unique(gaspar[,-1])

# this is to add genus and family column to the reduced landmark dataset
yep <- dplyr::right_join(gaspar2, subset3, by = "Genus")# Join matching rows from gaspar2 to subset3
rownames(yep) <- yep$genus_species
yep <- yep[dimnames(targets2)[[3]],]

#replacinG Scaridae with Labridae
yep$Family <- gsub("Scaridae","Labridae", as.character(yep$Family))

#replace genus name by family name for plotting
dimnames(targets2)[[3]] <- yep$Genus

# Generalized procruste analysis on the landmark data
Y.gpa<-gpagen(targets2, print.progress = T) 
All_tps <- Y.gpa$coords


################################################################################
# collection of morphological traits
# creation of intersection and new points
# head
#eye diameter
################################################################################  

# creation of additional points (intersection points) that will serve to create different length measures
# for more detail and figure see script XXX
back <- do.call(rbind,lapply(1:dim(All_tps)[3],function(i){
  line.line.intersection(P1=c(max(All_tps[,1,i]), max(All_tps[,2,i])), 
                         P2=c(max(All_tps[,1,i]), min(All_tps[,2,i])),
                         P3 = c(min(All_tps[,1,i]),min(All_tps[,2,i])), 
                         P4 = c(min(All_tps[,1,i])+3,min(All_tps[,2,i])), interior.only = FALSE)})) 

front <-  do.call(rbind,lapply(1:dim(All_tps)[3],function(i){
  line.line.intersection(P1=c(min(All_tps[,1,i]),min(All_tps[,2,i])), 
                         P2=c(min(All_tps[,1,i]),max(All_tps[,2,i])), 
                         P3=c(back[i,1],back[i,2]),P4=c(back[i,1]-3,back[i,2]), 
                         interior.only = FALSE)}))

top <- do.call(rbind, lapply(1:dim(All_tps)[3],function(i){
  line.line.intersection(P1=c(min(All_tps[,1,i]),max(All_tps[,2,i])), 
                         P2=c(max(All_tps[,1,i]),max(All_tps[,2,i])), 
                         P3=c(back[i,1],back[i,2]),P4=c(back[i,1],back[i,2]+5),
                         interior.only = FALSE)}))

hl <- do.call(rbind, lapply(1:dim(All_tps)[3],function(i){
  line.line.intersection(P1=c(front[i,1],front[i,2]), 
                         P2=c(back[i,1],back[i,2]),P3=c(All_tps[11,1,i],All_tps[11,2,i]), 
                         P4=c(All_tps[11,1,i],All_tps[11,2,i]-3),interior.only = FALSE)}))

pOD <- do.call(rbind, lapply(1:dim(All_tps)[3],function(i){
  line.line.intersection(P1=c(front[i,1],front[i,2]),P2=c(back[i,1],back[i,2]), 
                         P3=c(All_tps[12,1,i],All_tps[12,2,i]), 
                         P4=c(All_tps[12,1,i],All_tps[12,2,i]-3),interior.only=FALSE)}))

pAL <- do.call(rbind, lapply(1:dim(All_tps)[3],function(i){
  line.line.intersection(P1=c(front[i,1],front[i,2]),P2=c(back[i,1],back[i,2]), 
                         P3=c(All_tps[8,1,i],All_tps[8,2,i]), 
                         P4=c(All_tps[8,1,i],All_tps[8,2,i]-3),interior.only = FALSE)}))


dl_start <- do.call(rbind, lapply(1:dim(All_tps)[3],function(i){
  line.line.intersection(P1=c(front[i,1],front[i,2]),P2=c(back[i,1],back[i,2]), 
                         P3=c(All_tps[3,1,i],All_tps[3,2,i]),P4=c(All_tps[3,1,i],All_tps[3,2,i]-5)
                         ,interior.only=FALSE)}))

dl_end <- do.call(rbind, lapply(1:dim(All_tps)[3],function(i){
  line.line.intersection(P1=c(front[i,1],front[i,2]),P2=c(back[i,1],back[i,2]), 
                         P3=c(All_tps[4,1,i],All_tps[4,2,i]),P4=c(All_tps[4,1,i],All_tps[4,2,i]-5), 
                         interior.only=FALSE)}))

af_end <-  do.call(rbind, lapply(1:dim(All_tps)[3],function(i){
  line.line.intersection(P1=c(front[i,1],front[i,2]),P2=c(back[i,1],back[i,2]), 
                         P3=c(All_tps[7,1,i],All_tps[7,2,i]),P4=c(All_tps[7,1,i],All_tps[7,2,i]-5), 
                         interior.only = FALSE)}))

ped_top <- do.call(rbind, lapply(1:dim(All_tps)[3],function(i){
  line.line.intersection(P1=c(top[i,1],top[i,2]),P2=c(back[i,1],back[i,2]), 
                         P3=c(All_tps[4,1,i],All_tps[4,2,i]),P4=c(All_tps[4,1,i]+5,All_tps[4,2,i]),
                         interior.only = FALSE)}))


ped_bottom <- do.call(rbind, lapply(1:dim(All_tps)[3],function(i){
  line.line.intersection(P1=c(top[i,1],top[i,2]),P2=c(back[i,1],back[i,2]), 
                         P3=c(All_tps[7,1,i],All_tps[7,2,i]),P4=c(All_tps[7,1,i]+5,All_tps[7,2,i]),
                         interior.only = FALSE)}))

prePBL <- do.call(rbind, lapply(1:dim(All_tps)[3],function(i){
  line.line.intersection(P1=c(front[i,1],front[i,2]),P2=c(back[i,1],back[i,2]), 
                         P3=c(All_tps[9,1,i],All_tps[9,2,i]),P4=c(All_tps[9,1,i],All_tps[9,2,i]-5),
                         interior.only = FALSE)}))


# euclidean distances between landmarks and/or intersection points
sp_HL <- sapply(1:dim(All_tps)[3],function(i){ dist(rbind(front[i,], hl[i,]),method="euclidean")})
# head
sp_ED <-  sapply(1:dim(All_tps)[3],function(i){dist(rbind(All_tps[12,,i],All_tps[13,,i]),method="euclidean")})
sp_PreOD <- sapply(1:dim(All_tps)[3],function(i){dist(rbind(front[i,],pOD[i,]), method="euclidean")})
sp_ED_ratio <- sp_ED/sp_HL #eye diameter/head length
sp_PreOD_ratio <- sp_PreOD/sp_HL #pre_eye/HL 
sp_JL <-  sapply(1:dim(All_tps)[3],function(i){dist(rbind(All_tps[2,,i],All_tps[1,,i]),method="euclidean")})
sp_JL_ratio <- sp_JL/sp_HL #jaw length/head length

# trunk
sp_SL <- sapply(1:dim(All_tps)[3],function(i){dist(rbind(front[i,],back[i,]),method="euclidean")})
sp_BH <- sapply(1:dim(All_tps)[3],function(i){dist(rbind(top[i,],back[i,]),method="euclidean")})
sp_HL_ratio <- sp_HL/sp_SL #head lenght/standard length
sp_ELO_ratio <- sp_SL/sp_BH

# tail
sp_PH <- sapply(1:dim(All_tps)[3],function(i){dist(rbind(ped_top[i,],ped_bottom[i,]), method = "euclidean")})
sp_PL <- sapply(1:dim(All_tps)[3],function(i){dist(rbind(dl_end[i,],back[i,]),method="euclidean")})
sp_PH_ratio <- sp_PH/sp_BH #peduncle heigh/ body height
sp_PELO_ratio <- sp_PL/sp_PH

# fins
sp_DL <- sapply(1:dim(All_tps)[3],function(i){dist(rbind(dl_start[i,],dl_end[i,]),method = "euclidean")})
sp_PreDL <- sapply(1:dim(All_tps)[3],function(i){dist(rbind(front[i,],dl_start[i,]),method = "euclidean")})
sp_PreDL_ratio <- sp_PreDL/sp_SL
sp_DL_ratio <- sp_DL/sp_SL #dorsal fin length / SL


sp_AL <- sapply(1:dim(All_tps)[3],function(i){dist(rbind(pAL[i,],af_end[i,]),method = "euclidean")})
sp_AL_ratio <- sp_AL/sp_SL

sp_Pre_AL <- sapply(1:dim(All_tps)[3],function(i){dist(rbind(front[i,],pAL[i,]),method="euclidean")})
sp_Pre_AL_ratio <- sp_Pre_AL/sp_SL

sp_PBL <- sapply(1:dim(All_tps)[3],function(i){dist(rbind(All_tps[9,,i],All_tps[10,,i]),method = "euclidean")})
sp_PBL_ratio <- sp_PBL/sp_SL

sp_PrePBL <- sapply(1:dim(All_tps)[3],function(i){dist(rbind(front[i,],prePBL[i,]),method = "euclidean")})
sp_PrePBL_ratio <- sp_PrePBL/sp_SL


database <- cbind(yep,EyeD=sp_ED,JawL=sp_JL,HeadL=sp_HL,PreED=sp_PreOD,PeduncleH=sp_PH,
                  PeduncleL=sp_PL,DorsalL=sp_DL,AnalL=sp_AL,PectoralBL=sp_PBL,Pre_Anal=sp_Pre_AL,
                  SL=sp_SL,BH=sp_BH,Pre_PectoralBL=sp_PrePBL,Pre_DorsalL=sp_PreDL,ratio_ELO=sp_ELO_ratio
                  ,ratio_PeduncleELO=sp_PELO_ratio,ratio_EyeD=sp_ED_ratio,ratio_JawL=sp_JL_ratio,
                  ratio_HeadL=sp_HL_ratio,ratio_PeduncleH=sp_PH_ratio,ratio_DorsalL=sp_DL_ratio,
                  ratio_AnalL=sp_AL_ratio,ratio_PectoralBL=sp_PBL_ratio,ratio_PreED=sp_PreOD_ratio,
                  ratio_Pre_AL=sp_Pre_AL_ratio,ratio_Pre_PectoralBL=sp_PrePBL_ratio,ratio_Pre_DorsalL=sp_PreDL_ratio)
rownames(database) <- database$genus_species

#adding family information
database <- dplyr::arrange(database, as.numeric(database$sample_ID))
TC_lm_data <- list(landmark=All_tps, data=database)


saveRDS(TC_lm_data, file="results/TC_dataset_landmarks_metadata.RDS")
############################################################### end of R script. 

