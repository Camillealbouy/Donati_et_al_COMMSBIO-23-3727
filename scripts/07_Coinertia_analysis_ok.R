################################################################################
###     Continuity in morphological disparity in tropical reef fishes        ###
###                       across evolutionary scales                         ###
### Dr. Giulia Donati, giulia.donati@usys.ethz.ch &                          ###
### Dr Camille Albouy calbouy@ethz.ch                                        ###
###                                                                          ###
###      Coinertia analysis between the inter and intra specific data sets   ###
################################################################################

### Library loading
library("ade4")

### Selected species
Selc_sp <- c("Caranx_melampygus","Chaetodon_trifasciatus","Chromis_atripectoralis",
             "Chromis_ternatensis", "Chromis_weberi","Ctenochaetus_striatus",
             "Dascyllus_aruanus","Dascyllus_carneus","Dascyllus_trimaculatus",
             "Gomphosus_caeruleus","Halichoeres_hortulanus","Hemigymnus_fasciatus",
             "Lutjanus_kasmira","Myripristis_violacea","Oxymonacanthus_longirostris"
             ,"Parupeneus_macronemus","Pseudanthias_squamipinnis")

Selc_family <- c("Acanthuridae","Carangidae","Chaetodontidae","Holocentridae",
                 "Labridae","Lutjanidae","Monacanthidae","Mullidae","Pomacentridae"
                 ,"Serranidae")

################################################################################
###                    Intra specific data set preparation                   ###
################################################################################
Intra <- readRDS("results/Intra_specific_dataset.RDS")
Intra <- Intra$intra_variability_pop

Family <- c("Pomacentridae","Carangidae","Chaetodontidae","Pomacentridae","Pomacentridae",
            "Pomacentridae","Acanthuridae","Pomacentridae","Pomacentridae","Pomacentridae",
            "Labridae","Labridae", "Labridae","Lutjanidae","Lethrinidae","Lethrinidae","Holocentridae",
            "Acanthuridae", "Labridae","Monacanthidae","Pinguipedidae",
            "Mullidae","Serranidae","Caesonidae","Zanclidae")

Intra <- cbind(Family,Intra)                
Intra <- Intra[-which(Intra$Family=="Zanclidae" | Intra$Family=="Pinguipedidae"),]

Intra <- Intra[order(Intra$Family),]
Intra <- cbind(nb_fam = c(1,2,1,1,1,1,1,2,3,4,1,2,1,1,1,1,2,3,4,5,6,7,1) ,Intra)
Intra$Intra_fam_name <- paste(Intra$Family,Intra$nb_fam, sep="_")
Intra <- Intra[order(Intra$Intra_fam_name),]
colnames(Intra) <- paste(colnames(Intra),"Intra",sep="_")

################################################################################
###                   Inter specific data set preparation                    ###
################################################################################

Inter <- readRDS("results/Inter_specific_dataset.RDS")$inter_variability_obs

Inter$Family <- rownames(Inter)
Inter <- Inter[order(Inter$Family),]
rownames(Inter) <- NULL

# repeat rows to gain same number as Intra data set
families <- c("Acanthuridae", rep("Labridae", 3), "Lethrinidae", 
              rep("Pomacentridae", 6))
to_add <- do.call(rbind, lapply(families, function(fam) Inter[Inter$Family == fam, ]))

Inter <- rbind(Inter,to_add)
Inter <- Inter[order(Inter$Family),]
nb_fam <- c(1,2,1,1,1,1,1,2,3,4,1,2,1,1,1,1,2,3,4,5,6,7,1) 

Inter <- cbind(Inter,nb_fam)
Inter$Inter_fam_name <- paste(Inter$Family,Inter$nb_fam, sep="_")
Inter <- Inter[order(Inter$Inter_fam_name),]
colnames(Inter) <- paste(colnames(Inter),"Inter",sep="_")

################################################################################
###                        Co-inertia analyses                               ###
################################################################################

### Data set compilation
CoI_dataset <- cbind(Intra,Inter)
CoI_dataset  <- CoI_dataset[Selc_sp,c(17:29,46:58)]

### Intraspecific PCA
pca_intra <- dudi.pca(CoI_dataset[,grep("Intra",colnames(CoI_dataset))],
                      scannf=F,nf=4)

### Interspecific PCA
pca_inter <- dudi.pca(CoI_dataset[,grep("Inter",colnames(CoI_dataset))],
                      scannf=F,nf=4)

### Co-inertia analysis
coin1 <- coinertia(pca_inter, pca_intra,nf=,scannf=F)
rv1 <- RV.rtest(pca_inter$tab, pca_intra$tab, 99)

### Eigen values calculation 
eigs_inter <- round(pca_inter$eig/sum(pca_inter$eig)*100,2)
eigs_intra <- round(pca_intra$eig/sum(pca_intra$eig)*100,2)

### Final results compilation
v_obs <- c(eigs_inter[1], eigs_inter[2],eigs_intra[1],eigs_intra[2],round(coin1$RV,3),pval=rv1$pvalue)
names(v_obs) <-c("inter_PCAxis1","inter_PCAxis2","intra_PCAxis1","intra_PCAxis2","RV_coeff","p.value")

################################################################################ 
######################## end of script #########################################
################################################################################ 
