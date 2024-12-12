################################################################
####
####
####
##################################################################

### Libary loading
lib <- c("gtools","ade4","cluster","phytools","caper","lavaan","FactoMineR",
         "factoextra","piecewiseSEM")
sapply(lib,library, character.only=TRUE)

### Functions loading
get_res_funct <- function(Data){
  MAtgo <- quasieuclid(daisy(x=data.frame(Data),metric=c("gower")))
  PCOA <- dudi.pco(MAtgo,nf=10,scannf=FALSE)
  c(CorPCOA1_Beta=cor(PCOA$li[,1],Data_G[,"Jst"]),CorPCOA1_Gamma=cor(PCOA$li[,1],Data_G[,"Jt"]))
} #end of get_res_funct

### Data loading and data modification
Data_G <- read.csv2("data/Tab_Abd_TRAITS_WIO_18sp.csv",sep=",",dec=",",row.names=2)[,-c(2:19)]
Data_G <-  Data_G[-15,c(5:7,21:48)]

### Genetic data loading
F_data <- readRDS("data/F_dataset_2.RDS")

Data_G <- data.frame(F_data[rownames(Data_G),],Data_G)

Data_G$Diet <- ordered(Data_G$Diet,levels= c("PK","IM","IS","OM","FC"))
Data_G$Body_size <- as.vector(Data_G$Body_size)
Data_G$PLD <- as.vector(Data_G$PLD)
Data_G$Abundance <- as.vector(log(Data_G$Abundance))

Data_PCOA <-Data_G[,c("Body_size","PLD","Home_range","Reproduction","Schooling","Diet","Abundance")]

###############################################################
### Relation between body size with JT and PLD with Jst     ###
###############################################################

### INtervariability loading 
Data_inter <- readRDS("results/Inter_specific_dataset.RDS")$inter_variability_obs


### Searching for the best PCOA that represent the correlation between beta and gamma
listcom <-lapply(1:7,function(i){combinations(7,i,1:7)})

Ressquar <- lapply(1: length(listcom),function(i){
  do.call(rbind,lapply(1:nrow(listcom[[i]]),function(y){
    Res <- get_res_funct(Data_PCOA[,listcom[[i]][y,]])
    Nm <- colnames(Data_PCOA)[listcom[[i]][y,]]
    c(Nm,Res)
  }))
})

### The best PCOA is for two traits  Ressquar[[2]][1,]

########################################################################
### Perform the phyloPCA; library("FactoMineR")
Data_PCA <- Data_PCOA[,c(1:2)]; rownames(Data_PCA) <- Sp2Prune 
PCA <- PCA(Data_PCA)

tiff("Figures/Extended_figure3.tiff",width = 16 ,height = 14 ,units = "cm",res=300,compression = "lzw")

rownames(Data_PCA) <- gsub("_"," ",rownames(Data_PCA))
res.pca <- prcomp(Data_PCA, scale = TRUE)

factoextra::fviz_pca_biplot(res.pca,labelsize = 3,repel = TRUE,alpha.var ="contrib")+
  theme(panel.border = element_rect(color = "black",size=0.4,fill=NA))
                                                                                               
dev.off()

################################################################################
###                                                                          ###
###     Relation between the Gamma/beta and the main traits by accounting    ###
###                            for the phylogeny                             ###
###                                                                          ###
################################################################################

### Get the phylogenetic tree from the rabosky tree
### Tree Rabosky loading and modification
Rabos_Tree <- read.tree("data/phylo/actinopt_12k_treePL.tre")
Rabos_Tree$tip.label[which(Rabos_Tree$tip.label=="Gomphosus_varius")] <- "Gomphosus_caeruleus"
Sp2Prune <- gsub(" ","_",as.character(Data_G$Taxon))
Prun_Rabos_T <- drop.tip(Rabos_Tree,setdiff(Rabos_Tree$tip.label,Sp2Prune)) 

Test <- Prun_Rabos_T 
Test$tip.label <- as.character(Data_G$Family)
Fam_tree  <- drop.tip(Test, which(duplicated(Test$tip.label)==T))

### Data preparation for the pgls
PGLS_dataIn <- data.frame(Species=rownames(PCA$ind$coord),Beta_D=Data_G[,"Jst"],
                          Gamma_D=Data_G[,"Jt"],PCA_1=PCA$ind$coord[,1],
                          PCA_2=PCA$ind$coord[,2],BS=Data_G[,"Body_size"],
                          Pld=Data_G[,"PLD"])

LBeta <- phylosig(tree=Prun_Rabos_T,method="lambda",x=setNames(PGLS_dataIn$Beta_D,PGLS_dataIn[,1]))
LGama <- phylosig(tree=Prun_Rabos_T,method="lambda",x=setNames(PGLS_dataIn$Gamma_D,PGLS_dataIn[,1]))
PGLS_data <- comparative.data(Prun_Rabos_T,PGLS_dataIn,names.col=Species)

### PGLS traits vs Gamma and beta
summary(lm(Gamma_D~BS,PGLS_dataIn))
summary(pgls(Gamma_D~BS,PGLS_data,lambda=LGama$lambda))

summary(lm(Beta_D~Pld,PGLS_dataIn))
summary(pgls(Beta_D~Pld,PGLS_data,lambda=LBeta$lambda))

### PGLS with PPCA-1
summary(lm(Gamma_D~PCA_1,PGLS_dataIn))
summary(lm(Beta_D~PCA_1,PGLS_dataIn))

summary(lm(Beta_D~Pld,PGLS_dataIn))
summary(pgls(Beta_D~PCA_1,PGLS_data,lambda=LBeta$lambda))

### PGLS with PPCA-2
summary(lm(Gamma_D~PCA_2,PGLS_dataIn))
summary(lm(Beta_D~PCA_2,PGLS_dataIn))

### Contribution of  BS and PLD to the first PCA axis
cor.test(PCA$ind$coord[,1],Data_G[,"Body_size"])
cor.test(PCA$ind$coord[,1],Data_G[,"PLD"])

##########################################################################
### SEM Intra
Data_PCoA_TV <- Data_G[,grep("ratio",colnames(Data_G))]
MAtgoTV <- quasieuclid(daisy(x=data.frame(Data_PCoA_TV),metric=c("gower")))
PCA_intra  <- dudi.pco(MAtgoTV,nf=2,scannf=FALSE)$li
#PCA_intra <- phyl.pca(tree=Prun_Rabos_T,Y=Data_PCoA_TV,method="BM",mode ="cov")

Data_inter <- Data_inter[-c(2,7),grep("ratio",colnames(Data_inter))]
MAtgoTV_inter <- quasieuclid(daisy(x=data.frame(Data_inter),metric=c("gower")))
PCA_inter <- dudi.pco(MAtgoTV_inter,nf=4,scannf=FALSE)$li
#PCA_inter <- phyl.pca(tree=Fam_tree,Y=Data_inter,method="BM",mode ="cov")$S
PCA_inter <- PCA_inter[as.character(Data_G$Family),]

SEM_dataIn <- data.frame(Species=Data_G$Taxon,Alpha=Data_G[,"Js"],Beta=Data_G[,"Jst"],Gamma=Data_G[,"Jt"],
                         BS=Data_G[,"Body_size"],PLD=Data_G[,"PLD"],Abd=Data_G[,"Abundance"],
                         PCA1=PCA$ind$coord[,1],PCA2=PCA$ind$coord[,2],
                         PCoA_inter=PCA_inter[,1],PCoA_interII=PCA_inter[,2],
                         Data_G[,grep("ratio",colnames(Data_G))],Nb_Sp=Data_G$Nb_sp_per_family,
                         PCoA_intra= PCA_intra[,1],PCoA_intraII=PCA_intra[,2])

### SEM Intra with synthetic axes
SEM_dataIn[,2:27] <- scale(SEM_dataIn[,2:27]) 
model <- psem(lm(PCA1~Beta+Gamma,data=SEM_dataIn),
              lm(Beta  ~ PCoA_inter+PCoA_intra,data=SEM_dataIn),
              lm(Gamma ~ PCoA_inter+PCoA_intra,data=SEM_dataIn),
              lm(PCoA_intra ~ PCoA_inter+PCA1,data=SEM_dataIn),
              lm(PCoA_inter ~ PCA1,data=SEM_dataIn)
)
plot(model,node_attrs = list(cex=0.2,cex.axis=0.3))


################################################################################ 
######################## end of script #########################################
################################################################################