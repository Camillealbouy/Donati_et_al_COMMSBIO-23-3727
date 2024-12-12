################################################################################
# SES between observed and random intra-sp morpholigical trait var
#
################################################################################

### Libary loading
lib <- c("spatstat.utils","cowplot","dplyr","vioplot","png")
sapply(lib,library, character.only=TRUE)

### Function loading
source("scripts/functions/SES_functions.R")

################################################################################
###                    Intraspecific data set analyses                       ###
################################################################################
### Data loading
#data <- readRDS("../data/micro/micro_data3.RDS") 
data <- readRDS("results/Intra_specific_dataset.RDS")

### Observed data preparation
observed_data <- data$intra_variability_pop
observed_data$genus_species <- rownames(observed_data)
rownames(observed_data) <- NULL
observed_data$category <- "obs"

### Random data preparation
random_data <- data$intra_variability_pop_random

all_random <- do.call(rbind,lapply(random_data,function(x){
  x$genus_species <-  rownames(x)
  rownames(x) <- NULL; x$category <- "null"
  x
}))

### Final data set observed + random data
ensemble <- rbind(observed_data, all_random[,colnames(observed_data)])

################################################################################
# plotting histograms with obs and null per species and saving SES tables
Dat2plot_exp <- split(ensemble, ensemble$genus_species)

Res_intra <- list()
for(j in 1:length(Dat2plot_exp)) {
  sp <- names(Dat2plot_exp)[j]
  data <- Dat2plot_exp[[j]]
  n <- unique(Dat2plot_exp[[j]]$n)
  col_n <- colnames(data)[grep("ratio",colnames(data))]
  
  ###  Indicators calculation for the traits distribution
  SES_df  <- data.frame(do.call(rbind,lapply(1:length(col_n), function(x){
    Mean_obs <- mean(data[,col_n[x]][which(data$category=="obs")])
    Mean_null <- mean(data[,col_n[x]][which(data$category=="null")])
    SES <- (Mean_obs-Mean_null)/sd(data[,col_n[x]][which(data$category=="null")])
    Q <- quantile(data[,col_n[x]][which(data$category=="null")],c(0.025,0.97)) 
    pval <- 2*pnorm(-abs(SES))
    return(c(SES=SES,Q,Pval=pval))
  })))
  
  SES_df$traits <- col_n
  SES_df$genus_species <- sp
  SES_df$n <- n
  colnames(SES_df) <- c("SES","q2.5","q97.5","Pval","trait","species","n")
  Res_intra[[j]] <- SES_df
} # end of j

Res_intra_SES <- do.call(rbind,Res_intra)

########## mean SES per trait ##################################################
trait_split <- split(Res_intra_SES , Res_intra_SES$trait)

Res <- do.call(rbind, lapply(trait_split, function(x){
  median_SES <- median(x$SES) 
  mean_SES <-  mean(x$SES) 
  c(median_SES=median_SES, mean_SES=mean_SES)
}))

### Data preparation for the Figure
yop <- ensemble
yop$genus_species1<- yop$genus_species

four_sp <- c("Ctenochaetus_striatus","Chaetodon_trifasciatus","Hemigymnus_fasciatus",
             "Chromis_atripectoralis","Pseudanthias_squamipinnis")

yep <- yop[which(yop$genus_species%in%four_sp),]

yep$genus_species <- gsub("Ctenochaetus_striatus","Acanthuridae",yep$genus_species)
yep$genus_species <- gsub("Chaetodon_trifasciatus","Chaetodontidae",yep$genus_species)
yep$genus_species <- gsub("Hemigymnus_fasciatus","Labridae",yep$genus_species)
yep$genus_species <- gsub("Chromis_atripectoralis","Pomacentridae",yep$genus_species)
yep$genus_species <- gsub("Pseudanthias_squamipinnis","Serranidae",yep$genus_species)

yep[grep("ratio",colnames(yep))] <- lapply(yep[grep("ratio",colnames(yep))], range01)

Dat2pLot  <- split(yep, f=yep$genus_species1)

Intra_elo <- lapply(four_sp,function(x) get_violplot(sp=x,ratio="ratio_ELO",data=Dat2pLot))
Intra_PeduncleH <-lapply(four_sp,function(x) get_violplot(sp=x,ratio="ratio_PeduncleH",data=Dat2pLot))
Intra_Jaw <- lapply(four_sp,function(x) get_violplot(sp=x,ratio="ratio_JawL",data=Dat2pLot))


################################################################################
####                 Interspecific data set analyses                        ####
################################################################################

### Data loading
#Data_Inter <- readRDS("Data/macro/macro_data3.RDS") 
Data_Inter <- readRDS("results/Inter_specific_dataset.RDS")

### Observed data preparation
observed_data <- Data_Inter$inter_variability_obs
observed_data$category <- "obs"
observed_data$Family <- rownames(observed_data)
observed_data <- observed_data[,-which(colnames(observed_data)=="nb_genus")]
rownames(observed_data) <- NULL

### Random data preparation
random_data <- Data_Inter$inter_variability_random

all_random_inter <- do.call(rbind,lapply(random_data,function (x){
  x$Family <-  rownames(x)
  rownames(x) <- NULL; x$category <- "null"
  x
}))

### Final data set observed + random data
ensemble_inter <- rbind(observed_data, all_random_inter[, colnames(observed_data)])

################################################################################
Dat2plot_inter <- split(ensemble_inter, ensemble_inter$Family)

Res_inter <- list()
for(j in 1:length(Dat2plot_inter)) {
  Fam <- names(Dat2plot_inter)[j]
  data <- Dat2plot_inter[[j]]
  n <- unique(Dat2plot_inter[[j]]$n)
  col_n <- colnames(data)[grep("ratio",colnames(data))]
  
  ###  Indicators calculation for the traits distribution
  SES_df  <- data.frame(do.call(rbind,lapply(1:length(col_n), function(x){
    Mean_obs <- mean(data[,col_n[x]][which(data$category=="obs")])
    Mean_null <- mean(data[,col_n[x]][which(data$category=="null")])
    SES <- (Mean_obs-Mean_null)/sd(data[,col_n[x]][which(data$category=="null")])
    Q <- quantile(data[,col_n[x]][which(data$category=="null")],c(0.025,0.97)) 
    pval <- 2*pnorm(-abs(SES))
    return(c(SES=SES,Q,Pval=pval))
  })))
  
  SES_df$traits <- col_n
  SES_df$genus_species <- Fam
  SES_df$n <- n
  colnames(SES_df) <- c("SES","q2.5","q97.5","Pval","trait","Family","n")
  Res_inter[[j]] <- SES_df
} # end of j

Res_inter_SES <- do.call(rbind,Res_inter)

### Data preparation for the Figure
Selected_Fam <- c("Acanthuridae","Chaetodontidae","Labridae","Pomacentridae","Serranidae")

yep_inter <- ensemble_inter[which(is.element(ensemble_inter$Family,Selected_Fam)),]
yep_inter[grep("ratio",colnames(yep_inter))] <- lapply(yep_inter[grep("ratio",colnames(yep_inter))], range01)

Dat2pLot  <- split(yep_inter, f=yep_inter$Family)
Five_fam <- names(Dat2pLot)

Inter_elo <- lapply(Five_fam,function(x) get_violplotInter(Fam=x,ratio="ratio_ELO",data=Dat2pLot))
Inter_PeduncleH <-lapply(Five_fam,function(x) get_violplotInter(Fam=x,ratio="ratio_PeduncleH",data=Dat2pLot))
Inter_Ana <- lapply(Five_fam,function(x) get_violplotInter(Fam=x,ratio="ratio_AnalL",data=Dat2pLot))


####################### FIGURE #################################################
Selc_sp <- c("Ctenochaetus_striatus","Chaetodon_trifasciatus","Hemigymnus_fasciatus",
             "Chromis_atripectoralis","Pseudanthias_squamipinnis","Chromis_ternatensis",
             "Chromis_weberi")
Selc_family <- c("Acanthuridae","Chaetodontidae","Labridae","Pomacentridae","Serranidae")

### Small picture loading 
Pict_Fnm <- list.files("figures/Pictures/Red/")
PictF <- lapply(Pict_Fnm, function(i){readPNG(paste("figures/Pictures/Red/",i,sep=""))})
names(PictF) <- Pict_Fnm

Pict_nm <- list.files("figures/Pictures/Blue/")
Pict <- lapply(Pict_nm, function(i){readPNG(paste("figures/Pictures/Blue/",i,sep=""))})
names(Pict) <- Pict_nm 

#############
tiff(filename="Figures/Figures_3_VF.tiff",width=17,height=20,units="cm",compression="lzw",res=300)

par(mfcol=c(3,2),mar=c(2,2.4,1.2,0.5),lwd=0.5)

### Inter traits variation plot
#---- Body elongation 
plot(rnorm(100),type="n",xlim=c(0,5),ylim=c(0,1),axes=F,ylab="Body elongation (SL/BD)",xlab="",mgp=c(1.35,0,0),cex.lab=0.9)
axis(side=2,las=2,cex.axis=0.6,mgp=c(2,0.5,0),tcl=-0.25)
axis(side=1, labels=c("","","","",""),cex.axis=0.6,mgp=c(2,0.15,0),tcl=-0.25,at=seq(0.5,4.5,1))
wex=0.5;lwd=0.35;cex=0.2

for (i in 1:length(Inter_elo)){
  vioplot(Inter_elo[[i]]$Values,col="black",add=T,side="left",horizontal=F,at=i-0.5,drawRect=F,wex=wex,lwd=lwd,cex=cex)
  segments(x0=i-0.5,x1=i-1,y0=Inter_elo[[i]]$obs,y1=Inter_elo[[i]]$obs,lty=2,col="orange")
}
box(lwd=1)
mtext("a",side=2,at=1.07,las=2,cex=1.2,line = 1.2)

#---- Caudal peduncle
par(mar=c(2.5,2.4,0.5,0.5))
plot(rnorm(100),type="n",xlim=c(0,5),ylim=c(0,1),axes=F,xlab="",ylab="Caudal peduncle h. (PH/BD)",mgp=c(1.35,0,0),cex.lab=0.8)
axis(side=2,las=2,cex.axis=0.6,mgp=c(2,0.5,0),tcl=-0.25)
axis(side=1, labels=c("","","","",""),cex.axis=0.6,mgp=c(2,0.15,0),tcl=-0.25,at=seq(0.5,4.5,1))

for (i in 1:length(Inter_PeduncleH)){
  vioplot(Inter_PeduncleH[[i]]$Values,col="black",add=T,side="left",horizontal=F,at=i-0.5,drawRect=F,wex=wex,lwd=lwd,cex=cex)
  segments(x0=i-0.5,x1=i-1,y0=Inter_PeduncleH[[i]]$obs,y1=Inter_PeduncleH[[i]]$obs,lty=2,col="orange")
}
box(lwd=1)

#---- Ana Fin length
par(mar=c(3,2.4,0,0.5))
plot(rnorm(100),type="n",xlim=c(0,5),ylim=c(0,1),axes=F,ylab="Anal fin length (AL/SL)",xlab="",mgp=c(1.35,0,0),cex.lab=0.9)
axis(side=2,las=2,cex.axis=0.6,mgp=c(2,0.5,0),tcl=-0.25)
axis(side=1, labels=Five_fam,at=seq(0.5,4.5,1),cex.axis=0.6,mgp=c(2,0.15,0),tcl=-0.25)

for (i in 1:length(Intra_Jaw)){
  vioplot(Inter_Ana[[i]]$Values,col="black",add=T,side="left",horizontal=F,at=i-0.5,drawRect=F,wex=wex,lwd=lwd,cex=cex)
  segments(x0=i-0.5,x1=i-1,y0=Inter_Ana[[i]]$obs,y1=Inter_Ana[[i]]$obs,lty=2,col="orange")
}
box(lwd=1)

### Intra traits variation plot
#--------------------
par(mar=c(2,2.4,1.2,0.5))
plot(rnorm(100),type="n",xlim=c(0,5),ylim=c(0,1),axes=F,ylab="Body elongation (SL/BD)",xlab="",mgp=c(1.35,0,0),cex.lab=0.9)
axis(side=1, labels=c("","","","",""),cex.axis=0.6,mgp=c(2,0.15,0),tcl=-0.25,at=seq(0.5,4.5,1))
axis(side=2,las=2,cex.axis=0.6,mgp=c(2,0.5,0),tcl=-0.25)
wex=0.5;lwd=0.35;cex=0.2

for (i in 1:length(Intra_elo)){
  vioplot(Intra_elo[[i]]$Values,col="black",add=T,side="left",horizontal=F,at=i-0.5,drawRect=F,wex=wex,lwd=lwd,cex=cex)
  segments(x0=i-0.5,x1=i-1,y0=Intra_elo[[i]]$obs,y1=Intra_elo[[i]]$obs,lty=2,col="#3B77A2")
}
box(lwd=1)
mtext("b",side=2,at=1.07,las=2,cex=1.2,line = 1.2)
#--------------------
par(mar=c(2.5,2.4,0.5,0.5))
plot(rnorm(100),type="n",xlim=c(0,5),ylim=c(0,1),axes=F,xlab="",ylab="Caudal peduncle h. (PH/BD)",mgp=c(1.35,0,0),cex.lab=0.9)
axis(side=2,las=2,cex.axis=0.6,mgp=c(2,0.5,0),tcl=-0.25)
axis(side=1, labels=c("","","","",""),cex.axis=0.6,mgp=c(2,0.15,0),tcl=-0.25,at=seq(0.5,4.5,1))

for (i in 1:length(Intra_PeduncleH)){
  vioplot(Intra_PeduncleH[[i]]$Values,col="black",add=T,side="left",horizontal=F,at=i-0.5,drawRect=F,wex=wex,lwd=lwd,cex=cex)
  segments(x0=i-0.5,x1=i-1,y0=Intra_PeduncleH[[i]]$obs,y1=Intra_PeduncleH[[i]]$obs,lty=2,col="#3B77A2")
}
box(lwd=1)

#--------------------
par(mar=c(3,2.4,0,0.5))
plot(rnorm(100),type="n",xlim=c(0,5),ylim=c(0,1),axes=F,ylab="Jaw length (JH/HL)",xlab="",mgp=c(1.35,0,0),cex.lab=0.9)
axis(side=2,las=2,cex.axis=0.6,mgp=c(2,0.5,0),tcl=-0.25)
axis(side=1, labels = c("C.Striatus","C. trifasciatus", "H. fasciatus","C. atripectoralis", "P. squamipinnis"),
     at=seq(0.5,4.5,1),cex.axis=0.6,mgp=c(2,0.15,0),tcl=-0.25)

for (i in 1:length(Intra_Jaw)){
  vioplot(Intra_Jaw[[i]]$Values,col="black",add=T,side="left",horizontal=F,at=i-0.5,drawRect=F,wex=wex,lwd=lwd,cex=cex)
  segments(x0=i-0.5,x1=i-1,y0=Intra_Jaw[[i]]$obs,y1=Intra_Jaw[[i]]$obs,lty=2,col="#3B77A2")
}
box(lwd=1)
dev.off()

#################################################################################
### Data preparation for the supMat Figure
Genet_Data <- readRDS("data/F_dataset_2.RDS")
rownames(Genet_Data) <- Genet_Data$Taxon


yop <- ensemble
yop$genus_species1<- yop$genus_species

other_sp <-  c("Caranx_melampygus","Chaetodon_trifasciatus","Chromis_atripectoralis","Chromis_ternatensis", "Chromis_weberi",
                "Ctenochaetus_striatus","Dascyllus_aruanus","Dascyllus_carneus","Dascyllus_trimaculatus","Gomphosus_caeruleus",
                "Halichoeres_hortulanus","Hemigymnus_fasciatus","Lutjanus_kasmira","Myripristis_violacea","Oxymonacanthus_longirostris"
                          ,"Parupeneus_macronemus","Pseudanthias_squamipinnis")
ab_sp <-  c("Car_mel","Cha_tri","Chr_atr","Chro_ter", "Chr_web","Cte_str","Das_aru","Das_car","Das_tri","Gom_cae",
            "Hal_hor","Hem_fas","Lut_kas","Myr_vio","Oxy_lon","Par_mac","Pse_squ")

yep <- yop[which(yop$genus_species%in%other_sp),]
yep$Family <-sapply(1:length(yep$genus_species), function(i) Genet_Data[yep$genus_species[i],"Family"])
yep[grep("ratio",colnames(yep))] <- lapply(yep[grep("ratio",colnames(yep))], range01)

Dat2pLot  <- split(yep, f=yep$genus_species1)
ratio <- colnames(yep)[grep("ratio",colnames(yep))]

Intra_res <- lapply(1:length(ratio),function(i){
  lapply(other_sp,function(x) get_violplot(sp=x,ratio=ratio[i],data=Dat2pLot))
})
   
names <- c("Body elongation (SL/BD)","Caudal peduncle elongation (PL/PH)","Eye-head size relation (ED/HL)","Jaw-head length relation(JL/HL)",
           "Head bogy length relation (HL/SL)", "Caudal peduncle-body heigth relation (PH/BD)","Stand Dorsal fin length (DL/SL)",
           "StandAnal fin length (AL/SL) ","Stand pectoral fin length (Pec-L/SL)", "Eye-head position relation (Pre-ED/HL)",
           "Stand Pre-Anal length (Pre-ED/HL)","Stand Pre-pectoral length (Pre-Pec-L/SL)",
           "Stand Pre-Dorsal length (Pre-DL/SL)")

wex=0.5;lwd=0.35;cex=0.2


tiff(filename="Figures/Figures_SM_1_2.tiff",width=24,height=12,units="cm",compression="lzw",res=300)
par(mfrow=c(2,4),mar=c(2,2.4,1.5,0.5),lwd=0.5)

### Inter traits variation plot
for (j in 1:8){
  if(j==5 | j==6 | j==7 | j==8){par(mar=c(3,2.4,1,0.5))}
  #if(j==9 | j==10 | j==11 | j==12){par(mar=c(3,2.4,0.5,0.5))}
  
  Intra <- Intra_res[[j]]
  plot(rnorm(100),type="n",xlim=c(0,17),ylim=c(0,1),axes=F,ylab=names[j],xlab="",mgp=c(1.45,0,0),cex.lab=0.8)
  axis(side=2,las=2,cex.axis=0.75,mgp=c(2,0.5,0),tcl=-0.25)
  axis(side=1, labels=rep("",17),cex.axis=0.6,mgp=c(2,0.15,0),tcl=-0.25,at=seq(0.5,16.5,1),las=2)
  if(j==5|j==6|j==7|j==8){
    axis(side=1, labels=ab_sp,cex.axis=0.7,mgp=c(2,0.4,0),tcl=-0.25,at=seq(0.5,16.5,1),las=2)
  }
  
  for (i in 1:length(Intra)){
    vioplot(Intra[[i]]$Values,col="black",add=T,side="left",horizontal=F,at=i-0.5,drawRect=F,wex=wex,lwd=lwd,cex=cex)
    segments(x0=i-0.5,x1=i-1,y0=Intra[[i]]$obs,y1=Intra[[i]]$obs,lty=2,col="orange",lwd=0.9)
  }
  box(lwd=1)
  mtext(letters[j],side=2,at=1.09,las=2,cex=1,line=1.2)
}
dev.off()


tiff(filename="Figures/Figures_SM_2_2.tiff",width=24,height=12,units="cm",compression="lzw",res=300)
par(mfrow=c(2,4),mar=c(2,2.4,1.5,0.5),lwd=0.5)
### Inter traits variation plot
for (j in 9:13){
  #if(j==5 | j==6 | j==7 | j==8){par(mar=c(3,2.4,1,0.5))}
  if(j==13){par(mar=c(3,2.4,0.5,0.5))}
  
  Intra <- Intra_res[[j]]
  plot(rnorm(100),type="n",xlim=c(0,17),ylim=c(0,1),axes=F,ylab=names[j],xlab="",mgp=c(1.45,0,0),cex.lab=0.8)
  axis(side=2,las=2,cex.axis=0.75,mgp=c(2,0.5,0),tcl=-0.25)
  axis(side=1, labels=rep("",17),cex.axis=0.6,mgp=c(2,0.15,0),tcl=-0.25,at=seq(0.5,16.5,1),las=2)
  if(j==10|j==11|j==12|j==13){
    axis(side=1, labels=ab_sp,cex.axis=0.7,mgp=c(2,0.4,0),tcl=-0.25,at=seq(0.5,16.5,1),las=2)
  }
  
  for (i in 1:length(Intra)){
    vioplot(Intra[[i]]$Values,col="black",add=T,side="left",horizontal=F,at=i-0.5,drawRect=F,wex=wex,lwd=lwd,cex=cex)
    segments(x0=i-0.5,x1=i-1,y0=Intra[[i]]$obs,y1=Intra[[i]]$obs,lty=2,col="orange",lwd=0.9)
  }
  box(lwd=1)
  mtext(letters[j],side=2,at=1.09,las=2,cex=1,line=1.2)
}
dev.off()

################################################################################

Selected_Fam <- c("Acanthuridae","Carangidae","Chaetodontidae","Holocentridae","Labridae","Lutjanidae","Monacanthidae"
                                 ,"Mullidae","Pomacentridae","Serranidae")

yep_inter <- ensemble_inter[which(is.element(ensemble_inter$Family,Selected_Fam)),]
yep_inter[grep("ratio",colnames(yep_inter))] <- lapply(yep_inter[grep("ratio",colnames(yep_inter))], range01)

Dat2pLot  <- split(yep_inter, f=yep_inter$Family)
Five_fam <- names(Dat2pLot)

Inter_res <- lapply(1:length(ratio),function(i){
  lapply(Five_fam,function(x) get_violplotInter(Fam=x,ratio=ratio[i],data=Dat2pLot))
})


tiff(filename="Figures/Figures_SM_INter_1_2.tiff",width=24,height=12,units="cm",compression="lzw",res=300)
par(mfrow=c(2,4),mar=c(2,2.4,1.5,0.5),lwd=0.5)

### Inter traits variation plot
for (j in 1:8){
  if(j==5 | j==6 | j==7 | j==8){par(mar=c(3.7,2.4,0.3,0.5))}
 
  Inter<- Inter_res[[j]]
  plot(rnorm(100),type="n",xlim=c(0,10),ylim=c(0,1),axes=F,ylab=names[j],xlab="",mgp=c(1.35,0,0),cex.lab=0.8)
  axis(side=2,las=2,cex.axis=0.6,mgp=c(2,0.5,0),tcl=-0.25)
  axis(side=1, labels=rep("",10),cex.axis=0.6,mgp=c(2,0.15,0),tcl=-0.25,at=seq(0.5,9.5,1),las=2)
  if(j==5|j==6|j==7|j==8){
    axis(side=1, labels=names(Dat2pLot),cex.axis=0.55,mgp=c(2,0.4,0),tcl=-0.25,at=seq(0.5,9.5,1),las=2)
  }
  
  for (i in 1:length(Inter)){
    vioplot(Inter[[i]]$Values,col="black",add=T,side="left",horizontal=F,at=i-0.5,drawRect=F,wex=wex,lwd=lwd,cex=cex)
    segments(x0=i-0.5,x1=i-1,y0=Inter[[i]]$obs,y1=Inter[[i]]$obs,lty=2,col="Blue",lwd=0.9)
  }
  box(lwd=1)
  mtext(letters[j],side=2,at=1.09,las=2,cex=1,line=1.2)
}
dev.off()



tiff(filename="Figures/Figures_SM_INter_2_2.tiff",width=24,height=12,units="cm",compression="lzw",res=300)
par(mfrow=c(2,4),mar=c(2,2.4,1.5,0.5),lwd=0.5)

### Inter traits variation plot
for (j in 9:13){
  if(j==13){par(mar=c(3.7,2.4,0.3,0.5))}
  
  Inter<- Inter_res[[j]]
  plot(rnorm(100),type="n",xlim=c(0,10),ylim=c(0,1),axes=F,ylab=names[j],xlab="",mgp=c(1.35,0,0),cex.lab=0.8)
  axis(side=2,las=2,cex.axis=0.6,mgp=c(2,0.5,0),tcl=-0.25)
  axis(side=1, labels=rep("",10),cex.axis=0.6,mgp=c(2,0.15,0),tcl=-0.25,at=seq(0.5,9.5,1),las=2)
  if(j==10|j==11|j==12|j==13){
    axis(side=1, labels=names(Dat2pLot),cex.axis=0.55,mgp=c(2,0.4,0),tcl=-0.25,at=seq(0.5,9.5,1),las=2)
  }
  
  for (i in 1:length(Inter)){
    vioplot(Inter[[i]]$Values,col="black",add=T,side="left",horizontal=F,at=i-0.5,drawRect=F,wex=wex,lwd=lwd,cex=cex)
    segments(x0=i-0.5,x1=i-1,y0=Inter[[i]]$obs,y1=Inter[[i]]$obs,lty=2,col="Blue",lwd=0.9)
  }
  box(lwd=1)
  mtext(letters[j],side=2,at=1.09,las=2,cex=1,line=1.2)
}
dev.off()
