################################################################################
###     Continuity in morphological disparity in tropical reef fishes        ###
###                       across evolutionary scales                         ###
### Dr. Giulia Donati, giulia.donati@usys.ethz.ch &                          ###
### Dr Camille Albouy calbouy@ethz.ch                                        ###
###                                                                          ###
################################################################################
###      Investigating association between intra and interspecific           ###
###          datasets for each traits without Carangidae                     ###
################################################################################

### Library
lib_vect <- c("ade4","ape","nlme","geiger","caper","phytools","phylobase","dplyr","png")
sapply(lib_vect,library,character.only=TRUE)

### Species and families of interest
Selc_sp <- c("Chaetodon_trifasciatus","Chromis_atripectoralis","Chromis_ternatensis",
             "Chromis_weberi","Ctenochaetus_striatus","Dascyllus_aruanus",
             "Dascyllus_carneus","Dascyllus_trimaculatus","Gomphosus_caeruleus",
             "Halichoeres_hortulanus","Hemigymnus_fasciatus","Lutjanus_kasmira",
             "Myripristis_violacea","Oxymonacanthus_longirostris","Parupeneus_Internemus"
             ,"Pseudanthias_squamipinnis")

Selc_family <- c("Acanthuridae","Chaetodontidae","Holocentridae","Labridae","Lutjanidae",
                 "Monacanthidae","Mullidae","Pomacentridae","Serranidae")

### function 2 load
source("Scripts/functions/Plot_functions.R")

################################################################################
### Correlating observed inter-specific variability values - temporal corr.  ###
################################################################################ 
### Data compilation at the species level
Intra <- readRDS("results/Intra_specific_dataset.RDS")
Intra <- Intra$intra_variability_pop[Selc_sp,]

Family <- c("Chaetodontidae","Pomacentridae","Pomacentridae","Pomacentridae",
            "Acanthuridae","Pomacentridae","Pomacentridae","Pomacentridae",
            "Labridae","Labridae","Labridae","Lutjanidae","Holocentridae",
            "Monacanthidae","Mullidae","Serranidae")

Intra <- cbind(Family=Family,nb_fam=as.vector(table(Family)[Family]),Intra)
colnames(Intra) <- paste(colnames(Intra),"Intra",sep="_")

### Data compilation at the family level
Inter <- readRDS("results/Inter_specific_dataset.RDS")$inter_variability_obs
Inter <- Inter[Selc_family,]
Inter <- cbind(Family=rownames(Inter),Inter)

# repeat rows to gain same number as Intra data set
Inter <- Inter[Intra$Family_Intra,]
Inter <- cbind(Inter,nb_fam=Intra$nb_fam_Intra)
colnames(Inter) <- paste(colnames(Inter),"Inter",sep="_")

temporal <- cbind(Intra,Inter)

################################################################################
###.                             PGLS                                        ###
################################################################################
### Rabosky et al. 2018 super-tree
Rabosky <- read.tree("Data/phylo/actinopt_12k_treePL.tre")
Rabosky$tip.label <- gsub("Gomphosus_varius","Gomphosus_caeruleus",as.character(Rabosky$tip.label))
Rabosky$tip.label <- gsub(Rabosky$tip.label[grep("Acanthurus",Rabosky$tip.label)][1],
                          "Acanthurus_tractus",as.character(Rabosky$tip.label))

Rabosky_toPrune <- rownames(temporal)
TreeRabosky_sp <- drop.tip(Rabosky,setdiff(Rabosky$tip.label,Rabosky_toPrune)) 

### Combine the data and the tree into an object
test <- temporal
test <- cbind(Species=rownames(test),test)
comp.data <- comparative.data(TreeRabosky_sp,test,names.col="Species",vcv.dim=2,warn.dropped=T)

################################################################################
###                      Plots: inter intra corr                             ###
################################################################################
ratio <- colnames(test)[grep("ratio_",colnames(test))]
Intra_traits <- ratio[grep("_Intra",as.character(ratio))]
Inter_traits <- gsub("_Intra","_Inter", as.character(Intra_traits))

mains_end <- c("Body elongation","Caudal peduncle elongation",
               "Eye size","Jaw length","Head length", "Caudal peduncle height",
               "Dorsal fin size","Anal fin size","Pectoral fin size","Eye position",
               "Anal fin insertion","Pectoral fin insertion","Dorsal fin insertion")

modelpgls <- list()
for(i in 1:length(Intra_traits)){
  modelpgls[[i]] <- pgls(test[,Intra_traits[i]]~test[,Inter_traits[i]],comp.data) 
  modelpgls[[i]] <- summary(modelpgls[[i]])  
} # en of i 

pvalues <- sapply(1:length(modelpgls),function(i) modelpgls[[i]]$coefficients[2,4])
names(modelpgls) <- names(pvalues) <- Intra_traits

length(which(pvalues<=0.05))

## Definition des param??tres graphiques
warm <- rainbow(length(unique(test$Family_Inter)), start=rgb2hsv(col2rgb('red'))[1], 
                end=rgb2hsv(col2rgb('yellow'))[1])
col.gp_Inter <-warm[as.factor(test$Family_Inter)] 

cool <- rainbow(25,start=rgb2hsv(col2rgb('cyan'))[1],end=rgb2hsv(col2rgb('blue'))[1])
col.gp_Intra <-cool[as.factor(test$Species)]

################################################################################
########################  Plot figure 3  #######################################
################################################################################

vect <- c(6,8,4)
Intra_traitspl <- Intra_traits[vect]
Inter_traitspl <- Inter_traits[vect]
mains <- mains_end[vect]

modelpgls_s <-  modelpgls[vect]
pvalues <-  round(sapply(1:length(modelpgls_s),function(i) modelpgls_s[[i]]$coefficients[2,4]),2)
pvalues[which(pvalues<0.05)] <- "< 0.05"

r2 <- sapply(1:length(modelpgls_s),function(i) modelpgls_s[[i]]$r.squared)

pos_list <- list(pos=c(4,1,1,1,1,3,4,4,2,2,1,3,3,4,4,1,3,3,1,3,1,3),
                 pos=c(1,2,1,1,1,3,1,3,3,1,1,3,3,4,4,3,3,3,4,3,1,3),
                 pos=c(4,4,1,1,1,2,3,1,3,3,2,1,3,4,2,1,3,1,4,3,1,3))
main_list <- c(" (PL/PH)"," (PH/BD)"," (JL/HL)")
xlab_list <- c("Intraspecific variability","Intraspecific variability","Intraspecific variability")
ylab_list <- c("Interspecific variability","","")
txt_ls <- c("a","b","c"); txt_pos <- c("2.28","1.96","2.2")

tiff("figures/Figure_3_WOCARend.tiff",width=10,height=5, units="in",res=300,compression="lzw")
#jpeg("figures/Figure_3_WOCARend.jpg",width=10,height=5, units="in",res=300)
  layout(mat=rbind(c(1,1,2,2,3,3),c(1,1,2,2,3,3),c(4,4,5,5,5,5))) 
  par(mar=c(1,1,1,1),cex.lab=0.69)
  rvect <- c(0.04,0.035,0.034)
  par(mar=c(2.5,2.5,3,1.5))
  
  for(i in 1:3){
    plot(rnorm(1000),col="white",type="n",axes=F,xlab=xlab_list[[i]],ylab=ylab_list[[i]],
         main="",cex.lab=1,mgp=c(1.4,1,0),cex.main=1.2,
         xlim=c(min(test[,Intra_traitspl[i]])-0.1,max(test[,Intra_traitspl[i]])+0.1),
         ylim=c(min(test[,Intra_traitspl[i]])-0.1,max(test[,Inter_traitspl[i]])+0.1))
    title(paste(mains[i],main_list[i],sep=""),cex.main=1.2,adj=0.5,line=0.4)
    ac <- max(test[,Inter_traitspl[i]])/max(test[,Intra_traitspl[i]])
    

    for (j in 1:length(col.gp_Inter)){
      halfCircle(x=test[j,Intra_traitspl[i]],y=test[j,Inter_traitspl[i]],start=pi/4
                 ,aspect_coef=ac,end=-3*pi/4,r=rvect[i],col=col.gp_Intra[j],lwd=0.2)
      halfCircle(x=test[j,Intra_traitspl[i]],y=test[j,Inter_traitspl[i]],lwd=0.2,
                 start=-3*pi/4,end=pi/4,aspect_coef=ac,r=rvect[i],col=col.gp_Inter[j])
    } # end of j
    
    abline(lm(test[,Inter_traitspl[i]]~test[,Intra_traitspl[i]]),col="grey31",lty=2)
    mtext(txt_ls[i],side=2,line=0,at=txt_pos[i],cex=1.1,bg="white",las=2,font=2)
    
    ## Affiche les noms des familles
    text(test[,Intra_traitspl[i]], test[,Inter_traitspl[i]],labels=substring(test$Family_Inter,1,3),
         pos=pos_list[[i]] ,cex=1,offset=0.65)
    
    ## Construction des axes
    axis(side=1,line=0,las=1,cex.axis=0.9,lwd=0.35,tcl=-0.25,bg="white",mgp=c(1,0.35,0))
    axis(side=2,line=0,las=2,cex.axis=0.9,lwd=0.35,tcl=-0.25,bg="white",mgp=c(1,0.40,0))
    box(lwd=0.35)
    mylabel <- as.expression(bquote(R^2 == .(round(r2[i],2)) ))
    legend("bottomright",legend=c(mylabel,paste("PGLS P-value ",pvalues[i],sep="")),cex=1,bty="n")
  }
  
  Datasp_fam <- cbind(col.gp_Inter,col.gp_Intra,test$Family_Inter,test$Species)
  
  Data_fam <- cbind(Datasp_fam[,3],Datasp_fam[,1])
  Data_fam <- Data_fam[order(Data_fam[,1]),]
  
  par(mar=c(0,3,1,2))
  plot(NULL,xaxt='n',yaxt='n',bty='n',ylab='',xlab='',xlim=0:1,ylim=0:1)
  legend(x=0,y=0.95,legend=unique(Data_fam[,1]),pch=19,bty="n",col=unique(unique(Data_fam[,2])),cex=1.2,ncol=2)
  mtext("Families", line=-1.2, at=0.15,cex=0.9,font=2)
  
  col_end <- cbind(gsub("_"," ", Datasp_fam[,4]),Datasp_fam[,2])
  col_end <- col_end[order(col_end[,1]),]
  plot(NULL,xaxt='n',yaxt='n',bty='n',ylab='',xlab='',xlim=0:1,ylim=0:1)
  legend(x=0,y=0.95,legend=col_end[,1],pch=19,bty="n",col=col_end[,2],cex=1,ncol=3)
  mtext("Species",line=-1.2,at=0.06,cex=0.9,font=2)

dev.off()

###############################################################################
#######################  Plot figure annexe ###################################
############################################################################### 
vect <- c(6,8,4)
Intra_traitspl <- Intra_traits[-vect]
Inter_traitspl <- Inter_traits[-vect]
mains_New <- mains_end[-vect]

# Extract P-values and R-squared values from the PGLS models.
modelpgls_s <-  modelpgls[-vect]
pvalues <-  round(sapply(1:length(modelpgls_s),function(i) modelpgls_s[[i]]$coefficients[2,4]),2)
pvalues[which(pvalues<0.05)] <- "< 0.05"
r2 <- sapply(1:length(modelpgls_s),function(i) modelpgls_s[[i]]$r.squared)

# Define positional and label lists for plotting.
pos_list <- list(pos=c(4,1,4,3,1,3,4,3,3,3,1,3,3,4,4,3,1,1,1,3,1,3),
                 pos=c(1,2,1,1,1,3,1,3,3,1,1,3,3,4,4,3,3,3,4,3,1,3),
                 pos=c(4,4,1,1,1,2,3,1,3,3,2,1,3,4,2,1,3,1,4,3,1,3),
                 pos=c(4,1,1,1,1,3,4,4,2,2,1,3,3,4,4,1,3,3,1,3,1,3),
                 pos=c(1,2,1,1,1,3,1,3,3,1,1,3,3,4,4,3,3,3,4,3,1,3),
                 pos=c(4,4,1,1,1,2,3,1,3,3,2,1,3,4,2,1,3,1,4,3,1,3),
                 pos=c(4,1,1,1,1,3,4,4,2,2,1,3,3,4,4,1,3,3,1,3,1,3),
                 pos=c(1,2,1,1,1,3,1,3,3,1,1,3,3,4,4,3,3,3,4,3,1,3),
                 pos=c(4,4,1,1,1,2,3,1,3,3,2,1,3,4,2,1,3,1,4,3,1,3),
                 pos=c(4,4,1,1,1,2,3,1,3,3,2,1,3,4,2,1,3,1,4,3,1,3))

main_list <- c(" (SL/BD)"," (PL/PH)"," (ED/HL)"," (JL/HL)"," (HL/SL)"," (Pre-DL/SL)", 
               " (Pec-L/SL)"," (Pre-ED/HL)"," (AL/SL)"," (Pre-pecL/SL)"," (DL/SL)")
xlab_list <- c("","","","","","Intra-specific variability","Intra-specific variability",
               "Intraspecific variability","Intraspecific variability","Intraspecific variability")
ylab_list <- c("Interspecific variability","","","","","Interspecific variability","","","","")
txt_ls <- letters[1:10]; 
txt_pos <- c("24","6.4","1.65","1.2","2.12",
             "0.465","2.4","1.95","1.1","1.8")

xlim <- list(c(0,3),c(0,2.8),c(0,0.8),c(0,0.55),c(0,0.6)
             ,c(-0.01,0.2),c(0,1.5),c(0,0.6),c(0,0.5),c(0,0.5))
ylimA <- c(1,0.1,0.1,0.1,0.1,0,0.1,0.1,0.1,0.1)
xlab_pos <- c(0.1,0.75,0.35,0.3,0.3)

tiff("figures/Extended_data_figure4_WOcar.tiff",width=18,height=9,units="cm",res=300,compression="lzw")
#jpeg("Figures/Sup_mat/Extended_data_figure4_WOcar.jpg",width=10,height=5, units="in",res=300)
layout(mat=rbind(c(1,1,2,2,3,3,4,4,5,5),c(1,1,2,2,3,3,4,4,5,5),
                 c(6,6,7,7,8,8,9,9,10,10),c(6,6,7,7,8,8,9,9,10,10),
                 c(11,11,11,11,12,12,12,12,12,12))) 

rvect <- c(0.63,0.16,0.038,0.028,0.05,0.013,0.05,0.05,0.03,0.05)
par(mar=c(1.5,2,1.5,0.5))
for(i in 1:10){
  plot(rnorm(1000),col="white",type="n",axes=F,xlab="",ylab=ylab_list[[i]],
       main="",cex.lab=1/2,mgp=c(1,0.5,0),
       xlim= c(xlim[[i]][1],xlim[[i]][2]),
       ylim=c(min(test[,Intra_traitspl[i]])-ylimA[i],max(test[,Inter_traitspl[i]])+0.1))
  title(paste(mains_New[i],main_list[i],sep=""),cex.main=0.55,adj=0.5,line=0.4)
  
  mtext(xlab_list[i],side=1,line=1,at=xlab_pos[i],cex=0.35,bg="white")
  ac <- max(test[,Inter_traitspl[i]])/max(test[,Intra_traitspl[i]])
  
  ## Positionnement des points
  cool <- rainbow(25, start=rgb2hsv(col2rgb('cyan'))[1],end=rgb2hsv(col2rgb('blue'))[1])
  col.gp_Intra <- cool[as.factor(test$Species)]
  
  for (j in 1:length(col.gp_Inter)){
    halfCircle(x=test[j,Intra_traitspl[i]],y=test[j,Inter_traitspl[i]],start=pi/4
               ,aspect_coef=ac,end=-3*pi/4,r=rvect[i],col=col.gp_Intra[j],lwd=0.2)
    halfCircle(x=test[j,Intra_traitspl[i]],y=test[j,Inter_traitspl[i]],lwd=0.2,
               start=-3*pi/4,end=pi/4,aspect_coef=ac,r=rvect[i],col=col.gp_Inter[j])
  } # end of j
  
  abline(lm(test[,Inter_traitspl[i]]~test[,Intra_traitspl[i]]),col="grey31",lty=2,lwd=0.5)
  mtext(txt_ls[i],side=2,line=0,at=txt_pos[i],cex=1.1/2,bg="white",las=2,font=2)
  
  ## Affiche les noms des familles
  text(test[,Intra_traitspl[i]],test[,Inter_traitspl[i]],labels=substring(test$Family_Inter,1,3),
       pos=pos_list[[i]] ,cex=1/2.5,offset=0.65/2)
  
  ## Construction des axes
  axis(side=1,line=0,las=1,cex.axis=0.9/2,lwd=0.35,tcl=-0.25,bg="white",mgp=c(0.8,0.01,0))
  axis(side=2,line=0,las=2,cex.axis=0.9/2,lwd=0.35,tcl=-0.25,bg="white",mgp=c(0.5,0.3,0))
  box(lwd=0.35)
  
  mylabel <- as.expression(bquote(R^2 == .(round(r2[i],2)) ))
  legend("bottomright",legend=c(mylabel,paste("PGLS P-value ",pvalues[i],sep="")),cex=0.4,bty="n")
} # end of i 

Datasp_fam <- cbind(col.gp_Inter,test$Family_Inter)
Datasp_fam <- Datasp_fam[order(Datasp_fam[,2]),]
Datasp_sp <- cbind(col.gp_Intra,test$Species)

par(mar=c(0,3,1,2))
plot(NULL,xaxt='n',yaxt='n',bty='n',ylab='',xlab='',xlim=0:1,ylim=0:1)
legend(x=0,y=0.95,legend=unique(Datasp_fam[,2]),pch=19,bty="n",col=unique(Datasp_fam[,1]),cex=1.2/2,ncol=2)
mtext("Families",line=-0.6,at=0,cex=0.9/2,font=2)

col_end <- cbind(gsub("_"," ", Datasp_sp[,2]),Datasp_sp[,1])
col_end <- col_end[order(col_end[,1]),]
plot(NULL,xaxt='n',yaxt='n',bty='n',ylab='',xlab='',xlim=0:1,ylim=0:1)
legend(x=0,y=0.95,legend=col_end[,1],pch=19,bty="n",col=col_end[,2],cex=1.2/2,ncol=3)
mtext("Species",line=-0.6,at=0,cex=0.9/2,font=2)

dev.off()
################################################################################
######################## end of script #########################################
################################################################################
