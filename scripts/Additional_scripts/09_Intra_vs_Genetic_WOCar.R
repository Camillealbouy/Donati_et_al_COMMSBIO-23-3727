################################################################################
###.          Continuity in morphological disparity in tropical reef         ###
###                     fishes across evolutionary scales                    ###
###                       without Carangidae species                         ###
### Dr. Giulia Donati, giulia.donati@usys.ethz.ch &                          ###
### Dr Camille Albouy calbouy@ethz.ch                                        ###
###                                                                          ###
################################################################################
###  Figure 4: Intra-specific variability vs genetic diversity               ###
################################################################################

### Library
lib_vect <- c("ade4","ape","nlme","geiger","caper","phytools","phylobase","dplyr","png")
sapply(lib_vect,library,character.only=TRUE)

### Select the species
Selc_sp <- c("Chaetodon_trifasciatus","Chromis_atripectoralis","Chromis_ternatensis", 
             "Chromis_weberi","Ctenochaetus_striatus","Dascyllus_aruanus","Dascyllus_carneus",
             "Dascyllus_trimaculatus","Gomphosus_caeruleus","Halichoeres_hortulanus",
             "Hemigymnus_fasciatus","Lutjanus_kasmira","Myripristis_violacea","Oxymonacanthus_longirostris"
             ,"Parupeneus_macronemus","Pseudanthias_squamipinnis")

### Select the family              
Family <- c("Chaetodontidae","Pomacentridae","Pomacentridae","Pomacentridae",
            "Acanthuridae","Pomacentridae","Pomacentridae","Pomacentridae",
            "Labridae","Labridae","Labridae","Lutjanidae","Holocentridae",
            "Monacanthidae","Mullidae","Serranidae")

################################################################################
###                       Data compilation at the species level.             ###
################################################################################
Intra <- readRDS("results/Intra_specific_dataset.RDS")
Intra <- Intra$intra_variability_pop[Selc_sp,]

Intra <- cbind(Family=Family,nb_fam=as.vector(table(Family)[Family]),Intra)
colnames(Intra) <- paste(colnames(Intra),"Intra",sep="_")

Genet_Data <- readRDS("data/F_dataset_2.RDS")
Genet_Data <-  Genet_Data[Genet_Data$Species %in% Selc_sp,]

################################################################################
###           Pgl modèle between JT and intra specific varibility            ###
################################################################################
### Traits and Genetic data Preparation
F_data <- data.frame(Genet_Data,Intra,row.names=rownames(Intra),Sp_names=rownames(Genet_Data))

### Phylogenetic Data
Rabosky <- read.tree("data/phylo/actinopt_12k_treePL.tre")
Rabosky$tip.label <- gsub("Gomphosus_varius","Gomphosus_caeruleus",Rabosky$tip.label)
Rabosky$tip.label <- gsub(Rabosky$tip.label[grep("Acanthurus",Rabosky$tip.label)][1],
                          "Acanthurus_tractus",Rabosky$tip.label)
pruned.treeRabosky_sp <- drop.tip(Rabosky,setdiff(Rabosky$tip.label, rownames(Intra))) 

comp.data <- comparative.data(pruned.treeRabosky_sp,F_data,names.col=Species,
                              vcv.dim=2,warn.dropped=T)

################################################################################
###           Pgl modèle between JT and intra specific varibility            ###
################################################################################
ratio <- colnames(F_data)[grep("ratio_",colnames(F_data))]
Intra_traits <- ratio[grep("_Intra",as.character(ratio))]

modelpgls <- lapply(1:length(Intra_traits),function(i){
  summary(pgls(F_data$Jt~F_data[,Intra_traits[i]],comp.data)) 
})
names(modelpgls) <- Intra_traits

pvalues <- sapply(1:length(modelpgls),function(i) modelpgls[[i]]$coefficients[2,4])
names(pvalues) <- Intra_traits
length(which(pvalues<=0.05))

mains <- c("Body elongation","Caudal peduncle elongation","Eye size","Jaw length",
           "Head length", "Caudal peduncle height","Dorsal fin size","Anal fin size",
           "Pectoral fin size","Eye position","Anal fin insertion",
           "Pectoral fin insertion","Dorsal fin insertion")

################################################################################
###        Draw the plot Figure 4 without Carangidae main plot               ###   
################################################################################

### Selected traits 
vect <- c(8,2,13)
Intra_traitspl <- Intra_traits[vect]
mains_s <- mains[vect]

modelpgls_s <-  modelpgls[vect]
pvalues_s <- pvalues[vect]
pvalues_s[which(round(pvalues_s,3)<0.05)] <- "< 0.05"

r2 <- sapply(1:length(modelpgls_s),function(i) modelpgls_s[[i]]$r.squared)


pos_list <- list(pos=c(3,3,1,4,1,3,4,4,4,2,1,3,3,4,4,1,3),
                 pos=c(1,2,1,1,1,3,1,3,3,1,1,3,3,4,4,3,3),
                 pos=c(4,4,2,4,1,2,3,1,4,2,4,1,3,3,2,1,4))
main_list <- c(" (AL/SL)"," (PL/PH)"," (DL/SL)")
xlab_list <- c("Intraspecific variability","Intraspecific variability","Intraspecific variability")
expylab  <-as.expression(bquote('J'['T'] ))
ylab_list <- c(expylab,"","")
txt_ls <- c("a","b","c"); txt_pos <- rep("1.47",3)
xlab_pos <- c(0.325,1.25,0.3)

### Colour definition
cool <- rainbow(25, start=rgb2hsv(col2rgb('cyan'))[1],end=rgb2hsv(col2rgb('blue'))[1])
col.gp_Intra <- cool[as.factor(F_data$Species)]

warm <- rainbow(length(unique(F_data$Family)), start=rgb2hsv(col2rgb('red'))[1], 
                end=rgb2hsv(col2rgb('yellow'))[1])
col.gp_macro <-warm[F_data$Family] 


tiff("figures/Extended_figure_10_woCAr.tiff",width=10,height=4,units="in",res=300,compression="lzw")
#jpeg("Figures/Figure_4/Figure_4_WCAR.jpg",width=10,height=5, units="in",res=300)
layout(mat=rbind(c(1,1,2,2,3,3),c(1,1,2,2,3,3),c(1,1,2,2,3,3),c(1,1,2,2,3,3),c(4,4,4,4,4,4))) 

par(mar=c(2.5,2.8,3,1.5),cex.lab=0.6)
for(i in 1:3){
  plot(rnorm(1000),col="white",type="n",axes=F,xlab="",ylab=ylab_list[[i]],
       main="",cex.lab=1,mgp=c(1.8,1,0),
       xlim=c(min(F_data[,Intra_traitspl[i]])-0.1,max(F_data[,Intra_traitspl[i]])+0.2),
       ylim=c(min(F_data$Jt-0.01),max(F_data$Jt)+0.01))
  
  title(paste(mains_s[i],main_list[i],sep=""),cex.main=1.2,adj=0.5,line=0.4)
  mtext(xlab_list[i],side=1,line=1.5,at=xlab_pos[i],cex=0.75,bg="white")
  
  points(y=Genet_Data$Jt,x=F_data[,Intra_traitspl[i]],bg=col.gp_macro,pch=21,col="black")
  abline(lm(Genet_Data$Jt~F_data[,Intra_traitspl[i]]),col="grey31",lty=2)
  text(F_data[,Intra_traitspl[i]],Genet_Data$Jt,labels=F_data$Sp_names,pos=pos_list[[i]],cex=0.7,offset=0.6)
  
  axis(side=1,line=0,las=1,cex.axis=0.9,lwd=0.35,tcl=-0.25,bg="white",mgp=c(1,0.35,0))
  axis(side=2,line=0,las=2,cex.axis=0.9,lwd=0.35,tcl=-0.25,bg="white",mgp=c(1,0.40,0))
  box(lwd=0.5)
  mtext(txt_ls[i],side=2,line=1,at=txt_pos[i],cex=1.1,bg="white",las=2,font=2)
  
  mylabel <- as.expression(bquote(R^2 == .(round(r2[i],2)) ))
  legend("bottomright",legend=c(mylabel,paste("PGLS P-value ",pvalues_s[i],sep="")),cex=0.8,bty="n")
}

Datasp_fam <- cbind(col.gp_macro,as.character(F_data$Family))
Datasp_fam <- Datasp_fam[order(Datasp_fam[,2]),]

par(mar=c(0,3,1,2))
plot(NULL,xaxt='n',yaxt='n',bty='n',ylab='',xlab='',xlim=0:1,ylim=0:1)
legend(x=0.25,y=0.95,legend=unique(Datasp_fam[,2]),pch=19,bty="n",col=unique(unique(Datasp_fam[,1])),
       cex=0.9,ncol=5)
mtext("Families", line=-0.7, at=0.5,cex=0.65,font=2)
dev.off()

################################################################################
####                         Plot figure annexe                             ####
################################################################################
vect <- c(8,2,13)
Intra_traitspl <- Intra_traits[-vect]
mains_s <- mains[-vect]

modelpgls_s <-  modelpgls[-vect]
pvalues_s <- round(sapply(1:length(modelpgls_s),function(i) modelpgls_s[[i]]$coefficients[2,4]),2)
pvalues_s[which(round(pvalues_s,3)<0.05)] <- "< 0.05"

r2 <- sapply(1:length(modelpgls_s),function(i) modelpgls_s[[i]]$r.squared)

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

main_list <- c(" (SL/BD)"," (ED/HL)"," (JH/HL)"," (HL/SL)"," (PH/BD)"," (Pre-DL/SL)", 
               " (Pec-L/SL)"," (Pre-ED/HL)"," (Pre-AL/SL)"," (Pre-pecL/SL)")
xlab_list <- c("","","","","","Intraspecific variability","Intraspecific variability",
               "Intraspecific variability","Intraspecific variability","Intraspecific variability")

expylab  <-as.expression(bquote('J'['T'] ))
ylab_list <- c(expylab,"","","","",expylab,"","","","")
txt_ls <- letters[1:10];  txt_pos <- rep("1.48",10)
ylimA <- c(1,0.1,0.1,0.1,0.1,0,0.1,0.1,0.1,0.1)
xlab_pos <- c(0.1,0.75,0.35,0.3,0.3)

tiff("figures/Figure_4_W0Car.tiff",width=18,height=9, units="cm",res=300,compression ="lzw")
#jpeg("figures/Extended_figure_11_woCAr.tiff",width=18,height=9, units="cm",res=300)

layout(mat=rbind(c(1,1,2,2,3,3,4,4,5,5),c(1,1,2,2,3,3,4,4,5,5),
                 c(6,6,7,7,8,8,9,9,10,10),c(6,6,7,7,8,8,9,9,10,10),
                 c(11,11,11,11,11,11,11,11,11,11))) 

par(mar=c(1.5,2,1.5,0.5))
for(i in 1:10){
  plot(rnorm(1000),col="white",type="n",axes=F,xlab="",ylab=ylab_list[[i]],
       main="",cex.lab=0.7,mgp=c(1.1,1,0),
       xlim=c(min(F_data[,Intra_traitspl[i]])-0.1,max(F_data[,Intra_traitspl[i]])+0.2),
       ylim=c(min(F_data$Jt-0.05),max(F_data$Jt)+0.01))
  
  title(paste(mains_s[i],main_list[i],sep=""),cex.main=0.55,adj=0.5,line=0.4)
  mtext(xlab_list[i],side=1,line=0.9,at=xlab_pos[i],cex=0.4,bg="white")
  
  points(y=Genet_Data$Jt,x=F_data[,Intra_traitspl[i]],bg=col.gp_macro,pch=21,col="black",
         cex=0.6,lwd=0.3)
  
  abline(lm(Genet_Data$Jt~F_data[,Intra_traitspl[i]]),col="grey31",lty=2,lwd=0.4)
  text(F_data[,Intra_traitspl[i]],Genet_Data$Jt,labels=F_data$Sp_names,pos=pos_list[[i]],cex=0.4,offset=0.3)
  
  axis(side=1,line=0,las=1,cex.axis=0.4,lwd=0.35,tcl=-0.25,bg="white",mgp=c(1,0.1,0))
  axis(side=2,line=0,las=2,cex.axis=0.4,lwd=0.35,tcl=-0.25,bg="white",mgp=c(1,0.40,0))
  box(lwd=0.5)
  mtext(txt_ls[i],side=2,line=0,at=txt_pos[i],cex=0.6,bg="white",las=2,font=2)
  
  ## Test de correlation and legend
  mylabel <- as.expression(bquote(R^2 == .(round(r2[i],2)) ))
  legend("bottomright",legend=c(mylabel,paste("PGLS P-value ",pvalues_s[i],sep="")),cex=0.45,bty="n")
}

Datasp_fam <- cbind(col.gp_macro,as.character(F_data$Family))
Datasp_fam <- Datasp_fam[order(Datasp_fam[,2]),]

par(mar=c(0,3,1,2))
plot(NULL,xaxt='n',yaxt='n',bty='n',ylab='',xlab='',xlim=0:1,ylim=0:1)
legend(x=0.18,y=0.92,legend=unique(Datasp_fam[,2]),pch=19,bty="n",col=unique(unique(Datasp_fam[,1])),
       cex=0.7,ncol=5)
mtext("Families", line=-0.9, at=0.5,cex=0.6,font=2)
dev.off()

################################################################################ 
######################## end of script #########################################
################################################################################