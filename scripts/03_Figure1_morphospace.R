################################################################################
###     Continuity in morphological disparity in tropical reef fishes        ###
###                       across evolutionary scales                         ###
###                                                                          ###
### Dr. Giulia Donati, giulia.donati@usys.ethz.ch &                          ###
### Dr Camille Albouy calbouy@ethz.ch                                        ###
###                                                                          ###  
### Script to plot the figure 1 for and  based on the PlotTangentspace       ###
###              function plots for morphological variation                  ###
###  between species in families and between individuals in species.         ###
################################################################################


# Library loading ##############################################################
lib_vect <- c("geomorph","rgl","ape","shapes","abind","splancs","FD","scales","base",
              "rredlist","dplyr","reshape2","curl","tibble","png")
sapply(lib_vect,library,character.only=TRUE)

### Function loading
source("scripts/functions/morphological_functions.R")

### Selected species
Selc_sp <- c("Caranx_melampygus","Chaetodon_trifasciatus","Chromis_atripectoralis",
             "Chromis_ternatensis","Chromis_weberi","Ctenochaetus_striatus",
             "Dascyllus_aruanus","Dascyllus_carneus","Dascyllus_trimaculatus",
             "Gomphosus_caeruleus","Halichoeres_hortulanus","Hemigymnus_fasciatus",
             "Lutjanus_kasmira","Myripristis_violacea","Oxymonacanthus_longirostris",
             "Parupeneus_Inter_sp nemus","Pseudanthias_squamipinnis")

Selc_family <- c("Acanthuridae","Carangidae","Chaetodontidae","Holocentridae",
                 "Labridae","Lutjanidae","Monacanthidae","Mullidae","Pomacentridae",
                 "Serranidae")

### Small picture loading 
Pict_nm <- list.files("figures/Pictures/Red/")
Pict <- lapply(Pict_nm, function(i){readPNG(paste("figures/Pictures/Red/",i,sep=""))})
names(Pict) <- Pict_nm 

Pict_Fnm <- list.files("figures/Pictures/Blue/")
PictF <- lapply(Pict_Fnm, function(i){readPNG(paste("figures/Pictures/Blue/",i,sep=""))})
names(PictF) <- Pict_Fnm 

###----------------- Data compilation at the species level ------------------###
### Data loading
Intra_sp <- readRDS("results/intra_sp_morpho_traits.RDS")

### Data selection
Ld2selc <- do.call(c,sapply(Selc_sp, function(i) grep(i,as.character(Intra_sp$data$genus_species))))
All_tps_Intra_sp <- Intra_sp$landmark
All_tps_Intra_sp  <- All_tps_Intra_sp[-c(6,16),,]# match Thomas claverie Landmarks

# assigning a color vector to each species
dimnames(All_tps_Intra_sp)[[3]] <- Intra_sp$data$genus_species
All_tps_Intra_sp <- All_tps_Intra_sp[,,Ld2selc]

# Generate a color palette from cyan to blue
cool <- rainbow(17, start = rgb2hsv(col2rgb('cyan'))[1], end = rgb2hsv(col2rgb('blue'))[1])
groups_mi <- cool[as.factor(dimnames(All_tps_Intra_sp)[[3]])]
names(groups_mi) <- dimnames(All_tps_Intra_sp)[[3]]
col_sp_mi <- unique(groups_mi)

### Get the shapes profile for species data
Shape_Min <- get_shapes(A=All_tps_Intra_sp)

###------------------ Data compilation at the family level ------------------###
### Data loading
Inter_sp <- readRDS("results/TC_dataset_landmarks_metadata.RDS")

### Data selection
All_tps_Inter_sp  <- Inter_sp $landmark
dimnames(All_tps_Inter_sp )[[3]] <- Inter_sp $data$Family
Famselc <- do.call(c,sapply(Selc_family, function(i) 
                   grep(i,as.character(Inter_sp $data$Family))))
All_tps_Inter_sp  <- All_tps_Inter_sp [,,Famselc]

### Assigning a color vector to each family
warm <- rainbow(10, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
groups_ma <- warm[as.factor(dimnames(All_tps_Inter_sp )[[3]])] 
names(groups_ma) <- dimnames(All_tps_Inter_sp)[[3]]
col_sp <- unique(groups_ma)

### Get the shapes profile for family data
Shape_Mac <- get_shapes(A=All_tps_Inter_sp)

############################# Plotting Intra_sp and Inter_sp  on the same plot #########################

tiff("figures/Figure_1.tiff",width=16,height=12,units="cm",res=300)
layout(mat=rbind(c(1,1,2,2),c(1,1,2,2),c(1,1,2,2),c(3,3,4,4)))

par(mar=c(3,3,2,1))
plotTangentSpace_tot(A=All_tps_Inter_sp ,gp=dimnames(All_tps_Inter_sp )[[3]],col_sp=unique(col.gp_Inter_sp ),
                     xlim=c(-0.26,0.22),ylim=c(-0.15,0.2),groups=groups_ma,legend=unique(dimnames(All_tps_Inter_sp )[[3]]),
                     flip=F,warpgrids=F,cex.pt=0.2)
mtext("a",side=2,line=1.4,at=0.34,cex=1.2,bg="white",las=2)

rasterImage(image=Pict[[1]],xleft=-0.12,ybottom=0.17,xright= -0.07,ytop=0.2)
text("Acanthuridae",x=-0.1,y=0.21,font=3,cex=0.5)
rasterImage(image=Pict[[2]],xleft=0.12,ybottom=0.15,xright=0.20,ytop=0.20)# carangue
text("Carangidae",x=0.158,y=0.14,font=3,cex=0.5)
rasterImage(image=Pict[[5]],xleft=0.155,ybottom=0.07,xright=0.215,ytop=0.1)
text("Labridae",x=0.186,y=0.063,font=3,cex=0.5)
rasterImage(image=Pict[[8]],xleft=0.18,ybottom=-0.10,xright=0.23,ytop=-0.07)
text("Serranidae",x=0.205,y=-0.104,font=3,cex=0.5)
rasterImage(image=Pict[[3]],xleft=-0.26,ybottom=-0.07,xright=-0.22,ytop=-0.025)# chat
text("Chaetodontidae",x=-0.24,y=-0.08,font=3,cex=0.5)
rasterImage(image=Pict[[4]],xleft=-0.026,ybottom=-0.155,xright=0.026,ytop=-0.12)#
text("Holocentridae",x=0,y=-0.159,font=3,cex=0.5)

arrows(min(Shape_Mac$pcdata[,1]),(0.7*max(Shape_Mac$pcdata[,2])),min(Shape_Mac$pcdata[,1]),0,length=0.05,lwd=0.8)
arrows(max(Shape_Mac$pcdata[,1]),(0.9*min(Shape_Mac$pcdata[,2])),max(Shape_Mac$pcdata[,1]),0,length=0.05,lwd=0.8)

plotTangentSpace_tot(All_tps_Intra_sp,gp=dimnames(All_tps_Intra_sp)[[3]],polyall=T,col_sp=unique(col.gp_Intra_sp),
                     xlim=c(-0.15,0.20),ylim=c(-0.1,0.15),groups=groups_mi,warpgrids=F,legend=unique(dimnames(All_tps_Intra_sp)[[3]]),
                     cex.pt=0.2)
mtext("b",side=2,line=1.4,at=0.25,cex=1.2,bg="white",las=2)

rasterImage(image=PictF[[1]],xleft=-0.05,ybottom=-0.09,xright= -0.01,ytop=-0.07)
text("C. weberi",x=-0.03,y=-0.095,font=3,cex=0.5)
rasterImage(image=PictF[[2]],xleft=-0.05,ybottom=0.14,xright=-0.01,ytop=0.17)
text("C. striatus",x=-0.036,y=0.18,font=3,cex=0.5)
rasterImage(image=PictF[[3]],xleft=-0.15,ybottom=-0.02,xright=-0.13,ytop=-0)
text("D. carneus",x=-0.14,y=-0.026,font=3,cex=0.5)
rasterImage(image=PictF[[4]],xleft=-0.10,ybottom=-0.05,xright=-0.08,ytop=-0.03)
text("D. aruanus",x=-0.09,y=-0.055,font=3,cex=0.5)
rasterImage(image=PictF[[7]],xleft=0.04,ybottom=-0.125,xright=0.09,ytop=-0.095)
text("P. Inter_sp nemus",x=0.07,y=-0.13,font=3,cex=0.5)
rasterImage(image=PictF[[5]],xleft=0.125,ybottom=-0.095,xright=0.18,ytop=-0.07)
text("G. caeruleus",x=0.16,y=-0.097,font=3,cex=0.5)
rasterImage(image=PictF[[6]],xleft=0.175,ybottom=0.08,xright=0.21,ytop=0.1)
text("O. longirostris",x=0.185,y=0.105,font=3,cex=0.5)

arrows(min(Shape_Min$pcdata[,1]), (0.7 * max(Shape_Min$pcdata[,2])), min(Shape_Min$pcdata[,1]),0,length=0.05,lwd=0.8)
arrows(max(Shape_Min$pcdata[,1]), (0.7 * min(Shape_Min$pcdata[,2])), max(Shape_Min$pcdata[,1]),0,length=0.05,lwd=0.8)

par(mar=c(0,3,0,2))
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='',xlim=0:1,ylim=0:1)
legend(x=0.25,y=1,legend=unique(names(groups_ma)),pch=19,bty="n",col=unique(groups_ma),cex=.8,ncol=2)
#mtext("Species", line=-38, at=-0.26,cex=1)

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend(x=0,y=1,legend=gsub("_"," ",unique(names(groups_mi))),pch=19,bty="n",col=unique(groups_mi),cex=.7,ncol=2)
#mtext("Families", line=-38, at=-0.26,cex=1)
par(fig=c(0.001,0.25,0.8,0.9),new=T,mar=c(0,0,0,0))
geomorph:::tps(Shape_Mac$ref, Shape_Mac$shape.min,20,sz=0.5,grid.lwd=0.5)

par(fig=c(0.33,0.48,0.3,0.6),new=T,mar=c(0,0,0,0))
geomorph:::tps(Shape_Mac$ref, Shape_Mac$shape.max,20,sz=0.5,grid.lwd=0.5)

par(fig=c(0.5,0.75,0.8,0.9),new=T,mar=c(0,0,0,0))
geomorph:::tps(Shape_Min$ref,Shape_Min$shape.min,20,sz=0.5,grid.lwd=0.5)

par(fig=c(0.83,0.98,0.3,0.5),new=T,mar=c(0,0,0,0))
geomorph:::tps(Shape_Min$ref,Shape_Min$shape.max,20,sz=0.5,grid.lwd=0.5)

dev.off()
################################################################ end of script