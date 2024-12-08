#function to load to read tps data
read.tps <-  function(data) {
  a <- readLines(data)
  LM <- grep("LM", a)
  ID.ind <- grep("ID", a)  
  images <- basename(gsub("(IMAGE=)(.*)", "\\2", a[ID.ind - 1]))
  skip <- LM
  nrows <- as.numeric(gsub("(LM=)([0-9])", "\\2", grep("LM", a, value=T)))
  l <- length(LM)
  landmarks <- vector("list", l)
  
  for (i in 1:l) {
    cat("i=",i,"\n")
    landmarks[i] = list(data.frame(
      read.table(file=data, header=F, skip=LM[i],
                 nrows=nrows[i], col.names=c("X", "Y")),
      IMAGE = images[i],
      ID = read.table(file=data, header=F, skip=ID.ind[i]-1, 
                      nrows=1, sep="=", col.names="ID")[2,],
      Scale = read.table(file=data, header=F, skip=ID.ind[i],
                         nrows=1, sep="=")[,2]))
  }
  do.call(rbind, landmarks)
}

# Function to extract Landmark data and curve points data from .TPS per image
readland.tps <- function (file="/Volumes/LaCie/MSc_Camille/Project_morpho/Maldives/TPS/Amphiprion_clarkii_MV.TPS", specID = c("ID","IMAGE"), 
                          readcurves=T,warnmsg=T) {
  ignore.case = TRUE
  specID <- match.arg(specID)
  tpsfile <- scan(file = file, what = "char", sep = "\n", quiet = TRUE)
  lmdata <- grep("LM=", tpsfile, ignore.case)
  if (length(lmdata !=0)) {
    nland <- as.numeric(sub("LM=", "", tpsfile[lmdata], ignore.case))
    k <- 2
  }
  if (length(lmdata) == 0) {
    lmdata <- grep("LM3=", tpsfile, ignore.case)
    nland <- as.numeric(sub("LM3=", "", tpsfile[lmdata], ignore.case))
    k <- 3
  }
  if(any(nland == 0)){ stop("No landmark data for some specimens.") }
  n <- nspecs <- length(lmdata)
  if (max(nland) - min(nland) != 0) {
    stop("Number of landmarks not the same for all specimens.")
  }
  p <- nland[1]
  imscale <- as.numeric(sub("SCALE=", "", tpsfile[grep("SCALE", 
                                                       tpsfile, ignore.case)], ignore.case))
  if (is.null(imscale)) {
    imscale = array(1, nspecs)
  }
  if (warnmsg == T) {
    if (length(imscale) != nspecs) {
      print("Not all specimens have scale. Assuming landmarks have been previously scaled.")
    }
  }
  if (length(imscale) != nspecs) {
    imscale = array(1, nspecs)
  }
  crvs <- grep("CURVES=", tpsfile, ignore.case)
  if(length(crvs)>0){
    if (readcurves == TRUE && length(crvs) == 0){ stop("No CURVES= field present in file") } 
    ncurve <- as.numeric(sub("CURVES=", "", tpsfile[crvs], ignore.case))
    ncurvepts <- as.numeric(sub("POINTS=", "", tpsfile[grep("POINTS=", tpsfile, ignore.case)], ignore.case))
    if (max(ncurve) - min(ncurve) != 0) {
      stop("Number of curves not the same for all specimens.") }
    if (warnmsg == T && readcurves==T) {print(paste("Landmarks 1:", p, " are fixed landmarks.", sep=""))
      print(paste("Landmarks ", p+1, ":", p+sum(ncurvepts[1:ncurve[1]]), " are semilandmarks.", sep=""))}
    p <- nland[1] + sum(ncurvepts[1:ncurve[1]]) 
  }    
  tmp <- tpsfile[-(grep("=", tpsfile))]
  options(warn = -1)
  tmp <- matrix(as.numeric(unlist(strsplit(tmp,"\\s+"))),ncol = k, byrow = T)
  
  if (warnmsg == T) {
    if (sum(which(is.na(tmp) == TRUE)) > 0) {
      print("NOTE.  Missing data identified.")
    }
  }
  coords <- aperm(array(t(tmp), c(k, p, n)), c(2, 1, 3))
  imscale <- aperm(array(rep(imscale, p * k), c(n, k, p)), 
                   c(3, 2, 1))
  coords <- coords * imscale
  if (readcurves==F){coords<-coords[1:nland,,] 
  if(n==1) coords <- array(coords, c(nland,k,n))}
  if (specID == "None") {
    if (warnmsg == T) {print("No Specimen names extracted")
    }
  }
  if (specID == "IMAGE") {
    imageID <- (sub("IMAGE=", "", tpsfile[grep("IMAGE=", tpsfile, ignore.case)], 
                    ignore.case))
    if (length(imageID) != 0) {
      imageID <- sub(".jpg", "", imageID, ignore.case)
      imageID <- sub(".tif", "", imageID, ignore.case)
      imageID <- sub(".bmp", "", imageID, ignore.case)
      imageID <- sub(".tiff", "", imageID, ignore.case)
      imageID <- sub(".jpeg", "", imageID, ignore.case)
      imageID <- sub(".jpe", "", imageID, ignore.case)
      dimnames(coords)[[3]] <- as.list(imageID)
      if (warnmsg == T) {
        print("Specimen names extracted from line IMAGE=")
      }
    }
    if (length(imageID) == 0) {
      if (warnmsg == T) {
        print("No name given under IMAGE=. Specimen names not extracted")
      }
    } 
  }
  if (specID == "ID") {
    ID <- sub("ID=", "", tpsfile[grep("ID=", tpsfile, ignore.case)], ignore.case)
    if (length(ID) == 0) {
      if(warnmsg ==T){
        print("No name given under ID=. Specimen names not extracted")
      }
    }
    if (length(ID) != 0) {
      dimnames(coords)[[3]] <- as.list(ID)
      if (warnmsg == T) {
        print("Specimen names extracted from line ID=")
      }
    }
  }
  return(coords = coords)                    
} # end of readland.tps
# Function de transformation des donn??es brutes en donn??es 3D
transfor_array = function (Data=essai,ncol=2,nrow=18,nbdim=115){
  ncol = ncol
  nrow = nrow 
  nbdim = nbdim
  
  w=array(NA,c(nrow,ncol,nbdim))
  
  # remplissage du tableau 
  
  ind=unique(Data$Individu)
  
  for (i in 1:nbdim){
    cat("i =",i,"\n")
    z=subset(Data, Data$Individu == ind[i])[,4:5]
    w [,1,i] = z[,1]
    w[,2,i] = z[,2]
    
  }  
  return(w)
} # transform_array
################################################################################
get_shape_by_lag <- function (res_pro=a$rotated,nb_ld_marks=18,nb_ind_init=1,nb_ind_final=42){
  g=rep(0,nb_ld_marks)
  g2=rep(0,nb_ld_marks)
  
  for (i in nb_ind_init: nb_ind_final){
    g= cbind(g,res_pro[,1,i])
    g2=cbind(g2,res_pro[,2,i])
    
  }
  
  g=g[,-1]
  g2=g2[,-1]
  
  coordx= apply(g,1,mean)
  coordy= apply(g2,1,mean)
  
  return(data.frame(mean_coordx = coordx, mean_coordy=coordy))
}
################################################################################
centsiz = function(M) {p<-dim(M)[1]
size<-sqrt(sum(apply(M,2,var))*(p-1))
list("centroid_size"=size,"scales"=M/size)
}
################################################################################
mshape <- function(A){apply(A,c(1,2),mean)}
################################################################################
pPsup <- function(M1, M2){
  k<-ncol (M1)
  Z1<-trans1(centsiz(M1)[[2]])
  Z2<-trans1(centsiz(M2)[[2]])
  
  sv<-svd(t(Z2)%*%Z1)
  U<-sv$v;V<-sv$u;
  Delt<-sv$d
  sig<-sign(det (t(Z1)%*%Z2))
  Delt[k]<-sig*abs(Delt[k]) ;V[,k]<-sig*V[,k]
  Gam<-U%*%t(V)
  
  beta<-sum(Delt)
  
  list(Mp1=Z1%*%Gam,Mp2=Z2,rotation=Gam,DP=sqrt(sum(ild2(Z1%*%Gam,Z2)^2)),rho=acos(beta))
  
}
################################################################################
ild2 <- function(M1,M2){sqrt(apply((M1-M2)^2, 1,sum))}
################################################################################
orp <- function(A){
  p <- dim(A)[1]; k <- dim(A)[2];
  n <- dim(A)[3]
  Y1 <- as.vector(centsiz(mshape(A))[[2]])
  
  oo<-as.matrix(rep(1,n))%*%Y1
  I<-diag(1,k*p)
  mat<-matrix(NA,n,k*p)
  
  for(i in 1:n){
    mat[i,]<-as.vector(A[,,i])}
    Xp <- mat%*%(I-(Y1%*%t(Y1)))
    Xp1<- Xp+oo
    array(t(Xp1),dim=c(p,k,n))
}
################################################################################
pgpa <- function(A) {p<-dim(A)[1];k<-dim(A)[2];n<-dim(A)[3]
temp2<-temp1<-array(NA, dim=c(p,k,n)); Siz<-numeric(n)
for (i in 1:n) {
  Acs<-centsiz(A[,,i])
  Siz[i]<-Acs[[1]]
  temp1[,,i]<-trans1(Acs[[2]])}
  Qm1<-dist(t(matrix(temp1,k*p,n)))
  Q<-sum(Qm1); iter<-0
while (abs(Q)>0.001)
{for(i in 1:n){
  M<-mshape(temp1[,,-i])
  temp2[,,i]<-pPsup(temp1[,,i],M)[[1]]}
  Qm2<-dist(t(matrix(temp2,k*p,n)))
  Q<-sum(Qm1)-sum(Qm2)
  Qm1<-Qm2
  iter=iter+1
  temp1<-temp2}
list("rotated"=temp2,"it.number"=iter,"Q"=Q,"intereucl.dist"=Qm2,"mshape"=centsiz(mshape(temp2))[[2]],"cent.size"=Siz)}
################################################################################
trans1<-function(M){scale(M,scale=F)}
################################################################################
## Functions for graphs
shape_habitat_plot = function(data=data,main=main,col=col) {
  plot(rnorm(100),col="white",xlim=c(-0.15,0.25),ylim=c(-0.15,0.1),type="n",axes=TRUE,xlab="",ylab="", main=main)
  
  #Landmarks
  for (i in 1:17) points(data[i,1],data[i,2],col=col,pch=20,cex=1)
  
  segments(data[1,1],data[1,2],data[7,1],data[7,2],col=col)
  segments(data[2,1],data[2,2],data[11,1],data[11,2],col=col)
  
  #Curve 1
  segments(data[1,1],data[1,2],data[18,1],data[18,2], col=col)
  for (n in 18:29){
    segments(data[n,1],data[n,2],data[n+1,1],data[n+1,2], col=col)
  } 
  segments(data[30,1],data[30,2],data[2,1],data[2,2], col=col)
  
  #Curve 2
  segments(data[2,1],data[2,2],data[31,1],data[31,2], col=col)
  for (n in 31:47){
    segments(data[n,1],data[n,2],data[n+1,1],data[n+1,2], col=col)
  } 
  segments(data[48,1],data[48,2],data[5,1],data[5,2], col=col)
  
  #Curve 3
  segments(data[11,1],data[11,2],data[49,1],data[49,2], col=col)
  for (n in 49:60){
    segments(data[n,1],data[n,2],data[n+1,1],data[n+1,2], col=col)
  } 
  segments(data[61,1],data[61,2],data[1,1],data[1,2], col=col)
  
  
  
}

### Flip geomorph
PlotTangentSpace_flip <- function (A, axis1 = 1, axis2 = 2, warpgrids = TRUE, mesh = NULL, 
                                   label = NULL, groups = NULL, legend = FALSE, flip=TRUE, ...) 
{
  if (length(dim(A)) != 3) {
    stop("Data matrix not a 3D array (see 'arrayspecs').")
  }
  if (any(is.na(A)) == T) {
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').")
  }
  dots <- list(...)
  retx <- dots$retx
  if (is.null(retx)) 
    retx <- TRUE
  scale. <- dots$scale.
  if (is.null(scale.)) 
    scale. <- FALSE
  center <- dots$center
  if (is.null(center)) 
    center <- TRUE
  tol <- dots$tol
  k <- dim(A)[2]
  p <- dim(A)[1]
  n <- dim(A)[3]
  ref <- mshape(A)
  x <- two.d.array(A)
  if (is.null(tol)) {
    d <- prcomp(x)$sdev^2
    cd <- cumsum(d)/sum(d)
    cd <- length(which(cd < 1))
    if (length(cd) < length(d)) 
      cd <- cd + 1
    if (length(d) > 2) 
      tol <- max(c(d[cd]/d[1], 0.005))
    else tol <- 0
  }
  pc.res <- prcomp(x, center = center, scale. = scale., retx = retx, 
                   tol = tol)
  pcdata <- pc.res$x
  if(flip==TRUE){
    pcdata[,axis2] <- (pcdata[, axis2])*-1
  }
  if (warpgrids == FALSE) {
    if (legend == TRUE) {
      layout(t(matrix(c(1, 1, 2, 1, 1, 1, 1, 1, 1), 3, 
                      3)))
    }
    plot(pcdata[, axis1],  pcdata[, axis2], asp = 1, pch = 21, 
         bg = "black", cex = 2, xlab = paste("PC ", axis1), 
         ylab = paste("PC ", axis2))
    if (!is.null(groups)) {
      points(pcdata[, axis1], pcdata[, axis2], pch = 21, 
             bg = groups, cex = 2)
    }
    segments(min(pcdata[, axis1]), 0, max(pcdata[, axis1]), 
             0, lty = 2, lwd = 1)
    segments(0, min(pcdata[, axis2]), 0, max(pcdata[, axis2]), 
             lty = 2, lwd = 1)
    if (length(label != 0)) {
      if (isTRUE(label)) {
        text(pcdata[, axis1], pcdata[, axis2], seq(1, 
                                                   n), adj = c(-0.7, -0.7))
      }
      else {
        text(pcdata[, axis1], pcdata[, axis2], label, 
             adj = c(-0.1, -0.7))
      }
    }
    if (!is.null(groups) && legend == TRUE) {
      plot.new()
      if (is.factor(groups)) {
        legend(0.5, 1, legend = unique(groups), pch = 19, 
               bty = "n", col = unique(groups))
      }
      else {
        legend(0.5, 1, legend = unique(names(groups)), 
               pch = 19, bty = "n", col = unique(groups))
      }
    }
  }
  
  shapes <- shape.names <- NULL
  for (i in 1:ncol(pcdata)) {
    pcaxis.min <- min(pcdata[, i])
    pcaxis.max <- max(pcdata[, i])
    pc.min <- pc.max <- rep(0, dim(pcdata)[2])
    pc.min[i] <- pcaxis.min
    pc.max[i] <- pcaxis.max
    pc.min <- as.matrix(pc.min %*% (t(pc.res$rotation))) + 
      as.vector(t(ref))
    pc.max <- as.matrix(pc.max %*% (t(pc.res$rotation))) + 
      as.vector(t(ref))
    shapes <- rbind(shapes, pc.min, pc.max)
    shape.names <- c(shape.names, paste("PC", i, "min", sep = ""), 
                     paste("PC", i, "max", sep = ""))
  }
  shapes <- arrayspecs(shapes, p, k)
  shapes <- lapply(seq(dim(shapes)[3]), function(x) shapes[, 
                                                           , x])
  names(shapes) <- shape.names
  if (warpgrids == TRUE) {
    if (k == 2) {
      layout(t(matrix(c(2, 1, 4, 1, 1, 1, 1, 1, 3), 3, 
                      3)))
    }
    plot(pcdata[, axis1], pcdata[, axis2], asp = 1, pch = 21, 
         bg = "black", cex = 2, xlab = paste("PC ", axis1), 
         ylab = paste("PC ", axis2))
    if (!is.null(groups)) {
      points(pcdata[, axis1], pcdata[, axis2], pch = 21, 
             bg = groups, cex = 2)
    }
    segments(min(pcdata[, axis1]), 0, max(pcdata[, axis1]), 
             0, lty = 2, lwd = 1)
    segments(0, min(pcdata[, axis2]), 0, max(pcdata[, axis2]), 
             lty = 2, lwd = 1)
    if (length(label != 0)) {
      if (isTRUE(label)) {
        text(pcdata[, axis1], pcdata[, axis2], seq(1, 
                                                   n), adj = c(-0.7, -0.7))
      }
      else {
        text(pcdata[, axis1], pcdata[, axis2], label, 
             adj = c(-0.1, -0.1))
      }
    }
    shape.min <- shapes[[which(names(shapes) == paste("PC", 
                                                      axis1, "min", sep = ""))]]
    shape.max <- shapes[[which(names(shapes) == paste("PC", 
                                                      axis1, "max", sep = ""))]]
    if (k == 2) {
      arrows(min(pcdata[, axis1]), (0.7 * max(pcdata[, 
                                                     axis2])), min(pcdata[, axis1]), 0, length = 0.1, 
             lwd = 2)
      arrows(max(pcdata[, axis1]), (0.7 * min(pcdata[, 
                                                     axis2])), max(pcdata[, axis1]), 0, length = 0.1, 
             lwd = 2)
      geomorph:::tps(ref, shape.min, 20)
      geomorph:::tps(ref, shape.max, 20)
    }
    if (!is.null(groups) && legend == TRUE) {
      plot.new()
      if (is.factor(groups)) {
        legend(0.5, 1, legend = unique(groups), pch = 19, 
               bty = "n", col = unique(groups))
      }
      else {
        legend(0.5, 1, legend = unique(names(groups)), 
               pch = 19, bty = "n", col = unique(groups))
      }
    }
    if (k == 3) {
      if (is.null(mesh) == TRUE) {
        open3d()
        mfrow3d(1, 2)
        plot3d(shape.min, type = "s", col = "gray", main = paste("PC ", 
                                                                 axis1, " negative"), size = 1.25, aspect = FALSE, 
               xlab = "", ylab = "", zlab = "", box = FALSE, 
               axes = FALSE)
        plot3d(shape.max, type = "s", col = "gray", main = paste("PC ", 
                                                                 axis1, " positive"), size = 1.25, aspect = FALSE, 
               xlab = "", ylab = "", zlab = "", box = FALSE, 
               axes = FALSE)
      }
      if (is.null(mesh) == FALSE) {
        open3d()
        mfrow3d(1, 2)
        cat(paste("\nWarping mesh to negative end of axis ", 
                  axis1, "\n", sep = ""))
        plotRefToTarget(ref, shape.min, mesh, method = "surface")
        title3d(main = paste("PC ", axis1, " negative"))
        next3d()
        cat(paste("\nWarping mesh to positive end of axis ", 
                  axis1, "\n", sep = ""))
        plotRefToTarget(ref, shape.max, mesh, method = "surface")
        title3d(main = paste("PC ", axis1, " positive"))
      }
    }
    layout(1)
  }
  out <- list(pc.summary = summary(pc.res), pc.scores = pcdata, 
              pc.shapes = shapes, sdev = pc.res$sdev, rotation = pc.res$rotation)
  class(out) = "plotTangentSpace"
  out
}


#================== Indicator function calculation ============================#
gb_Mric_function <- function(pcdata=pcdata,ax=c(axis1,axis2)){
  coord_d<-pcdata[,ax]
  rs <- nrow(coord_d)
  
  # Morphological Richness
  MRic <- round(convhulln(coord_d,"FA")$vol,10)
  
  # identity of vertices
  #vert0 <- melt(convhulln(coord_d,"Fx"))[,3]
  #vert2 <- unique(vert0)
  vert0<-convhulln(coord_d,"Fx TO 'vert.txt'")
  vert1<-scan("vert.txt",quiet=T)
  vert2<-(vert1+1)[-1]
  
  
  ME_vert_k<- row.names(coord_d)[vert2] 
  ME_vertices_cells <- ME_vert_k[!is.na(ME_vert_k)] # arendre
  
  # number of vertices
  nbME_vertices<-length(ME_vert_k)
  
  ### Calcul du FDiv
  # coord values of vertices
  trvertices<-coord_d[vert2,]
  
  # coordinates of the center of gravity of the vertices (B)
  B<-apply(trvertices,2,function(x) {mean(x, na.rm=TRUE)})
  
  # computing euclidian dstances to B (dB)
  dB<-apply(coord_d, 1, function(x) { (sum((x-B)^2) )^0.5} )
  
  # mean of dB values, deviations to mean 
  meandB<-mean(dB)
  devdB<-dB-meandB
  
  # computation of MDiv morphometric divergence
  MDiv<-round((sum(devdB)+meandB) / (sum(abs(devdB))+meandB) ,6)# arendre
  
  #computation of morphological Originality
  dist_sp<-as.matrix(dist(coord_d,method="euclidean")) ; dist_sp[which(dist_sp==0)]<-NA
  orig_sp<-apply(dist_sp, 1, min, na.rm=T)
  ls_original_sp <- names(which.min(orig_sp))
  mt_original_sp <- names(which.max(orig_sp))
  
  Orig_vect <- c(mean=mean(orig_sp,na.rm=T), min=min(orig_sp, na.rm=T), max=max(orig_sp, na.rm=T), se=sd(orig_sp)/sqrt(length(orig_sp)),ls_original_sp=ls_original_sp, mt_original_sp= mt_original_sp)
  
  #computing functional specialization 
  O<-apply(coord_d, 2, mean)
  speS<-apply(coord_d, 1, function(x) { (sum((x-O)^2) )^0.5} )
  Specialization <- c(mean=mean(speS,na.rm=T), min=min(speS, na.rm=T), max=max(speS, na.rm=T), se=sd(speS)/sqrt(length(speS)))
  
  list(data=data.frame(RS=rs, MRic=MRic,MDiv=MDiv,nb_vertices=nbME_vertices),specialization= Specialization, originality=Orig_vect)
  
} # 

#===================== Plot tangent space total ==============================#
plotTangentSpace_tot <- function (A,axis1=1,axis2=2,warpgrids=TRUE,mesh = NULL,  gp,
                                  label = NULL,groups=groups,legend=unique(dimnames(A)[[3]]), col_sp =unique(col.gp),
                                  pos1=c(0.11,0.28,0.60,0.84),s_class=FALSE,
                                  pos2=c(0.5,1,0,0.5),sizept_grid=0.8,cex.pt=0.6,type="n"
                                  ,xlim=c(-0.3,0.2),ylim=c(-0.18,0.2),Traits=NULL,txt.cex=1,
                                  col_Traits=c("#66C3E1","#294AE7"),polyall=TRUE,
                                  ax_calInd=c(1,3),cex.axis=0.2,lwd_axis=0.2,
                                  clabel=0.7,mgp_axisx=c(1.8,0.5,0),mgp_axisy=c(1.8,0.5,0),flip=F,...) {
  if (length(dim(A)) != 3) {
    stop("Data matrix not a 3D array (see 'arrayspecs').")
  }
  if (any(is.na(A)) == T) {
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').")
  }
  dots <- list(...)
  retx <- dots$retx
  if (is.null(retx)) 
    retx <- TRUE
  scale. <- dots$scale.
  if (is.null(scale.)) 
    scale. <- FALSE
  center <- dots$center
  if (is.null(center)) 
    center <- TRUE
  tol <- dots$tol
  k <- dim(A)[2]
  p <- dim(A)[1]
  n <- dim(A)[3]
  ref <- mshape(A)
  x <- two.d.array(A)
  if (is.null(tol)) {
    d <- prcomp(x)$sdev^2
    cd <- cumsum(d)/sum(d)
    cd <- length(which(cd < 1))
    if (length(cd) < length(d)) 
      cd <- cd + 1
    tol <- max(c(d[cd]/d[1], 0.005))
  }
  
  gp=dimnames(A)[[3]]
  pc.res <- prcomp(x, center = center, scale. = scale., retx = retx, tol = tol)
  pcdata <- pc.res$x
  if(flip==TRUE){pcdata[,2] <- (pcdata[, 2])*-1}
  a <- summary(pc.res)$importance
  eigenVal <- a[2,axis1:axis2]
  
  xlab1 <- paste("PC ",axis1," (",round(eigenVal[axis1]*100,2),"%)",sep="")
  ylab1 <- paste("PC ",axis2," (",round(eigenVal[axis2]*100,2),"%)",sep="")
  
  #### Morphological indicators calculation
  data_tot <- data.frame(x=pcdata[,axis1],y=pcdata[,axis2])
  vert0 <- convhulln(data_tot,"Fx TO 'vert.txt'")
  vert1 <- scan("vert.txt",quiet=T)
  vert2 <- (vert1+1)[-1]
  
  FE_vert_tot <- row.names(data_tot)[vert2]
  ResMorpsp <- gb_Mric_function(pcdata=pcdata,ax=ax_calInd)

  #############################################" 
  shapes <- shape.names <- NULL
  for (i in 1:ncol(pcdata)) {
    pcaxis.min <- min(pcdata[, i])
    pcaxis.max <- max(pcdata[, i])
    pc.min <- pc.max <- rep(0, dim(pcdata)[2])
    pc.min[i] <- pcaxis.min
    pc.max[i] <- pcaxis.max
    pc.min <- as.matrix(pc.min %*% (t(pc.res$rotation))) + as.vector(t(ref))
    pc.max <- as.matrix(pc.max %*% (t(pc.res$rotation))) + as.vector(t(ref))
    shapes <- rbind(shapes, pc.min, pc.max)
    shape.names <- c(shape.names, paste("PC", i, "min", sep = ""), paste("PC", i, "max", sep = ""))
  }
  shapes <- arrayspecs(shapes, p, k)
  shapes <- lapply(seq(dim(shapes)[3]), function(x) shapes[, , x])
  names(shapes) <- shape.names
  
  if (warpgrids == TRUE) {
    xlab1 <- paste("PC ",axis1," (",round(eigenVal[axis1]*100,2),"%)",sep="")
    ylab1 <- paste("PC ",axis2," (",round(eigenVal[axis2]*100,2),"%)",sep="")
    
    plot(pcdata[,axis1],pcdata[,axis2],xlim=xlim,ylim=ylim,asp=1,pch=21,bg="black",
         cex=cex.pt,xlab=xlab1,ylab=ylab1,type=type,yaxt="n",mgp=c(1.8,0.5,0),tcl=-0.25)
    axis(side=2,at=seq(-0.2,0.3,0.1),labels=seq(-0.2,0.3,0.1),las=2,tcl=-0.25,mgp=c(1.8,0.5,0), cex=1)
    
    if(warpgrids == TRUE && polyall == TRUE){
      polygon(data.frame(data_tot[FE_vert_tot,1],data_tot[FE_vert_tot,2]),col="#9A9B9F30")
      points(data_tot[FE_vert_tot,1],data_tot[FE_vert_tot,2],pch=21,bg="#9A9B9F",cex=cex.pt)
      
      labels_leg <- c(paste0("Vertice number = ",ResMorpsp$data$nb_vertices),
                      paste0("SR = ",ResMorpsp$data$RS), 
                      paste0("MRic = ",round(ResMorpsp$data$MRic,3)),
                      paste0("MDiv = ",round(ResMorpsp$data$MDiv,3)))
      #legend("topright",legend=labels_leg,title="Morphological index",cex=0.75,bg="#9A9B9F30")
      #legend("bottomleft", legend = unique(names(groups)), pch = 19, bty = "n", 
      #       col = unique(groups), cex=0.7,horiz=T)
    } # end of if
    
    if (!is.null(groups)) {
      reg <- unique(names(groups))
      
      MRic_Vect <- vector()
      for(i in 1:length(reg)){
        cat("i=",i,"\n")
        tokeep <- which(gp==reg[i])
        
        if(length(tokeep)<=2){
          
        } else {
          data <- data.frame(x=pcdata[tokeep, axis1], y=pcdata[tokeep, axis2])
          vert0<-convhulln(data,"Fx TO 'vert.txt'")
          vert1<-scan("vert.txt",quiet=T)
          vert2<-(vert1+1)[-1]
          
          MRic_Vect[i] <- round(convhulln(data,"FA")$vol,10)
          
          FE_vert_k<- row.names(data)[vert2]
          polygon(data.frame(data[FE_vert_k,1],data[FE_vert_k,2]),col=alpha(col_sp[i],0.4)) 
          points(data[FE_vert_k,1],data[FE_vert_k,2],pch=21,bg=1,cex=cex.pt)
        }
      }
      points(pcdata[, axis1], pcdata[, axis2], pch=19, bg=groups, col=groups, cex=0.4)
    }
    
    segments(min(pcdata[,axis1]),0,max(pcdata[,axis1]),0,lty=2,lwd=1)
    segments(0,min(pcdata[, axis2]),0,max(pcdata[, axis2]),lty=2,lwd=1)
    
    if(s_class==TRUE){
      col_class <- colorRampPalette(col_Traits)(length(levels(Traits)))
      s.class_mod(dfxy=data.frame(pcdata[,axis1], pcdata[,axis2]),fac=Traits, 
                  wt=rep(1,length(Traits)),xax=1,yax=2,cstar=0,cellipse=1.5,axesell=T, 
                  label=levels(Traits),clabel=0.7,cpoint=0,pch=20,col=col_class,xlim=xlim,
                  ylim=ylim,grid=F,addaxes=T,origin=c(0,0),cgrid=1.3,pixmap=NULL,
                  contour=NULL,area=NULL,add.plot=TRUE)
    } # end of if
    
    if (length(label != 0)) {
      if (isTRUE(label)) {text(pcdata[, axis1], pcdata[, axis2], seq(1,n), adj=c(-0.7, -0.7))}
      else {
        text(pcdata[, axis1], pcdata[, axis2],label, adj = c(-0.1, -0.1), cex=txt.cex)
      }
    }
    
    shape.min <- shapes[[which(names(shapes) == paste("PC", axis1, "min", sep = ""))]]
    shape.max <- shapes[[which(names(shapes) == paste("PC", axis1, "max", sep = ""))]]
    
    if (k == 2) {
      arrows(min(pcdata[, axis1]), (0.7 * max(pcdata[, axis2])), min(pcdata[, axis1]), 0, length = 0.1, lwd = 2)
      arrows(max(pcdata[, axis1]), (0.7 * min(pcdata[, axis2])), max(pcdata[, axis1]), 0, length = 0.1, lwd = 2)
      
      par(fig=pos1,new=TRUE,mar=c(5,0,0,0))
      geomorph:::tps(ref, shape.min, 20,sz = sizept_grid)
      par(fig=pos2,new=TRUE,mar=c(4,10,10,4))
      geomorph:::tps(ref, shape.max, 20,sz=sizept_grid)
    }
    
    if (!is.null(groups) && legend == TRUE) {
      par(new=T)
      if (is.factor(groups)) {
        legend(0.1, 0.1, legend = unique(groups), pch = 19, bty = "n", col = unique(groups))
      }
      else {
        legend(0.1,0.1, legend = unique(names(groups)), pch = 19, bty = "n", col = unique(groups))
      }
    }
    
    if (k == 3) {
      if (is.null(mesh) == TRUE) {
        open3d()
        mfrow3d(1, 2)
        plot3d(shape.min, type = "s", col = "gray", main = paste("PC ", axis1, " negative"), size = 1.25, aspect = FALSE, 
               xlab = "", ylab = "", zlab = "", box = FALSE, 
               axes = FALSE)
        plot3d(shape.max, type = "s", col = "gray", main = paste("PC ", axis1, " positive"), size = 1.25, aspect = FALSE, 
               xlab = "", ylab = "", zlab = "", box = FALSE, 
               axes = FALSE)
      }
      if (is.null(mesh) == FALSE) {
        open3d()
        mfrow3d(1, 2)
        cat(paste("\nWarping mesh to negative end of axis ", 
                  axis1, "\n", sep = ""))
        plotRefToTarget(ref, shape.min, mesh, method = "surface")
        title3d(main = paste("PC ", axis1, " negative"))
        next3d()
        cat(paste("\nWarping mesh to positive end of axis ", 
                  axis1, "\n", sep = ""))
        plotRefToTarget(ref, shape.max, mesh, method = "surface")
        title3d(main = paste("PC ", axis1, " positive"))
      }
    }
    layout(1)
  }
  
  
  if (warpgrids == FALSE) {
    
    xlab1 <- paste("PC ",axis1," (",round(eigenVal[axis1]*100,2),"%)",sep="")
    ylab1 <- paste("PC ",axis2," (",round(eigenVal[axis2]*100,2),"%)",sep="")
    
    plot(pcdata[,axis1],pcdata[,axis2],xlim=xlim,ylim=ylim,asp=1,pch=21,bg="black",
         cex=cex.pt,xlab=xlab1,ylab=ylab1,type=type,yaxt="n",mgp=c(1.8,0.5,0),tcl=-0.25,cex.axis=0.75,cex.lab=0.8)
    axis(side=2,at=seq(-0.2,0.3,0.1),labels=seq(-0.2,0.3,0.1),las=2,tcl=-0.25,mgp=c(1.8,0.5,0),cex.axis=0.75)
    box()
    
    if(warpgrids == FALSE && polyall == TRUE){
      
      polygon(data.frame(data_tot[FE_vert_tot,1],data_tot[FE_vert_tot,2]),col="#9A9B9F30",lwd=0.6)
      points(data_tot[FE_vert_tot,1],data_tot[FE_vert_tot,2],pch=21,bg="#9A9B9F",cex=cex.pt)
      
      labels_leg <- c(paste0("Vertice number = ",ResMorpsp$data$nb_vertices),
                      paste0("SR = ",ResMorpsp$data$RS), 
                      paste0("MRic = ",round(ResMorpsp$data$MRic,3)),
                      paste0("MDiv = ",round(ResMorpsp$data$MDiv,3)))
    } # end of if
    
    if (!is.null(groups)) {
      reg <- unique(names(groups))
      MRic_Vect <- vector()
      for(i in 1:length(reg)){
        cat("i=",i,"\n")
        tokeep <- which(gp==reg[i])
        
        if(length(tokeep)<=2){
          
        } else {
          data<- data.frame(x=pcdata[tokeep, axis1], y=pcdata[tokeep, axis2])
          vert0<-convhulln(data,"Fx TO 'vert.txt'")
          vert1<-scan("vert.txt",quiet=T)
          vert2<-(vert1+1)[-1]
          
          MRic_Vect[i] <- round(convhulln(data,"FA")$vol,10)
          FE_vert_k<- row.names(data)[vert2]
          polygon(data.frame(data[FE_vert_k,1],data[FE_vert_k,2]),col=paste0(unique(groups)[i],60),lwd=0.6)
          points(data[FE_vert_k,1],data[FE_vert_k,2],pch=21,bg=paste0(unique(groups)[i],60),cex=cex.pt)
        }
      }
      points(pcdata[, axis1], pcdata[, axis2], pch = 19, col = groups, cex = cex.pt,lwd=0.1)
    }
    
    segments(min(pcdata[,axis1]),0,max(pcdata[, axis1]),0,lty=2)
    segments(0, min(pcdata[,axis2]),0,max(pcdata[,axis2]),lty =2)
    
    if(s_class==TRUE){
      col_class <- colorRampPalette(col_Traits)(length(levels(Traits)))
      s.class_mod(dfxy=data.frame(pcdata[,axis1], pcdata[,axis2]),fac=Traits, 
                  wt=rep(1,length(Traits)),xax=1,yax=2,cstar=0,cellipse=1.5,axesell=T, 
                  label=levels(Traits),clabel=clabel,cpoint=0,pch=20,col=col_class,xlim=xlim,
                  ylim=ylim,grid=F,addaxes=T,origin=c(0,0),cgrid=1.3,pixmap=NULL,
                  contour=NULL,area=NULL,add.plot=TRUE)
    } # end of if
  } # end of warpgrids=False  
  
  out <- list(pc.summary = summary(pc.res), pc.scores = pcdata, 
              pc.shapes = shapes, sdev = pc.res$sdev, rotation = pc.res$rotation,MRic_Vect=(MRic_Vect/ResMorpsp$data$MRic))
  
  class(out) = "plotTangentSpace"
  out
} #  end of plotTangentSpace_tot 

########################################################
### get shapes
get_shapes  <- function (A,axis1=1,axis2=2,mesh = NULL,label = NULL,groups=groups,
                         legend=unique(dimnames(A)[[3]]),flip=F,...) {
  if (length(dim(A)) != 3) {
    stop("Data matrix not a 3D array (see 'arrayspecs').")
  }
  if (any(is.na(A)) == T) {
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').")
  }
  dots <- list(...)
  retx <- dots$retx
  if (is.null(retx)) 
    retx <- TRUE
  scale. <- dots$scale.
  if (is.null(scale.)) 
    scale. <- FALSE
  center <- dots$center
  if (is.null(center)) 
    center <- TRUE
  tol <- dots$tol
  k <- dim(A)[2]
  p <- dim(A)[1]
  n <- dim(A)[3]
  ref <- mshape(A)
  x <- two.d.array(A)
  if (is.null(tol)) {
    d <- prcomp(x)$sdev^2
    cd <- cumsum(d)/sum(d)
    cd <- length(which(cd < 1))
    if (length(cd) < length(d)) 
      cd <- cd + 1
    tol <- max(c(d[cd]/d[1], 0.005))
  }
  
  gp <- dimnames(A)[[3]]
  pc.res <- prcomp(x, center = center, scale. = scale., retx = retx, tol = tol)
  pcdata <- pc.res$x
  if(flip==TRUE){
    pcdata[,2] <- (pcdata[, 2])*-1
  }
  a <- summary(pc.res)$importance
  eigenVal <-a[2,axis1:axis2]
  
  #############################################" 
  shapes <- shape.names <- NULL
  for (i in 1:ncol(pcdata)) {
    pcaxis.min <- min(pcdata[, i])
    pcaxis.max <- max(pcdata[, i])
    pc.min <- pc.max <- rep(0, dim(pcdata)[2])
    pc.min[i] <- pcaxis.min
    pc.max[i] <- pcaxis.max
    pc.min <- as.matrix(pc.min %*% (t(pc.res$rotation))) + as.vector(t(ref))
    pc.max <- as.matrix(pc.max %*% (t(pc.res$rotation))) + as.vector(t(ref))
    shapes <- rbind(shapes, pc.min, pc.max)
    shape.names <- c(shape.names, paste("PC", i, "min", sep = ""), paste("PC", i, "max", sep = ""))
  }
  shapes <- arrayspecs(shapes, p, k)
  shapes <- lapply(seq(dim(shapes)[3]), function(x) shapes[, , x])
  names(shapes) <- shape.names
  
  shape.min <- shapes[[which(names(shapes) == paste("PC", axis1, "min", sep = ""))]]
  shape.max <- shapes[[which(names(shapes) == paste("PC", axis1, "max", sep = ""))]]
  
  list(shape.min=shape.min,shape.max=shape.max,pcdata=pcdata,ref=ref)
  
}
#######################################################

plotTangentSpace1 <- function (A, axis1 = 1, axis2 = 2, warpgrids = TRUE, mesh = NULL, 
                               label = NULL, groups = NULL, legend = FALSE,gp=gp, ...,pos1=c(0.11,0.28,0.60,0.84),pos2=c(0.5,1,0,0.5),sizept_grid=1.2) 
{
  if (length(dim(A)) != 3) {
    stop("Data matrix not a 3D array (see 'arrayspecs').")
  }
  if (any(is.na(A)) == T) {
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').")
  }
  dots <- list(...)
  retx <- dots$retx
  if (is.null(retx)) 
    retx <- TRUE
  scale. <- dots$scale.
  if (is.null(scale.)) 
    scale. <- FALSE
  center <- dots$center
  if (is.null(center)) 
    center <- TRUE
  tol <- dots$tol
  k <- dim(A)[2]
  p <- dim(A)[1]
  n <- dim(A)[3]
  ref <- mshape(A)
  x <- two.d.array(A)
  gp=dimnames(A)[[3]]
  
  if (is.null(tol)) {
    d <- prcomp(x)$sdev^2
    cd <- cumsum(d)/sum(d)
    cd <- length(which(cd < 1))
    if (length(cd) < length(d)) 
      cd <- cd + 1
    tol <- max(c(d[cd]/d[1], 0.005))
  }
  pc.res <- prcomp(x, center = center, scale. = scale., retx = retx, tol = tol)
  pcdata <- pc.res$x
  if (warpgrids == FALSE) {
    if (legend == TRUE) {
      layout(t(matrix(c(1, 1, 2, 1, 1, 1, 1, 1, 1), 3, 3)))
    }
    plot(pcdata[, axis1], pcdata[, axis2], asp = 1, pch = 21, 
         bg = "black", cex = 2, xlab = paste("PC ", axis1), 
         ylab = paste("PC ", axis2),...)
    if (!is.null(groups)) {
      points(pcdata[, axis1], pcdata[, axis2], pch = 21, 
             bg = groups, cex = 2)
    }
    segments(min(pcdata[, axis1]), 0, max(pcdata[, axis1]), 
             0, lty = 2, lwd = 1)
    segments(0, min(pcdata[, axis2]), 0, max(pcdata[, axis2]), 
             lty = 2, lwd = 1)
    if (length(label != 0)) {
      if (isTRUE(label)) {
        text(pcdata[, axis1], pcdata[, axis2], seq(1, 
                                                   n), adj = c(-0.7, -0.7))
      }
      else {
        text(pcdata[, axis1], pcdata[, axis2], label, 
             adj = c(-0.1, -0.7))
      }
    }
    if (!is.null(groups) && legend == TRUE) {
      plot.new()
      if (is.factor(groups)) {
        legend(0.5, 1, legend = unique(groups), pch = 19, 
               bty = "n", col = unique(groups))
      }
      else {
        legend(0.5, 1, legend = unique(names(groups)), 
               pch = 19, bty = "n", col = unique(groups))
      }
    }
  }
  
  shapes <- shape.names <- NULL
  for (i in 1:ncol(pcdata)) {
    pcaxis.min <- min(pcdata[, i])
    pcaxis.max <- max(pcdata[, i])
    pc.min <- pc.max <- rep(0, dim(pcdata)[2])
    pc.min[i] <- pcaxis.min
    pc.max[i] <- pcaxis.max
    pc.min <- as.matrix(pc.min %*% (t(pc.res$rotation))) + 
      as.vector(t(ref))
    pc.max <- as.matrix(pc.max %*% (t(pc.res$rotation))) + 
      as.vector(t(ref))
    shapes <- rbind(shapes, pc.min, pc.max)
    shape.names <- c(shape.names, paste("PC", i, "min", sep = ""), paste("PC", i, "max", sep = ""))
  }
  shapes <- arrayspecs(shapes, p, k)
  shapes <- lapply(seq(dim(shapes)[3]), function(x) shapes[, , x])
  names(shapes) <- shape.names
  
  if (warpgrids == TRUE) {
    plot(pcdata[, axis1], pcdata[, axis2], asp = 1, pch = 21, bg = "black", cex = 1, xlab = paste("PC ", axis1), ylab = paste("PC ", axis2))
    if (!is.null(groups)) {
      
      reg <- unique(names(groups))
      
      for(i in 1:length(reg)){
        cat("i=",i,"\n")
        tokeep <- which(gp==reg[i])
        
        if(length(tokeep)<=2){
          
        } else {
          data<- data.frame(x=pcdata[tokeep, axis1], y=pcdata[tokeep, axis2])
          vert0<-convhulln(data,"Fx TO 'vert.txt'")
          vert1<-scan("vert.txt",quiet=T)
          vert2<-(vert1+1)[-1]
          FE_vert_k<- row.names(data)[vert2]
          
          polygon(data.frame(data[FE_vert_k,1],data[FE_vert_k,2]),col=paste0(unique(groups)[i],60))
        }
      }
      points(pcdata[, axis1], pcdata[, axis2], pch = 21, bg = groups, cex = 1)
    }
    
    segments(min(pcdata[, axis1]), 0, max(pcdata[, axis1]),0, lty = 2, lwd = 1)
    segments(0,min(pcdata[, axis2]), 0, max(pcdata[, axis2]), lty = 2, lwd = 1)
    
    if (length(label != 0)) {
      if (isTRUE(label)) {text(pcdata[, axis1], pcdata[, axis2], seq(1,n), adj = c(-0.7, -0.7))}
      else {
        text(pcdata[, axis1], pcdata[, axis2], label, adj = c(-0.1, -0.1))
      }
    }
    shape.min <- shapes[[which(names(shapes) == paste("PC", axis1, "min", sep = ""))]]
    shape.max <- shapes[[which(names(shapes) == paste("PC", axis1, "max", sep = ""))]]
    if (k == 2) {
      arrows(min(pcdata[, axis1]), (0.7 * max(pcdata[, axis2])), min(pcdata[, axis1]), 0, length = 0.1, lwd = 2)
      arrows(max(pcdata[, axis1]), (0.7 * min(pcdata[, axis2])), max(pcdata[, axis1]), 0, length = 0.1, lwd = 2)
      
      par(fig=pos1,new=TRUE,mar=c(0,0,0,0))
      geomorph:::tps(ref, shape.min, 20,sz = sizept_grid)
      par(fig=pos2,new=TRUE,mar=c(0,0,0,0))
      geomorph:::tps(ref, shape.max, 20,sz=sizept_grid)
    }
    if (!is.null(groups) && legend == TRUE) {
      plot.new()
      if (is.factor(groups)) {
        legend(0.5, 1, legend = unique(groups), pch = 19, bty = "n", col = unique(groups))
      }
      else {
        legend(0.5, 1, legend = unique(names(groups)), pch = 19, bty = "n", col = unique(groups))
      }
    }
    if (k == 3) {
      if (is.null(mesh) == TRUE) {
        open3d()
        mfrow3d(1, 2)
        plot3d(shape.min, type = "s", col = "gray", main = paste("PC ", 
                                                                 axis1, " negative"), size = 1.25, aspect = FALSE, 
               xlab = "", ylab = "", zlab = "", box = FALSE, 
               axes = FALSE)
        plot3d(shape.max, type = "s", col = "gray", main = paste("PC ", 
                                                                 axis1, " positive"), size = 1.25, aspect = FALSE, 
               xlab = "", ylab = "", zlab = "", box = FALSE, 
               axes = FALSE)
      }
      if (is.null(mesh) == FALSE) {
        open3d()
        mfrow3d(1, 2)
        cat(paste("\nWarping mesh to negative end of axis ", 
                  axis1, "\n", sep = ""))
        plotRefToTarget(ref, shape.min, mesh, method = "surface")
        title3d(main = paste("PC ", axis1, " negative"))
        next3d()
        cat(paste("\nWarping mesh to positive end of axis ", 
                  axis1, "\n", sep = ""))
        plotRefToTarget(ref, shape.max, mesh, method = "surface")
        title3d(main = paste("PC ", axis1, " positive"))
      }
    }
    layout(1)
  }
  out <- list(pc.summary = summary(pc.res), pc.scores = pcdata, 
              pc.shapes = shapes, sdev = pc.res$sdev, rotation = pc.res$rotation)
  class(out) = "plotTangentSpace"
  out
}
