### function 2 load
get_violplot <- function(sp="Ctenochaetus_striatus",ratio="ratio_ELO",data=Dat2pLot){
  y <- data[sp][[1]]
  Values <- y[,ratio]
  obs <- Values[c(which(y$category=="obs"))]
  Values <- Values[-c(which(y$category=="obs"))]
  return(list(obs=obs,Values=Values))
}

get_violplotInter <- function(Fam="Acanthuridae",ratio="ratio_ELO",data=Dat2pLot){
  y <- data[Fam][[1]]
  Values <- y[,ratio]
  obs <- Values[c(which(y$category=="obs"))]
  Values <- Values[-c(which(y$category=="obs"))]
  return(list(obs=obs,Values=Values))
}

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
