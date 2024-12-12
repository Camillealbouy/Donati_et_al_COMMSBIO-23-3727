### function 2 load
halfCircle <- function(x,y,r,start=-pi/2,end=pi/2,nsteps=300,aspect_coef=3.6,...){
  rs <- seq(start,end,len=nsteps)
  
  if (start<0) {
    xc <- (x-r*cos(rs)/aspect_coef)
    yc <- y-r*sin(rs)
  } else{
    xc <- (x+r*cos(rs)/aspect_coef)
    yc <- y+r*sin(rs)
  }
  polygon(xc,yc,...)
}