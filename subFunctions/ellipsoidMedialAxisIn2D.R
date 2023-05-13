
#generate ellipsoid PDM without center
ellipsoidGenerator_2D_4 <- function(center,a,b,n){
  theta<-seq(0, 2*pi, length.out = n)
  # points<-center
  points<-c()
  for (i in 1:length(theta)) {
      x<-a*cos(theta[i])
      y<-b*sin(theta[i])
      points<-rbind(points,center+c(x,y))
  }
  return(points)
}


#number of boundary points
n<-41

myEllipse<-ellipsoidGenerator_2D_4(center = c(0,0),a = a,b = b,n = n)

a<-10
b<-5
medialAxis<-c()
theta<-seq(from=0,to=2*pi,length.out = n)
for (i in 1:n) {
  medialAxis<-rbind(medialAxis,c(((a^2-b^2)/a)*cos(theta[i]),0))
}

plot(myEllipse,pch=20,xlim = c(-a-1,a+1),ylim = c(-a-1,a+1),xlab = '',ylab = '')
par(new=TRUE)
plot(medialAxis,pch=20,xlim = c(-a-1,a+1),ylim = c(-a-1,a+1),xlab = '',ylab = '')
par(new=TRUE)
for (i in 1:n) {
  par(new=TRUE)
  plot(rbind(myEllipse[i,],medialAxis[i,]),type = 'l',xlim = c(-a-1,a+1),ylim = c(-a-1,a+1),xlab = '',ylab = '')
}
