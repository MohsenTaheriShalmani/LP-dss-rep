library("numDeriv") #find normal vectors

normalsOfSkeletalSheet <- function(centeredSkel) {
  
  # x y z
  x<-centeredSkel[,1]
  y<-centeredSkel[,2]
  z<-centeredSkel[,3]
  
  fit2 <- lm(z ~ poly(x, y, degree = 4 ,raw = TRUE), data=as.data.frame(cbind(z,x,y)))
  # summary(fit2)
  
  newDATA<-data.frame(x=x, y=y)
  surfacePointsZ<-predict(fit2, newdata = newDATA)

  medialPoints<-cbind(x,y,surfacePointsZ)
  
  #goodness of fit
  euclideanDistances<-array(NA,dim(centeredSkel)[1])
  for (i in 1:dim(centeredSkel)[1]) {
    euclideanDistances[i]<-norm(medialPoints[i,]-centeredSkel[i,],type = "2")
  }
  # sum(euclideanDistances)   #L1 norm
  # sum(euclideanDistances^2) #L2 norm
  cat("Goodness of fit is:",sum(euclideanDistances^2)/centroid.size(centeredSkel),"\n") #difference is too little
  
  # normals 
  coefs<-fit2$coefficients
  
  fitedFuction <- function(xy) {          #polynomial degree 4
    coefs[[1]]*(xy[1]^0)*(xy[2]^0)+
      coefs[[2]]*(xy[1]^1)*(xy[2]^0)+
      coefs[[3]]*(xy[1]^2)*(xy[2]^0)+
      coefs[[4]]*(xy[1]^3)*(xy[2]^0)+
      coefs[[5]]*(xy[1]^4)*(xy[2]^0)+
      coefs[[6]]*(xy[1]^0)*(xy[2]^1)+
      coefs[[7]]*(xy[1]^1)*(xy[2]^1)+
      coefs[[8]]*(xy[1]^2)*(xy[2]^1)+
      coefs[[9]]*(xy[1]^3)*(xy[2]^1)+
      coefs[[10]]*(xy[1]^0)*(xy[2]^2)+
      coefs[[11]]*(xy[1]^1)*(xy[2]^2)+
      coefs[[12]]*(xy[1]^2)*(xy[2]^2)+
      coefs[[13]]*(xy[1]^0)*(xy[2]^3)+
      coefs[[14]]*(xy[1]^1)*(xy[2]^3)+
      coefs[[15]]*(xy[1]^0)*(xy[2]^4)-
      xy[3]
  }
  
  medialVectors<-array(NA,dim = dim(medialPoints))
  medialVectors2<-array(NA,dim = dim(medialPoints)) #this is an estimation of normals of the lm surface
  medialNormals<-array(NA,dim = dim(medialPoints))
  for (i in 1:dim(medialPoints)[1]) {
    medialNormals[i,]<-convertVec2unitVec((-grad(fitedFuction, medialPoints[i,])))
    medialVectors[i,]<-medialNormals[i,]+medialPoints[i,]
    medialVectors2[i,]<-medialNormals[i,]+centeredSkel[i,]
  }
  
  result<-list(medialNormals=medialNormals, medialVectors=medialVectors, medialVectors2=medialVectors2)
  return(result)
}
