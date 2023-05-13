library(RiemBase)

frechetMean <- function(directions) {
  
  allDirTemp<-t(directions)
  data1 <- list()
  for (j in 1:dim(allDirTemp)[1]){
    data1[[j]] <-allDirTemp[j,]
  }
  data2 <- riemfactory(data1, name="sphere")
  ### Compute Fre'chet Mean
  out1<- rbase.mean(data2)
  meanFrechet<-as.vector(out1$x)
  
  return(meanFrechet)
  
}

#calculate initial mean frames where frames is a 3*3*nsample matrix.
initialMeanFrames <- function(frames, method= "Frechet") {
  
  if(method != "Frechet" & method != "PNS"){
    stop("Method should be specified as Frechet or PNS ! ")
  }
  
  # NB! for extremely concentrated data we use Mardia mean direction 
  sd1V1<-prcomp(t(frames[1,,]))$sdev[1]
  sd2V1<-prcomp(t(frames[1,,]))$sdev[2]
  sd1V2<-prcomp(t(frames[2,,]))$sdev[1]
  sd2V2<-prcomp(t(frames[2,,]))$sdev[2]
  
  if(sd1V1<1e-02 | sd2V1<1e-02){
    v1<-convertVec2unitVec(colMeans(t(frames[1,,])))
  }else if(method == "Frechet"){
    v1<-frechetMean(frames[1,,])
  }else if(method == "PNS"){
    sphereType<-kurtosisTestFunction(frames[1,,])
    v1<-pns(frames[1,,],sphere.type = sphereType)$PNS$mean
  }
  
  if(sd1V2<1e-02 | sd2V2<1e-02){
    v2<-convertVec2unitVec(colMeans(t(frames[2,,])))
  }else if(method == "Frechet"){
    v2<-frechetMean(frames[2,,])
  }else if(method == "PNS"){
    sphereType<-kurtosisTestFunction(frames[2,,])
    v2<-pns(frames[2,,],sphere.type = sphereType)$PNS$mean
  }
  
  
  #NB! tempV3 is not v3
  tempV3<-convertVec2unitVec(myCrossProduct(v1,v2))
  meanVec<-convertVec2unitVec((v1+v2)/2)
  
  #rotate tempV3 to the north and meanVec to (0,0,1)
  twoVec<-rbind(tempV3,meanVec)
  R1<-rotMat(twoVec[1,],c(0,0,1))
  twoVec2<-twoVec%*%t(R1)
  R2<-rotMat(twoVec2[2,],c(1,0,0))
  twoVec3<-twoVec2%*%t(R2)
  #rotate back 
  R_back<-solve(t(R2%*%R1))
  m1<-c(cos(pi/4),-sin(pi/4),0)
  m2<-c(cos(pi/4),sin(pi/4),0)
  meanV1<-m1%*%R_back
  meanV2<-m2%*%R_back
  meanV3<-convertVec2unitVec(myCrossProduct(meanV1,meanV2))
  
  initMeanFrame<-rbind(meanV1,meanV2,meanV3)
  
  return(initMeanFrame)
}


gradientDescent4meanFrame <- function(frames,
                                      Boost=TRUE,
                                      stepSize=0.001,
                                      threshold=1e-2,
                                      whileLimit=1000,
                                      method= "Frechet",
                                      plotting=FALSE) {
  
  
  perfectFrame<-rbind(c(0,0,1),c(1,0,0),c(0,1,0))
  mean_n<-frechetMean(frames[1,,])
  mean_b1<-frechetMean(frames[2,,])
  mean_b2<-frechetMean(frames[3,,])
  
  if(Boost==FALSE){
    
    #Starting point as aligned centroids
    R1<-rotMat(convertVec2unitVec(c(1,1,1)),convertVec2unitVec(mean_n+mean_b1+mean_b2))
    perfectFrame<-perfectFrame%*%t(R1)
    
  }else if(Boost==TRUE){
    
    #Boost by starting point as mean of n and b
    perfectFrame<-initialMeanFrames(frames, method= method)
    
  }else{
    stop("Please specify Boost as TRUE or FALSE!")
  }
  
  d<-1
  k<-1
  # open3d()
  if(plotting==TRUE){
    spheres3d(x = 0, y = 0, z = 0, radius = 1,col = "lightblue", alpha=0.1)
    vectors3d(1.5*diag(3), color="black", lwd=1,alpha=0.5)
    vectors3d(mean_n, color="blue", lwd=1)
    vectors3d(mean_b1, color="red", lwd=1)
    vectors3d(mean_b2, color="green", lwd=1) 
  }
  while (d > threshold & k < whileLimit) { # stop while loop after whileLimit iteration 
    step1<-convertVec2unitVec(perfectFrame[1,]+
                                stepSize*convertVec2unitVec(mean_n-perfectFrame[1,]))
    R1<-rotMat(perfectFrame[1,],step1)
    perfectFrame<-perfectFrame%*%t(R1)
    step2<-convertVec2unitVec(perfectFrame[2,]+
                                stepSize*convertVec2unitVec(mean_b1-perfectFrame[2,]))
    R2<-rotMat(perfectFrame[2,],step2)
    perfectFrame<-perfectFrame%*%t(R2)
    step3<-convertVec2unitVec(perfectFrame[3,]+
                                stepSize*convertVec2unitVec(mean_b2-perfectFrame[3,]))
    R3<-rotMat(perfectFrame[3,],step3)
    perfectFrame<-perfectFrame%*%t(R3)
    a<-acos(pmin(pmax((sum(perfectFrame[1,]*mean_n)) ,-1.0), 1.0))
    b<-acos(pmin(pmax((sum(perfectFrame[2,]*mean_b1)) ,-1.0), 1.0))
    c<-acos(pmin(pmax((sum(perfectFrame[3,]*mean_b2)) ,-1.0), 1.0))
    d<-sqrt(a^2+b^2+c^2)
    cat("delta is:", d, "\n")
    k<-k+1
    
    if(plotting==TRUE){
      vectors3d(perfectFrame[1,], color="blue", lwd=1)
      vectors3d(perfectFrame[2,], color="red", lwd=1)
      vectors3d(perfectFrame[3,], color="green", lwd=1)
    }
  }
  
  # rgl.close()
  return(perfectFrame)
}


