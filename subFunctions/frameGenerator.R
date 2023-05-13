
secondFrameVectorFunction <- function(p1,vertex,p2,normalVec) {
  
  u1<-convertVec2unitVec(p2-vertex)
  u2<-convertVec2unitVec(p1-vertex)
  
  v<-convertVec2unitVec(u1-u2)
  
  #projection of tempPoint on tangent space
  #sum(v*normalVec) is dot product that is the point distance from the plane
  projected_point <- v-(sum(v*normalVec)/norm(normalVec,type = "2"))*normalVec
  
  secondVec<-convertVec2unitVec(projected_point)
  
  return(secondVec)
}

# generate the frames
frameGenerator <- function(centeredSkel,medialNormals,
                           framesCenters,framesBackPoints,framesFronts) {
  
  numberOfFrames<-length(framesCenters)
  
  # NB!!! number of frames is equal to upSpoeksNumber
  framesFirstVec<-array(NA,dim = c(numberOfFrames,3))
  framesSecondVec<-array(NA,dim = c(numberOfFrames,3))
  framesThirdVec<-array(NA,dim = c(numberOfFrames,3))
  for (i in 1:numberOfFrames) {
    
    a<-framesBackPoints[i]
    b<-framesCenters[i]
    c<-framesFronts[i]
    
    framesFirstVec[b,]<-medialNormals[b,]
    
    if(c==Inf){
      framesSecondVec[b,]<-secondFrameVectorFunction(p1 = centeredSkel[a,],vertex = centeredSkel[b,],
                                                     p2 =convertVec2unitVec(centeredSkel[b,]-centeredSkel[a,])+centeredSkel[b,],
                                                     normalVec = medialNormals[b,])
    }else{
      framesSecondVec[b,]<-secondFrameVectorFunction(p1 = centeredSkel[a,],vertex = centeredSkel[b,],
                                                     p2 = centeredSkel[c,],normalVec = medialNormals[b,])
    }
    
    framesThirdVec[b,]<-convertVec2unitVec(myCrossProduct(framesFirstVec[b,],framesSecondVec[b,]))
  }
  
  result<-list(framesFirstVec=framesFirstVec,
               framesSecondVec=framesSecondVec,
               framesThirdVec=framesThirdVec)
  
  return(result)
}



source("subFunctions/normalsOfSkeletalSheetBySpline.R")

frameGenerator2 <- function(centeredSkel,
                            framesCenters,
                            framesBackPoints,
                            framesFronts,
                            numberOfFrames,
                            numberOfLayers,
                            numberOf2DspokePoints) {
  
  # normal vectors
  temp<-normalsOfSkeletalSheet(centeredSkel = centeredSkel)
  medialNormals<-temp$medialNormals
  
  framesFirstVec<-array(NA,dim = c(numberOfFrames,3))
  framesSecondVec<-array(NA,dim = c(numberOfFrames,3))
  framesThirdVec<-array(NA,dim = c(numberOfFrames,3))
  
  for (i in 1:numberOfFrames) {
    
    a<-framesBackPoints[i]
    b<-framesCenters[i]
    c<-framesFronts[i]
    
    framesFirstVec[b,]<-medialNormals[b,]
    
    if(b%%numberOfLayers>(numberOf2DspokePoints-1) | b%%numberOfLayers==0){
      framesSecondVec[b,]<-secondFrameVectorFunction(p1 = centeredSkel[a,],
                                                     vertex = centeredSkel[b,],
                                                     p2 = centeredSkel[c,],
                                                     normalVec = medialNormals[b,]) 
    }else{
      framesSecondVec[b,]<-secondFrameVectorFunction(p1 = centeredSkel[c,],
                                                     vertex = centeredSkel[b,],
                                                     p2 = centeredSkel[a,],
                                                     normalVec = medialNormals[b,])
    }
    
    framesThirdVec[b,]<-convertVec2unitVec(myCrossProduct(framesFirstVec[b,],framesSecondVec[b,]))
  }
  
  # frames
  result<-list(framesFirstVec=framesFirstVec,
               framesSecondVec=framesSecondVec,
               framesThirdVec=framesThirdVec)
  
  return(result)
  
}


frameGenerator4 <- function(centeredSkel,
                            framesCenters,
                            framesBackPoints,
                            framesFronts,
                            numberOfFrames,
                            numberOfLayers,
                            numberOf2DspokePoints,
                            framesCentersCrestRemoved) {
  
  # normal vectors
  temp<-normalsOfSkeletalSheet(centeredSkel = centeredSkel)
  medialNormals<-temp$medialNormals
  
  framesFirstVec<-array(NA,dim = c(numberOfFrames,3))
  framesSecondVec<-array(NA,dim = c(numberOfFrames,3))
  framesThirdVec<-array(NA,dim = c(numberOfFrames,3))
  
  for (i in 1:numberOfFrames) {
    
    a<-framesBackPoints[i]
    b<-framesCenters[i]
    c<-framesFronts[i]
    
    framesFirstVec[b,]<-medialNormals[b,]
    
    if(!b %in% framesCentersCrestRemoved){
      framesSecondVec[b,]<-c(0,0,0)
      framesThirdVec[b,]<-c(0,0,0)
    }else{
      if(b%%numberOfLayers>(numberOf2DspokePoints-1) | b%%numberOfLayers==0){
        framesSecondVec[b,]<-secondFrameVectorFunction(p1 = centeredSkel[a,],
                                                       vertex = centeredSkel[b,],
                                                       p2 = centeredSkel[c,],
                                                       normalVec = medialNormals[b,]) 
      }else{
        framesSecondVec[b,]<-secondFrameVectorFunction(p1 = centeredSkel[c,],
                                                       vertex = centeredSkel[b,],
                                                       p2 = centeredSkel[a,],
                                                       normalVec = medialNormals[b,])
      }
      
      framesThirdVec[b,]<-convertVec2unitVec(myCrossProduct(framesFirstVec[b,],framesSecondVec[b,]))
    }  
  }
  
  
  # frames
  result<-list(framesFirstVec=framesFirstVec,
               framesSecondVec=framesSecondVec,
               framesThirdVec=framesThirdVec)
  
  return(result)
  
}

# frame generator without spline fit
frameGenerator5 <- function(centeredSkel,
                            medialNormals,
                            framesCenters,
                            framesBackPoints,
                            framesFronts,
                            numberOfFrames,
                            numberOfLayers,
                            numberOf2DspokePoints,
                            framesCentersCrestRemoved) {
  
  
  framesFirstVec<-array(NA,dim = c(numberOfFrames,3))
  framesSecondVec<-array(NA,dim = c(numberOfFrames,3))
  framesThirdVec<-array(NA,dim = c(numberOfFrames,3))
  
  for (i in 1:numberOfFrames) {
    
    a<-framesBackPoints[i]
    b<-framesCenters[i]
    c<-framesFronts[i]
    
    framesFirstVec[b,]<-medialNormals[b,]
    
    if(!b %in% framesCentersCrestRemoved){
      framesSecondVec[b,]<-c(0,0,0)
      framesThirdVec[b,]<-c(0,0,0)
    }else{
      if(b%%numberOfLayers>(numberOf2DspokePoints-1) | b%%numberOfLayers==0){
        framesSecondVec[b,]<-secondFrameVectorFunction(p1 = centeredSkel[a,],
                                                       vertex = centeredSkel[b,],
                                                       p2 = centeredSkel[c,],
                                                       normalVec = medialNormals[b,]) 
      }else{
        framesSecondVec[b,]<-secondFrameVectorFunction(p1 = centeredSkel[c,],
                                                       vertex = centeredSkel[b,],
                                                       p2 = centeredSkel[a,],
                                                       normalVec = medialNormals[b,])
      }
      
      framesThirdVec[b,]<-convertVec2unitVec(myCrossProduct(framesFirstVec[b,],framesSecondVec[b,]))
    }  
  }
  
  
  # frames
  result<-list(framesFirstVec=framesFirstVec,
               framesSecondVec=framesSecondVec,
               framesThirdVec=framesThirdVec)
  
  return(result)
  
}


# frame generator without spline fit
frameGenerator7 <- function(skelPoints,
                            medialNormals,
                            framesCenters,
                            framesBackPoints,
                            framesFrontPoints) {
  
  numberOfFrames<-nrow(skelPoints)
  
  framesFirstVec<-array(NA,dim = c(numberOfFrames,3))
  framesSecondVec<-array(NA,dim = c(numberOfFrames,3))
  framesThirdVec<-array(NA,dim = c(numberOfFrames,3))
  
  for (i in 1:numberOfFrames) {
    
    if(is.na(framesCenters[i])){
      next
    }
    
    a<-framesBackPoints[i]
    b<-framesCenters[i]
    c<-framesFrontPoints[i]
    
    framesFirstVec[i,]<-medialNormals[b,]
    
    framesSecondVec[i,]<-secondFrameVectorFunction(p1 = skelPoints[a,],
                                                   vertex = skelPoints[b,],
                                                   p2 = skelPoints[c,],
                                                   normalVec = medialNormals[b,]) 
    
    framesThirdVec[i,]<-convertVec2unitVec(myCrossProduct(framesFirstVec[i,],framesSecondVec[i,]))
    
  }
  
  
  # frames
  result<-list(framesFirstVec=framesFirstVec,
               framesSecondVec=framesSecondVec,
               framesThirdVec=framesThirdVec)
  
  return(result)
  
}




