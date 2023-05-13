normalsOfSkeletalSheetByConnections <- function(skelPoints,
                                                framesCenters,
                                                backPointsIndices=backPointsIndices,
                                                frontPointsIndices=frontPointsIndices) {
  
  medialNormals<-array(NA,dim = dim(skelPoints))
  for (i in 1:length(framesCenters)) {
    if(is.na(framesCenters[i])){
      next
    }
    
    a<-backPointsIndices[i]
    b<-framesCenters[i]
    c<-frontPointsIndices[i]
    
    p1<-skelPoints[a,]
    p2<-skelPoints[b,]
    p3<-skelPoints[c,]
    
    v1<-convertVec2unitVec2(p1-p2)
    v2<-convertVec2unitVec2(p3-p2)
    
    u1<-convertVec2unitVec2(v2-v1)
    
    u2<-convertVec2unitVec2(myCrossProduct(v2,v1))
    
    normalTemp1<-convertVec2unitVec2(myCrossProduct(u2,u1))
    normalTemp2<-convertVec2unitVec2(myCrossProduct(u1,u2))
    
    if(normalTemp1[3]>normalTemp2[3]){
      medialNormals[i,]<-normalTemp1
    }else{
      medialNormals[i,]<-normalTemp2
    }
    
  }
  
  return(medialNormals)
  
}
