cutAndStretchSpokes <- function(BoundaryPDM,
                                SkeletalPDM,
                                meshPDM,
                                meshPolygonData) {
  numberOfPoints<-dim(meshPDM)[1]
  tipOfCuttedSpokes<-array(NA, dim = dim(BoundaryPDM))
  nTotalRadii<-dim(BoundaryPDM)[1]
  for (i in 1:nTotalRadii) {
    spokeDirectionGlobal<-convertVec2unitVec(BoundaryPDM[i,]-SkeletalPDM[i,])
    spokeTail<-SkeletalPDM[i,]
    intersections<-array(NA,dim = c(dim(meshPolygonData)[1],3))
    for (j in 1:dim(meshPolygonData)[1]) {
      p1<-polyMatrix[j,1]
      p2<-polyMatrix[j,2]
      p3<-polyMatrix[j,3]
      point1<-meshPDM[p1,]
      point2<-meshPDM[p2,]
      point3<-meshPDM[p3,]
      tempIntersection<-rayTriangleIntersection(rayOrigin = spokeTail,rayDirection = spokeDirectionGlobal,
                                                triangleVertex1 = point1,triangleVertex2 = point2,triangleVertex3 = point3)
      intersections[j,]<-tempIntersection
    }
    distances<-rep(Inf,dim(meshPolygonData)[1])
    for (k in 1:dim(meshPolygonData)[1]) {
      if(!is.na(intersections[k,1])){
        distances[k]<-norm(intersections[k,]-spokeTail,type = "2")
      }
    }
    
    tipOfCuttedSpokes[i,]<-intersections[which.min(distances),]
  }
  return(tipOfCuttedSpokes) #new boundary
}




cutAndStretchSpokes2D <- function(allSpokesTips,
                                  allSpokesTails,
                                  boundaryPoints1,
                                  boundaryPoints2) {
  
  # rayTriangleIntersection is in 3D so for 2D we consider the second point of the
  # triangle as c(colMeans(allSpokesTips),100) which is the vertex of the pyramid
  
  tipOfRays<-array(NA,dim = dim(allSpokesTips))
  for (i in 1:dim(allSpokesTips)[1]) {
    rayOrigin<-allSpokesTails[i,]
    rayDirection<-allSpokesTips[i,]-allSpokesTails[i,]
    intersections<-array(NA,dim = c(dim(boundaryPoints1)[1],3))
    p2<-c(colMeans(allSpokesTails),100) #p2 is fixed at the tip of pyramid
    for (j in 1:dim(boundaryPoints1)[1]) {
      p1<-c(boundaryPoints1[j,],0) #convert to 3D
      p3<-c(boundaryPoints2[j,],0) #convert to 3D
      tempIntersection<-rayTriangleIntersection(rayOrigin = c(rayOrigin,0), #convert to 3D
                                                rayDirection = c(rayDirection,0), #convert to 3D
                                                triangleVertex1 = p1,
                                                triangleVertex2 = p2,
                                                triangleVertex3 = p3) 
      
      intersections[j,]<-tempIntersection
    }
    if(sum(is.na(intersections[,1]))==dim(boundaryPoints1)[1]){
      tipOfRays[i,]<-allSpokesTips[i,]
    }else{
      distances<-rep(Inf,dim(boundaryPoints1)[1]) #find the closest intersect
      for (k in 1:dim(boundaryPoints1)[1]) {
        if(!is.na(intersections[k,1])){
          distances[k]<-norm(intersections[k,]-c(rayOrigin,0),type = "2")
        }
      }
      hitedLineNumber<-which.min(distances)
      tipOfRayTemp<-intersections[hitedLineNumber,]
      tipOfRays[i,]<-tipOfRayTemp[1:2] 
    }
    
  }
  
  return(tipOfRays)
}



cutAndStretchSpokes_4MultiObject_2D <- function(allSpokesTips,
                                                allSpokesTails,
                                                allMeshes2D) {
  
  tipOfRaysWithLabels<-c()
  for (i in 1:dim(allSpokesTips)[1]) {
    rayOrigin<-allSpokesTails[i,]
    rayDirection<-allSpokesTips[i,]-allSpokesTails[i,]
    intersection<-c()
    p2<-c(colMeans(allSpokesTails),100) #p2 is fixed at the tip of pyramid
    intersectionWithLable<-rep(NA,2+1)
    distanceInit<-Inf
    for (t in 1:length(allMeshes2D)) {
      boundaryPoints1 = allMeshes2D[[t]]
      boundaryPoints2 = rbind(allMeshes2D[[t]][2:dim(allMeshes2D[[t]])[1],],allMeshes2D[[t]][1,])
      for (j in 1:dim(boundaryPoints1)[1]) {
        p1<-c(boundaryPoints1[j,],0) #convert to 3D
        p3<-c(boundaryPoints2[j,],0) #convert to 3D
        tempIntersection<-rayTriangleIntersection(rayOrigin = c(rayOrigin,0), #convert to 3D
                                                  rayDirection = c(rayDirection,0), #convert to 3D
                                                  triangleVertex1 = p1,
                                                  triangleVertex2 = p2,
                                                  triangleVertex3 = p3) 
        
        if(!is.na(tempIntersection[1])){
          tempDis<-norm(tempIntersection-c(rayOrigin,0),type = "2")
          if(tempDis<distanceInit){
            distanceInit<-tempDis
            intersectionWithLable<-c(tempIntersection[1:dim(boundaryPoints1)[2]],t)
          }
        }
      } 
    }
    if(!is.na(intersectionWithLable[1])){
      tipOfRaysWithLabels<-rbind(tipOfRaysWithLabels,intersectionWithLable)
    }
  }
  
  return(tipOfRaysWithLabels)
}




cutAndStretchSpokes_4Bivalves_2D <- function(allSpokesTips,
                                             allSpokesTails,
                                             allMeshes2D) {
  
  tipOfRaysWithLabels<-c()
  for (i in 1:dim(allSpokesTips)[1]) {
    rayOrigin<-allSpokesTails[i,]
    rayDirection<-allSpokesTips[i,]-allSpokesTails[i,]
    intersection<-c()
    p2<-c(colMeans(allSpokesTails),100) #p2 is fixed at the tip of pyramid
    intersectionWithLable<-rep(NA,2+1)
    distanceInit<-Inf
    for (t in 1:length(allMeshes2D)) {
      boundaryPoints1 = allMeshes2D[[t]][1:(dim(allMeshes2D[[t]])[1]-1),]
      boundaryPoints2 = allMeshes2D[[t]][2:dim(allMeshes2D[[t]])[1],]
      for (j in 1:dim(boundaryPoints1)[1]) {
        p1<-c(boundaryPoints1[j,],0) #convert to 3D
        p3<-c(boundaryPoints2[j,],0) #convert to 3D
        tempIntersection<-rayTriangleIntersection(rayOrigin = c(rayOrigin,0), #convert to 3D
                                                  rayDirection = c(rayDirection,0), #convert to 3D
                                                  triangleVertex1 = p1,
                                                  triangleVertex2 = p2,
                                                  triangleVertex3 = p3) 
        
        if(!is.na(tempIntersection[1])){
          tempDis<-norm(tempIntersection-c(rayOrigin,0),type = "2")
          if(tempDis<distanceInit){
            distanceInit<-tempDis
            intersectionWithLable<-c(tempIntersection[1:dim(boundaryPoints1)[2]],t)
          }
        }
      } 
    }
    if(!is.na(intersectionWithLable[1])){
      tipOfRaysWithLabels<-rbind(tipOfRaysWithLabels,intersectionWithLable)
    }
  }
  
  return(tipOfRaysWithLabels)
}



cutAndStretchSpokes_4Bivalves_2D_NA_report <- function(allSpokesTips,
                                             allSpokesTails,
                                             allMeshes2D) {
  
  tipOfRaysWithLabels<-c()
  for (i in 1:dim(allSpokesTips)[1]) {
    rayOrigin<-allSpokesTails[i,]
    rayDirection<-allSpokesTips[i,]-allSpokesTails[i,]
    intersection<-c()
    p2<-c(colMeans(allSpokesTails),100) #p2 is fixed at the tip of pyramid
    intersectionWithLable<-rep(NA,2+1)
    distanceInit<-Inf
    for (t in 1:length(allMeshes2D)) {
      boundaryPoints1 = allMeshes2D[[t]][1:(dim(allMeshes2D[[t]])[1]-1),]
      boundaryPoints2 = allMeshes2D[[t]][2:dim(allMeshes2D[[t]])[1],]
      for (j in 1:dim(boundaryPoints1)[1]) {
        p1<-c(boundaryPoints1[j,],0) #convert to 3D
        p3<-c(boundaryPoints2[j,],0) #convert to 3D
        tempIntersection<-rayTriangleIntersection(rayOrigin = c(rayOrigin,0), #convert to 3D
                                                  rayDirection = c(rayDirection,0), #convert to 3D
                                                  triangleVertex1 = p1,
                                                  triangleVertex2 = p2,
                                                  triangleVertex3 = p3) 
        
        if(!is.na(tempIntersection[1])){
          tempDis<-norm(tempIntersection-c(rayOrigin,0),type = "2")
          if(tempDis<distanceInit){
            distanceInit<-tempDis
            intersectionWithLable<-c(tempIntersection[1:dim(boundaryPoints1)[2]],t)
          }
        }
      } 
    }
    if(!is.na(intersectionWithLable[1])){
      tipOfRaysWithLabels<-rbind(tipOfRaysWithLabels,intersectionWithLable)
    }else{
      tipOfRaysWithLabels<-rbind(tipOfRaysWithLabels,c(NA,NA,NA))
      next
    }
  }
  
  return(tipOfRaysWithLabels)
}




cutAndStretchSpokes_4Bivalves_2D_fast <- function(tempCircle,
                                                  tmeshUp2,
                                                  tmeshDown2,
                                                  thresholdDistance) {
  center<-colMeans(tempCircle)
  
  circleIn3D<-cbind(tempCircle,rep(0,nrow(tempCircle)))
  centerIn3D<-c(center,0)
  circleMesh3D<-as.mesh3d(circleIn3D)
  
  normalsTemp<-circleIn3D-matrix(rep(centerIn3D,nrow(circleIn3D)),ncol = 3,byrow = TRUE)
  
  circleMesh3D$normals<-rbind(apply(normalsTemp,FUN = convertVec2unitVec,MARGIN = 1),
                              rep(1,nrow(circleIn3D)))
  
  # shade3d(circleMesh3D, col="black")  #surface mech
  
  intersectionsUp<-vcgRaySearch(circleMesh3D,mesh = tmeshUp2)
  intersectionsDown<-vcgRaySearch(circleMesh3D,mesh = tmeshDown2)
  
  if(sum(intersectionsUp$quality)<2 | sum(intersectionsDown$quality)<2
     | min(c(intersectionsUp$distance,intersectionsDown$distance))<thresholdDistance){
    return(next)
  }else{
    bothUpDownIntersect<-which(intersectionsUp$quality==1 & intersectionsDown$quality==1)
    
    if(length(bothUpDownIntersect)>0){
      for (i in 1:length(bothUpDownIntersect)) {
        if(intersectionsUp$distance[bothUpDownIntersect[i]]<
           intersectionsDown$distance[bothUpDownIntersect[i]]){
          
          intersectionsUp$quality[bothUpDownIntersect[i]]<-1
          intersectionsDown$quality[bothUpDownIntersect[i]]<-0
        }else{
          intersectionsUp$quality[bothUpDownIntersect[i]]<-0
          intersectionsDown$quality[bothUpDownIntersect[i]]<-1
        }
      }
    }
    
    
    tipOfCuttedUpSpokes<-(t(intersectionsUp$vb)[intersectionsUp$quality==1,1:3])
    tipOfCuttedDownSpokes<-(t(intersectionsDown$vb)[intersectionsDown$quality==1,1:3])
    
    # plot3d(tipOfCuttedUpSpokes,type="s",radius = 0.2,col = "black",expand = 10,box=FALSE,add = TRUE)
    # plot3d(tipOfCuttedDownSpokes,type="s",radius = 0.2,col = "black",expand = 10,box=FALSE,add = TRUE)
    
    # for (i in 1:nrow(tipOfCuttedUpSpokes)) {
    #   plot(rbind(center,tipOfCuttedUpSpokes[i,1:2]),type = 'l',xlim = plotlim,ylim = plotlim,xlab = "",ylab = "")
    #   par(new=TRUE)
    # }
    # for (i in 1:nrow(tipOfCuttedDownSpokes)) {
    #   plot(rbind(center,tipOfCuttedDownSpokes[i,1:2]),type = 'l',xlim = plotlim,ylim = plotlim,xlab = "",ylab = "")
    #   par(new=TRUE)
    # }
    # plot(vert2points(tmeshUp2)[,1:2],type = 'l',xlim = plotlim,ylim = plotlim,xlab = "",ylab = "")
    # par(new=TRUE)
    # plot(vert2points(tmeshDown2)[,1:2],type = 'l',xlim = plotlim,ylim = plotlim,xlab = "",ylab = "")
    
    tipOfRaysWithLabels<-rbind(cbind(tipOfCuttedUpSpokes,rep(1,nrow(tipOfCuttedUpSpokes))),
                               cbind(tipOfCuttedDownSpokes,rep(2,nrow(tipOfCuttedDownSpokes))))
    
    
    return(tipOfRaysWithLabels)  
  }
  
}



