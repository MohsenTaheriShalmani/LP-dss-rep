# generate 3D up and down spokes
generate3DSpokes <- function(meshPDM,
                             medialPoints3D,
                             polyDegree3D=4,
                             k_close=40) {
  # find the average of k closet boundary points
  medialPoints3D_Dis2Boundary<-array(NA,dim = c(numberOfFrames,dim(meshPDM)[1]))
  for (i in 1:numberOfFrames) {
    for (j in 1:dim(meshPDM)[1]) {
      medialPoints3D_Dis2Boundary[i,j]<-norm(medialPoints3D[i,]-meshPDM[j,],type = "2")
    }
  }
  # choose the number of closest points
  averageOf_K_ClosestDistance3D<-rep(NA,numberOfFrames)
  for (i in 1:numberOfFrames) {
    tempArray<-sort(medialPoints3D_Dis2Boundary[i,])
    averageOf_K_ClosestDistance3D[i]<-mean(tempArray[1:k_close])
  }
  radiiAverage3D<-averageOf_K_ClosestDistance3D 
    
  
  smoothMesh<-c()
  numberOfPoints3Dsmoother<-3
  mesh4smooth<-meshPDM
  for (i in 1:dim(PolygonsCsv)[1]) {
    p1<-polyMatrix[i,1]
    p2<-polyMatrix[i,2]
    p3<-polyMatrix[i,3]
    point1<-mesh4smooth[p1,]
    point2<-mesh4smooth[p2,]
    point3<-mesh4smooth[p3,]
    
    smoothMesh<-rbind(smoothMesh,
                      generatePointsBetween3Points(point1 = point1,
                                                   point2 = point2,
                                                   point3 = point3,
                                                   numberOf2DspokePoints = numberOfPoints3Dsmoother))
  }
  
  # x y z
  x<-medialPoints3D[,1]
  y<-medialPoints3D[,2]
  z<-medialPoints3D[,3]
  
  fit4 <- lm(z ~ poly(x, y, degree = polyDegree3D ,raw = TRUE), data=as.data.frame(cbind(z,x,y)))
  
  xData5<-data.frame(x=smoothMesh[,1],y=smoothMesh[,2])
  surfacePointsZ5<-predict(fit4, newdata = xData5)
  
  new3D_mesh<-cbind(smoothMesh[,1:2],surfacePointsZ5)
  
  topPartMeshPDM<-smoothMesh[smoothMesh[,3]>surfacePointsZ5,]
  bottomPartMeshPDM<-smoothMesh[smoothMesh[,3]<surfacePointsZ5,]
  
  #choose the number Of closet points
  N<-round((dim(topPartMeshPDM)[1]/numberOfFrames)*2+1)
  tipOfUpSpokes3D<-array(NA,dim = dim(medialPoints3D))
  pb <- txtProgressBar(style = 3)
  for (i in 1:numberOfFrames) {
    setTxtProgressBar(pb,i/numberOfFrames)
    tempDis<-rep(Inf,dim(topPartMeshPDM)[1])
    for (j in 1:dim(topPartMeshPDM)[1]) {
      tempDis[j]<-norm(medialPoints3D[i,]-topPartMeshPDM[j,],type = "2")
    }
    n_ClosestPoints<-topPartMeshPDM[which.minn(tempDis,n = N),]
    tempDirections<-array(NA,dim = c(N,dim(meshPDM)[2]))
    for (k in 1:N) {
      tempDirections[k,]<-convertVec2unitVec(n_ClosestPoints[k,]-medialPoints3D[i,])
    }
    
    tipOfUpSpokes3D[i,]<-medialPoints3D[i,]+radiiAverage3D[i]*frechetMean(t(tempDirections))
  }
  close(pb)
  #choose the number Of closet points
  N<-round((dim(bottomPartMeshPDM)[1]/numberOfFrames)*2)
  tipOfDownSpokes3D<-array(NA,dim = dim(medialPoints3D))
  pb <- txtProgressBar(style = 3)
  for (i in 1:numberOfFrames) {
    setTxtProgressBar(pb,i/numberOfFrames)
    tempDis<-rep(Inf,dim(bottomPartMeshPDM)[1])
    for (j in 1:dim(bottomPartMeshPDM)[1]) {
      tempDis[j]<-norm(medialPoints3D[i,]-bottomPartMeshPDM[j,],type = "2")
    }
    n_ClosestPoints<-bottomPartMeshPDM[which.minn(tempDis,n = N),]
    tempDirections<-array(NA,dim = c(N,dim(meshPDM)[2]))
    for (k in 1:N) {
      tempDirections[k,]<-convertVec2unitVec(n_ClosestPoints[k,]-medialPoints3D[i,])
    }
    
    tipOfDownSpokes3D[i,]<- medialPoints3D[i,]+radiiAverage3D[i]*frechetMean(t(tempDirections))
  }
  close(pb)
  tipOfCuttedUpSpokes<-cutAndStretchSpokes(BoundaryPDM = tipOfUpSpokes3D,
                                           SkeletalPDM = medialPoints3D,
                                           meshPolygonData = PolygonsCsv,
                                           meshPDM = meshPDM)
  tipOfCuttedDownSpokes<-cutAndStretchSpokes(BoundaryPDM = tipOfDownSpokes3D,
                                             SkeletalPDM = medialPoints3D,
                                             meshPolygonData = PolygonsCsv,
                                             meshPDM = meshPDM)
  
  result<-list("tipOfCuttedUpSpokes"=tipOfCuttedUpSpokes,
               "tipOfCuttedDownSpokes"=tipOfCuttedDownSpokes)
  
  return(result)
}



