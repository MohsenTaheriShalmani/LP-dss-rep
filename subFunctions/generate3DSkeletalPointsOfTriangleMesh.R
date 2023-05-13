generate3DSkeletalPoints <- function(meshPDM,
                                     polyMatrix=polyMatrix,
                                     outward=outward,
                                     numberOf2DspokePoints=numberOf2DspokePoints,
                                     numberOfSpanialPoints=numberOfSpanialPoints,
                                     alpha1=alpha1,
                                     alpha2=alpha2,
                                     threshold4EachLineLength=threshold4EachLineLength,
                                     thresholdAngle=thresholdAngle,
                                     threshold4EachLineLength2D=threshold4EachLineLength2D,
                                     thresholdAngle2D=thresholdAngle2D,
                                     polyDegree3D=polyDegree3D,
                                     polyDegree2D=polyDegree2D){
  
  # pick ray direction and origin
  dimPoly<-dim(polyMatrix)[1]
  allOrigins<-array(NA,dim = dim(polyMatrix))
  allTargets<-array(NA,dim = dim(polyMatrix))
  allOriginNormals<-array(NA,dim = dim(polyMatrix))
  allTargetNormals<-array(NA,dim = dim(polyMatrix))
  pb <- txtProgressBar(style = 3)
  for (i in 1:dimPoly) {
    setTxtProgressBar(pb,i/dimPoly)
    p1<-polyMatrix[i,1]
    p2<-polyMatrix[i,2]
    p3<-polyMatrix[i,3]
    point1<-meshPDM[p1,]
    point2<-meshPDM[p2,]
    point3<-meshPDM[p3,]
    rayOrigin1<-(point1+point2+point3)/3
    normalOfTheTriangle<-unitNormalOfTriangle(point1,point2,point3)
    
    if(outward==TRUE){
      rayDirection1<-(-normalOfTheTriangle)
    }else if(outward==FALSE){
      rayDirection1<-normalOfTheTriangle
    }else{
      stop("Please clearify whether the normals are outward or not!")
    }
    
    allOrigins[i,]<-rayOrigin1
    allOriginNormals[i,]<-normalOfTheTriangle
    
    intersections<-array(NA,dim = c(dim(PolygonsCsv)[1],3))
    for (j in 1:dim(PolygonsCsv)[1]) {
      if(j==i){next}
      p1<-polyMatrix[j,1]
      p2<-polyMatrix[j,2]
      p3<-polyMatrix[j,3]
      point1<-meshPDM[p1,]
      point2<-meshPDM[p2,]
      point3<-meshPDM[p3,]
      tempIntersection<-rayTriangleIntersection(rayOrigin = rayOrigin1 ,
                                                rayDirection = rayDirection1,
                                                triangleVertex1 = point1,
                                                triangleVertex2 = point2,
                                                triangleVertex3 = point3)
      intersections[j,]<-tempIntersection
    }
    if(sum(is.na(intersections[,1]))==dim(PolygonsCsv)[1]){ #for cases that normal is outward
      
      rayDirection1<-(-rayDirection1)
      
      intersections<-array(NA,dim = c(dim(PolygonsCsv)[1],3))
      for (j in 1:dim(PolygonsCsv)[1]) {
        if(j==i){next}
        p1<-polyMatrix[j,1]
        p2<-polyMatrix[j,2]
        p3<-polyMatrix[j,3]
        point1<-meshPDM[p1,]
        point2<-meshPDM[p2,]
        point3<-meshPDM[p3,]
        tempIntersection<-rayTriangleIntersection(rayOrigin = rayOrigin1 ,
                                                  rayDirection = rayDirection1,
                                                  triangleVertex1 = point1,
                                                  triangleVertex2 = point2,
                                                  triangleVertex3 = point3)
        intersections[j,]<-tempIntersection
      }
    }
    distances<-rep(Inf,dim(PolygonsCsv)[1])
    for (k in 1:dim(PolygonsCsv)[1]) {
      if(!is.na(intersections[k,1])){
        distances[k]<-norm(intersections[k,]-rayOrigin1,type = "2")
      }
    }
    hitedTriangleNumber<-which.min(distances)
    tipOfShootedArrow<-intersections[hitedTriangleNumber,]
    
    p1<-polyMatrix[hitedTriangleNumber,1]
    p2<-polyMatrix[hitedTriangleNumber,2]
    p3<-polyMatrix[hitedTriangleNumber,3]
    point1<-meshPDM[p1,]
    point2<-meshPDM[p2,]
    point3<-meshPDM[p3,]
    
    normalOfTheHitedTriangle<-unitNormalOfTriangle(point1,point2,point3)
    
    allTargets[i,]<-tipOfShootedArrow
    allTargetNormals[i,]<-normalOfTheHitedTriangle
    
  }
  close(pb)
  
  allGeodesicDistances<-array(NA,dimPoly)
  for (i in 1:dimPoly) {
    allGeodesicDistances[i]<-acos(pmin(pmax(as.numeric(allOriginNormals[i,]%*%allTargetNormals[i,]) ,-1.0), 1.0))
  }
  
  allInternalLineLengths<-array(NA,dimPoly)
  for (i in 1:dimPoly) {
    allInternalLineLengths[i]<-norm(allTargets[i,]-allOrigins[i,],type = "2")
  }
  
  allMiddlepoints<-array(NA,dim = dim(polyMatrix))
  for (i in 1:dimPoly) {
    allMiddlepoints[i,]<-(allOrigins[i,]+allTargets[i,])/2
  }
  
  allClosestDistances2Boundary<-rep(NA,dimPoly)
  pb <- txtProgressBar(style = 3)
  for (i in 1:dimPoly) {
    setTxtProgressBar(pb,i/dimPoly)
    distances<-rep(Inf,numberOfPoints)
    for (j in 1:numberOfPoints) {
      distances[j]<-norm(allMiddlepoints[i,]-meshPDM[j,],type = "2")
    }
    allClosestDistances2Boundary[i]<-distances[which.min(distances)]
  }
  close(pb)
  
  selectedTriangles<-which(allGeodesicDistances>thresholdAngle &
                             allInternalLineLengths<2*threshold4EachLineLength*allClosestDistances2Boundary)
  
  
  selectedMiddlePoint<-allMiddlepoints[selectedTriangles,]
  
  
  # smoother
  # x y z
  x<-selectedMiddlePoint[,1]
  y<-selectedMiddlePoint[,2]
  z<-selectedMiddlePoint[,3]
  
  fit4 <- lm(z ~ poly(x, y, degree = polyDegree3D ,raw = TRUE), data=as.data.frame(cbind(z,x,y)))
  # summary(fit2)
  
  newDATA<-data.frame(x=x, y=y)
  surfacePointsZ<-predict(fit4, newdata = newDATA)
  
  medialPoints<-cbind(x,y,surfacePointsZ)
  
  projectedOnXYplane<-cbind(x,y)
  
  # Alpha-convex hull
  ahull.obj <- ahull(projectedOnXYplane, alpha = alpha1)
  
  indices<-ahull.obj$ashape.obj$alpha.extremes
  ahull.obj$ashape.obj$edges
  
  #boundaryPoints1 is connected to boundaryPoints2 to form the edges
  boundaryPoints1<-ahull.obj$ashape.obj$edges[,3:4]
  boundaryPoints2<-ahull.obj$ashape.obj$edges[,5:6]
  
  allNormalVecTips<-array(NA,dim = c(dim(boundaryPoints1)[1],2))
  allNormalVecTails<-array(NA,dim = c(dim(boundaryPoints1)[1],2))
  allNormalsDirections2D<-array(NA,dim = c(dim(boundaryPoints1)[1],2))
  for (i in 1:dim(boundaryPoints1)[1]) {
    dx_dy<-boundaryPoints2[i,]-boundaryPoints1[i,] #normal is c(-dy,dx) where dy=y2-y1 and dx=x2-x1
    normalTemp1<-convertVec2unitVec(c(-dx_dy[2],dx_dy[1]))
    normalTemp1
    
    normalTemp2<-(-normalTemp1)
    
    lineSegmentMean<-(boundaryPoints1[i,]+boundaryPoints2[i,])/2
    normalVecTip1<-lineSegmentMean+normalTemp1
    normalVecTip2<-lineSegmentMean+normalTemp2
    
    distances1<-rep(NA,dim(projectedOnXYplane)[1])
    distances2<-rep(NA,dim(projectedOnXYplane)[1])
    for (j in 1:dim(projectedOnXYplane)[1]) {
      distances1[j]<-norm(normalVecTip1-projectedOnXYplane[j,],type = "2")
      distances2[j]<-norm(normalVecTip2-projectedOnXYplane[j,],type = "2")
    }
    if(min(distances1)<min(distances2)){
      normalVecTip<-normalVecTip1
      normalTemp<-normalTemp1
    }else{
      normalVecTip<-normalVecTip2
      normalTemp<-normalTemp2
    }
    
    allNormalsDirections2D[i,]<-normalTemp
    allNormalVecTips[i,]<-normalVecTip
    allNormalVecTails[i,]<-lineSegmentMean
  }
  
  boundaryPoints2D<-rbind(boundaryPoints1,boundaryPoints2)
  
  
  # rayTriangleIntersection is in 3D so for 2D we consider the second point of the
  # triangle as c(colMeans(boundaryPoints1),100) which is the vertex of the pyramid
  tipOfRays<-array(NA,dim = dim(boundaryPoints1))
  allTargetNormals2DMesh<-array(NA,dim = dim(boundaryPoints1))
  pb <- txtProgressBar(style = 3)
  for (i in 1:dim(boundaryPoints1)[1]) {
    setTxtProgressBar(pb,i/dim(boundaryPoints1)[1])
    rayOrigin<-allNormalVecTails[i,]
    rayDirection<-allNormalVecTips[i,]-allNormalVecTails[i,]
    intersections<-array(NA,dim = c(dim(boundaryPoints1)[1],3))
    p2<-c(colMeans(boundaryPoints1),100) #p2 is fixed at the tip of pyramid
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
    if(sum(is.na(intersections[,1]))==dim(boundaryPoints1)[1]){ #for cases that normal is outward
      
      rayDirection<-(-rayDirection)
      
      intersections<-array(NA,dim = c(dim(boundaryPoints1)[1],3))
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
    }
    distances<-rep(Inf,dim(boundaryPoints1)[1]) #find the closest intersect
    for (k in 1:dim(boundaryPoints1)[1]) {
      if(!is.na(intersections[k,1])){
        distances[k]<-norm(intersections[k,]-c(rayOrigin,0),type = "2")
      }
    }
    distances[i]<-Inf #we do not consider the ray origin 
    hitedLineNumber<-which.min(distances)
    if(is.na(intersections[hitedLineNumber,1])){
      tipOfRays[i,]<-rayOrigin
    }else{
      tipOfRayTemp<-intersections[hitedLineNumber,]
      tipOfRays[i,]<-tipOfRayTemp[1:2] 
    }
    
    # normal of the hitted line segment
    dx_dy<-boundaryPoints2[hitedLineNumber,]-boundaryPoints1[hitedLineNumber,]
    
    normalTemp1<-convertVec2unitVec(c(-dx_dy[2],dx_dy[1]))
    normalTemp1
    
    normalTemp2<-(-normalTemp1)
    
    lineSegmentMean<-(boundaryPoints1[hitedLineNumber,]+boundaryPoints2[hitedLineNumber,])/2
    normalVecTip1<-lineSegmentMean+normalTemp1
    normalVecTip2<-lineSegmentMean+normalTemp2
    
    distances1<-rep(NA,dim(projectedOnXYplane)[1])
    distances2<-rep(NA,dim(projectedOnXYplane)[1])
    for (j in 1:dim(projectedOnXYplane)[1]) {
      distances1[j]<-norm(normalVecTip1-projectedOnXYplane[j,],type = "2")
      distances2[j]<-norm(normalVecTip2-projectedOnXYplane[j,],type = "2")
    }
    if(min(distances1)<min(distances2)){
      normalVecTip<-normalVecTip1
      normalTemp<-normalTemp1
    }else{
      normalVecTip<-normalVecTip2
      normalTemp<-normalTemp2
    }
    
    allTargetNormals2DMesh[i,]<-normalTemp
    
  }
  close(pb)
  
  # find the spine
  allGeodesicDistances2D<-array(NA,dim(boundaryPoints1)[1])
  for (i in 1:dim(boundaryPoints1)[1]) {
    allGeodesicDistances2D[i]<-acos(pmin(pmax(as.numeric(allNormalsDirections2D[i,]%*%allTargetNormals2DMesh[i,]) ,-1.0), 1.0))
  }
  
  allInternalLineLengths2D<-array(NA,dim(boundaryPoints1)[1])
  for (i in 1:dim(boundaryPoints1)[1]) {
    allInternalLineLengths2D[i]<-norm(tipOfRays[i,]-allNormalVecTails[i,],type = "2")
  }
  
  allMiddlepoints2D<-array(NA,dim = dim(boundaryPoints1))
  for (i in 1:dim(boundaryPoints1)[1]) {
    allMiddlepoints2D[i,]<-(tipOfRays[i,]+allNormalVecTails[i,])/2
  }
  
  allBoundaryPoints2D<-rbind(boundaryPoints1,boundaryPoints2)
  allClosestDistances2Boundary2D<-rep(NA,dim(allMiddlepoints2D)[1])
  for (i in 1:dim(allMiddlepoints2D)[1]) {
    distances<-rep(Inf,dim(allBoundaryPoints2D)[1])
    for (j in 1:dim(allBoundaryPoints2D)[1]) {
      distances[j]<-norm(allMiddlepoints2D[i,]-allBoundaryPoints2D[j,],type = "2")
    }
    allClosestDistances2Boundary2D[i]<-distances[which.min(distances)]
  }
  
  selectedLineSegments<-which(allGeodesicDistances2D>thresholdAngle2D & 
                                allInternalLineLengths2D<2*threshold4EachLineLength2D*
                                allClosestDistances2Boundary2D)
  
  selectedMiddlePoint2D<-allMiddlepoints2D[selectedLineSegments,]
  
  # x y z
  x<-selectedMiddlePoint2D[,1]
  y<-selectedMiddlePoint2D[,2]
  
  fit4_2D <- lm(y ~ poly(x, degree = polyDegree2D,raw = TRUE),
                data=as.data.frame(cbind(x,y))) #NB!!! raw=F calculate orthogonal polynomial
  
  xData1<-data.frame(x=x)
  curvePoints1<-predict(fit4_2D, newdata = xData1)
  
  medialPoints_2D<-cbind(x,curvePoints1)
  
  minVal<-min(medialPoints_2D[,1])
  maxVal<-max(medialPoints_2D[,1])
  
  # tail to head
  xValues<-seq(minVal, maxVal, length = numberOfSpanialPoints)
  # #head to tail
  # xValues<-seq(maxVal, minVal, length = numberOfSpanialPoints) #head to tail
  
  xData2<-data.frame(x=xValues)
  curvePoints2<-predict(fit4_2D, newdata = xData2)
  
  medialPointsEqualyDis_2D<-cbind(xValues,curvePoints2)
  
  # find the average of k closet boundary points
  medialPoints2D_Dis2Boundary<-array(NA,dim = c(dim(medialPointsEqualyDis_2D)[1],dim(allBoundaryPoints2D)[1]))
  for (i in 1:dim(medialPointsEqualyDis_2D)[1]) {
    for (j in 1:dim(allBoundaryPoints2D)[1]) {
      medialPoints2D_Dis2Boundary[i,j]<-norm(medialPointsEqualyDis_2D[i,]-allBoundaryPoints2D[j,],type = "2")
    }
  }
  
  # choose the number of closest points
  k<-10
  averageOf_K_ClosestDistance<-rep(NA,dim(medialPointsEqualyDis_2D)[1])
  for (i in 1:dim(medialPointsEqualyDis_2D)[1]) {
    tempArray<-sort(medialPoints2D_Dis2Boundary[i,])
    averageOf_K_ClosestDistance[i]<-mean(tempArray[1:k])
  }
  
  radiiAverage2D<-averageOf_K_ClosestDistance
  
  regeneratedInternalPoints<-c()
  for (i in 1:dim(medialPointsEqualyDis_2D)[1]) {
    a<-averageOf_K_ClosestDistance[i]
    b<-averageOf_K_ClosestDistance[i]
    regeneratedInternalPoints<-rbind(regeneratedInternalPoints,
                                     ellipsoidGenerator_2D(center = medialPointsEqualyDis_2D[i,]
                                                           ,a = a,b = b, 10, 1))
  }
  
  # find bounday points by alpha convex hull
  # Alpha-convex hull
  ahull.obj <- ahull(regeneratedInternalPoints, alpha = alpha2)
  
  indices<-ahull.obj$ashape.obj$alpha.extremes
  # ahull.obj$ashape.obj$edges
  
  #boundaryPoints1 is connected to boundaryPoints2 to form the edges
  boundaryPoints1<-ahull.obj$ashape.obj$edges[,3:4]
  boundaryPoints2<-ahull.obj$ashape.obj$edges[,5:6]
  
  #smooth the boundary by adding points between each pairs
  new2D_mesh<-c()
  numberOfPoints<- 50
  for (i in 1:dim(boundaryPoints1)[1]) {
    
    tempPoints<-generatePointsBetween2Points(boundaryPoints1[i,],
                                             boundaryPoints2[i,],
                                             numberOfPoints = numberOfPoints)
    
    new2D_mesh<-rbind(new2D_mesh,tempPoints)
    
  }
  
  # generate spokes 2D
  xData4<-data.frame(x=new2D_mesh[,1])
  curvePointsY4<-predict(fit4_2D, newdata = xData4)
  medialPointstest<-cbind(new2D_mesh[,1],curvePointsY4)
  
  topPart<-new2D_mesh[new2D_mesh[,2]>curvePointsY4,]
  bottomPart<-new2D_mesh[new2D_mesh[,2]<curvePointsY4,]
  
  
  #choose the number Of closet points
  N<-round((dim(topPart)[1]/dim(medialPointsEqualyDis_2D)[1])/2)
  # find spokes
  tipOfTopSpokes<-array(NA,dim = dim(medialPointsEqualyDis_2D))
  pb <- txtProgressBar(style = 3)
  for (i in 1:dim(medialPointsEqualyDis_2D)[1]) {
    setTxtProgressBar(pb,i/dim(medialPointsEqualyDis_2D)[1])
    tempDis<-rep(Inf,dim(topPart)[1])
    for (j in 1:dim(topPart)[1]) {
      tempDis[j]<-norm(medialPointsEqualyDis_2D[i,]-topPart[j,],type = "2")
    }
    n_ClosestPoints<-topPart[which.minn(tempDis,n = N),]
    tempDirections<-array(NA,dim = c(N,dim(topPart)[2]))
    for (k in 1:N) {
      tempDirections[k,]<-convertVec2unitVec(n_ClosestPoints[k,]-medialPointsEqualyDis_2D[i,])
    }
    
    tipOfTopSpokes[i,]<-medialPointsEqualyDis_2D[i,]+radiiAverage2D[i]*frechetMean(t(tempDirections))
  }
  close(pb)
  #choose the number Of closet points
  N<-round((dim(bottomPart)[1]/dim(medialPointsEqualyDis_2D)[1])/2)
  # find spokes
  tipOfBottomSpokes<-array(NA,dim = dim(medialPointsEqualyDis_2D))
  pb <- txtProgressBar(style = 3)
  for (i in 1:dim(medialPointsEqualyDis_2D)[1]) {
    setTxtProgressBar(pb,i/dim(medialPointsEqualyDis_2D)[1])
    tempDis<-rep(Inf,dim(bottomPart)[1])
    for (j in 1:dim(bottomPart)[1]) {
      tempDis[j]<-norm(medialPointsEqualyDis_2D[i,]-bottomPart[j,],type = "2")
    }
    n_ClosestPoints<-bottomPart[which.minn(tempDis,n = N),]
    tempDirections<-array(NA,dim = c(N,dim(bottomPart)[2]))
    for (k in 1:N) {
      tempDirections[k,]<-convertVec2unitVec(n_ClosestPoints[k,]-medialPointsEqualyDis_2D[i,])
    }
    
    tipOfBottomSpokes[i,]<-medialPointsEqualyDis_2D[i,]+radiiAverage2D[i]*frechetMean(t(tempDirections))
  }
  close(pb)
  
  extraSpokeHead<-(convertVec2unitVec(tipOfTopSpokes[1,]-medialPointsEqualyDis_2D[1,])+
                     convertVec2unitVec(tipOfBottomSpokes[1,]-medialPointsEqualyDis_2D[1,]))/2+
    medialPointsEqualyDis_2D[1,]
  
  extraSpokeTail<-(convertVec2unitVec(tipOfTopSpokes[numberOfSpanialPoints,]-medialPointsEqualyDis_2D[numberOfSpanialPoints,])+
                     convertVec2unitVec(tipOfBottomSpokes[numberOfSpanialPoints,]-medialPointsEqualyDis_2D[numberOfSpanialPoints,]))/2+
    medialPointsEqualyDis_2D[numberOfSpanialPoints,]
  twoHeadAndTailSpokes<-rbind(extraSpokeHead,extraSpokeTail)
  
  #cut and strech to reach the boundary
  tipOfCuttedTopSpokes2D<-cutAndStretchSpokes2D(allSpokesTips = tipOfTopSpokes,
                                                allSpokesTails = medialPointsEqualyDis_2D,
                                                boundaryPoints1 = boundaryPoints1,
                                                boundaryPoints2 = boundaryPoints2)
  tipOfCuttedBottomSpokes2D<-cutAndStretchSpokes2D(allSpokesTips = tipOfBottomSpokes,
                                                   allSpokesTails = medialPointsEqualyDis_2D,
                                                   boundaryPoints1 = boundaryPoints1,
                                                   boundaryPoints2 = boundaryPoints2)
  
  tipOfCuttedHeadAndTailSpokes2D<-cutAndStretchSpokes2D(allSpokesTips = twoHeadAndTailSpokes,
                                                        allSpokesTails = medialPointsEqualyDis_2D[c(1,numberOfSpanialPoints),],
                                                        boundaryPoints1 = boundaryPoints1,
                                                        boundaryPoints2 = boundaryPoints2)
  for (i in 1:dim(tipOfTopSpokes)[1]) {
    if(is.na(tipOfCuttedTopSpokes2D[i,1])){
      tipOfTopSpokes[i,]<-tipOfTopSpokes[i,]
    }else{
      tipOfTopSpokes[i,]<-tipOfCuttedTopSpokes2D[i,]
    }
  }
  for (i in 1:dim(tipOfBottomSpokes)[1]) {
    if(is.na(tipOfCuttedBottomSpokes2D[i,1])){
      tipOfBottomSpokes[i,]<-tipOfBottomSpokes[i,]
    }else{
      tipOfBottomSpokes[i,]<-tipOfCuttedBottomSpokes2D[i,]
    }
  }
  tipOfExtraSpokeHead<-tipOfCuttedHeadAndTailSpokes2D[1,]
  tipOfExtraSpokeTail<-tipOfCuttedHeadAndTailSpokes2D[2,]
  
  # spokes interpolation
  topSpokesPoints<-c()
  for (i in 1:dim(medialPointsEqualyDis_2D)[1]) {
    tempPoints<-generatePointsBetween2Points(medialPointsEqualyDis_2D[i,],
                                             tipOfTopSpokes[i,],
                                             numberOfPoints = numberOf2DspokePoints) 
    topSpokesPoints<-rbind(topSpokesPoints,tempPoints)
  }
  bottomSpokesPoints<-c()
  for (i in 1:dim(medialPointsEqualyDis_2D)[1]) {
    tempPoints<-generatePointsBetween2Points(medialPointsEqualyDis_2D[i,],
                                             tipOfBottomSpokes[i,],
                                             numberOfPoints = numberOf2DspokePoints) 
    bottomSpokesPoints<-rbind(bottomSpokesPoints,tempPoints)
  }
  extraSpokeHeadPoints<-generatePointsBetween2Points(medialPointsEqualyDis_2D[1,],
                                                     tipOfExtraSpokeHead,
                                                     numberOfPoints = numberOf2DspokePoints)
  extraSpokeTailPoints<-generatePointsBetween2Points(medialPointsEqualyDis_2D[numberOfSpanialPoints,],
                                                     tipOfExtraSpokeTail,
                                                     numberOfPoints = numberOf2DspokePoints)
  
  
  skeletalPoints_2D<-c()
  k1<-1
  k2<-1
  for (i in 1:(dim(topSpokesPoints)[1]/numberOf2DspokePoints)) {
    range1<-c(k1:(k1+numberOf2DspokePoints-1))
    range2<-c((k2+1):(k2+numberOf2DspokePoints-1))
    skeletalPoints_2D<-rbind(skeletalPoints_2D,topSpokesPoints[range1,])
    skeletalPoints_2D<-rbind(skeletalPoints_2D,bottomSpokesPoints[range2,])
    k1<-k1+numberOf2DspokePoints
    k2<-k2+numberOf2DspokePoints
  }
  
  #add two extra head and tail points
  skeletalPoints_2D<-rbind(extraSpokeHeadPoints[-1,],
                           skeletalPoints_2D,
                           extraSpokeTailPoints[-1,])
  
  
  newData_2d<-data.frame(x=skeletalPoints_2D[,1], y=skeletalPoints_2D[,2])
  surfacePointsZ2<-predict(fit4, newdata = newData_2d)
  medialPoints3D<-cbind(skeletalPoints_2D[,1],skeletalPoints_2D[,2],surfacePointsZ2)
  
  return(medialPoints3D)
}

