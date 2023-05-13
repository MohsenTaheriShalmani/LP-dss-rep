
generate3DSkeletalPoints2 <- function(tmesh=tmesh,
                                      plotting=TRUE,
                                      plotSubregions=TRUE,
                                      numberOf2DspokePoints=numberOf2DspokePoints,
                                      numberOfSpanialPoints=numberOfSpanialPoints,
                                      alpha1=alpha1,
                                      alpha2=alpha2,
                                      k_closestPoint2D,
                                      k_closestPoint3D,
                                      threshold4EachLineLength=threshold4EachLineLength,
                                      thresholdAngle=thresholdAngle,
                                      threshold4EachLineLength2D=threshold4EachLineLength2D,
                                      thresholdAngle2D=thresholdAngle2D,
                                      polyDegree3D=polyDegree3D,
                                      polyDegree2D=polyDegree2D,
                                      numberOfPoints2Dremesh=numberOfPoints2Dremesh,
                                      increase3DBoundaryPoint=FALSE){
  
  
  # check the mesh orientation because the normals must be outwards
  if (checkFaceOrientation(tmesh)==FALSE) {
    tmesh <- Morpho::invertFaces(tmesh)
  }
  
  #plot
  if(plotting==TRUE){
    open3d()
    # wire3d(tmesh, col="grey")  #wire mesh
    shade3d(tmesh, col="white",alpha=0.2)  #surface mesh
    plotNormals(tmesh)
  }
  
  #make a new mesh with inward normals
  tmesh2<-tmesh
  tmesh2$normals[1:3,]<-(-tmesh2$normals[1:3,])
  #move the mesh points a little inwards
  tmesh2$vb[1:3,]<-tmesh2$vb[1:3,]+0.5*tmesh2$normals[1:3,]
  
  #plot
  if(plotting==TRUE){
    open3d()
    shade3d(tmesh, col="white",alpha=0.2)  #surface mech
    shade3d(tmesh2, col="white",alpha=0.2)  #surface mech 
  }
  
  #calculate outward normals of the vertices
  tmesh2 <- vcgUpdateNormals(tmesh2)
  #make a new mesh with inward normals
  tmesh2$normals[1:3,]<-(-tmesh2$normals[1:3,])
  #find the intercetions of rays with internal normal directions
  intersections <- vcgRaySearch(tmesh2,mesh = tmesh)
  
  #plot
  if(plotting==TRUE){
    open3d()
    spheres3d(vert2points(intersections),col="blue",radius = 0.2) #plot intersections
    #NB!!! we have the information of normals at intersections intersections$normals !!!!
    shade3d(tmesh, col="white",alpha=0.2)  #surface mech
    plotNormals(intersections) 
  }
  
  allMiddlepoints<-(t(tmesh$vb[1:3,])+t(intersections$vb[1:3,]))/2
  
  #plot
  if(plotting==TRUE){
    open3d()
    shade3d(tmesh, col="white",alpha=0.2)  #surface mech
    spheres3d(allMiddlepoints,col="black",radius = 0.2) #plot intersections
  }
  
  allGeodesicDistances<-acos(pmin(pmax(as.numeric(rowSums(t(tmesh$normals[1:3,])*t(intersections$normals[1:3,]))))))
  
  tempMatrix<-t(tmesh$vb[1:3,])-t(intersections$vb[1:3,])
  allInternalLineLengths<-apply(tempMatrix, 1, myNorm)
  
  allClosestDistances2Boundary<-rep(NA,dim(allMiddlepoints)[1])
  pb <- txtProgressBar(style = 3)
  for (i in 1:dim(allMiddlepoints)[1]) {
    setTxtProgressBar(pb,i/dim(allMiddlepoints)[1])
    
    distances<-rep(Inf,dim(t(tmesh$vb[1:3,]))[1])
    tempMatrix1<-matrix(rep(allMiddlepoints[i,],dim(t(tmesh$vb[1:3,]))[1]),ncol = 3,byrow = T)
    tempMatrix2<-tempMatrix1-t(tmesh$vb[1:3,])
    distances<-apply(tempMatrix2, 1, myNorm)
    
    allClosestDistances2Boundary[i]<-distances[which.min(distances)]
  }
  close(pb)
  
  selectedTriangles<-which(allGeodesicDistances>thresholdAngle &
                             allInternalLineLengths<2*threshold4EachLineLength*allClosestDistances2Boundary)
  length(selectedTriangles)
  
  if(plotting==TRUE){
    open3d()
    plot3d(allMiddlepoints[selectedTriangles,],type="p",col = "black",expand = 10,box=FALSE,add = TRUE)
    # wire3d(tmesh, col="blue")  #wire mesh
    shade3d(tmesh, col="blue",alpha=0.2) 
  }
  
  ####################################################################################################
  ####################################################################################################
  # smoother
  
  sphereCenters<-allMiddlepoints[selectedTriangles,]
  
  selectedMiddlePoint<-sphereCenters
  
  # x y z
  x<-selectedMiddlePoint[,1]
  y<-selectedMiddlePoint[,2]
  z<-selectedMiddlePoint[,3]
  
  fit4 <- lm(z ~ poly(x, y, degree = polyDegree3D ,raw = TRUE), data=as.data.frame(cbind(z,x,y)))
  # summary(fit2)
  
  newDATA<-data.frame(x=x, y=y)
  surfacePointsZ<-predict(fit4, newdata = newDATA)
  
  medialPoints<-cbind(x,y,surfacePointsZ)
  
  #plot
  if(plotting==TRUE){
    open3d()
    plot3d(medialPoints,type="p",col = "black",expand = 10,box=FALSE,add = TRUE)
    shade3d(tmesh, col="blue",alpha=0.2) 
  }
  
  projectedOnXYplane<-cbind(x,y)

  #define plot limit
  plotlim <- c((min(projectedOnXYplane[,1])-10),(max(projectedOnXYplane[,1])+10))
  
  #plot
  if(plotting==TRUE){
    plotshapes(projectedOnXYplane)
    plot(projectedOnXYplane,pch=20,xlim = plotlim,ylim = plotlim,xlab = "",ylab = "") 
  }
  
  # Alpha-convex hull
  ahull.obj <- ahull(projectedOnXYplane, alpha = alpha1)
  
  #plot
  if(plotting==TRUE){
    plot(ahull.obj) 
  }
  
  indices<-ahull.obj$ashape.obj$alpha.extremes
  ahull.obj$ashape.obj$edges
  
  #plot
  if(plotting==TRUE){
    plot(projectedOnXYplane) # all points
    plotshapes(projectedOnXYplane[indices,]) #boundary points 
  }
  
  #boundaryPoints1 is connected to boundaryPoints2 to form the edges
  boundaryPoints1<-ahull.obj$ashape.obj$edges[,3:4]
  boundaryPoints2<-ahull.obj$ashape.obj$edges[,5:6]
  
  #plot 2D
  if(plotting==TRUE){
    for(i in 1:dim(boundaryPoints1)[1]){
      plot(rbind(boundaryPoints1[i,],boundaryPoints2[i,]),type="l",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "")
      par(new=TRUE)
    } 
  }
  
  # plot 2D mesh and its normal vectors
  allNormalVecTips<-array(NA,dim = c(dim(boundaryPoints1)[1],2))
  allNormalVecTails<-array(NA,dim = c(dim(boundaryPoints1)[1],2))
  allNormalsDirections2D<-array(NA,dim = c(dim(boundaryPoints1)[1],2))
  for (i in 1:dim(boundaryPoints1)[1]) {
    dx_dy<-boundaryPoints2[i,]-boundaryPoints1[i,] #normal is c(-dy,dx) where dy=y2-y1 and dx=x2-x1
    
    normalTemp<-convertVec2unitVec(c(-dx_dy[2],dx_dy[1]))
    
    lineSegmentMean<-(boundaryPoints1[i,]+boundaryPoints2[i,])/2
    normalVecTip<-lineSegmentMean+normalTemp
    
    #check the tip of the normal is inside the alpha convex hall
    normalTipIsInHull<-inahull(ahull.obj, normalVecTip)
    
    if(normalTipIsInHull==FALSE){
      normalTemp<-(-normalTemp)
      normalVecTip<-lineSegmentMean+normalTemp
    }
    
    allNormalsDirections2D[i,]<-normalTemp
    allNormalVecTips[i,]<-normalVecTip
    allNormalVecTails[i,]<-lineSegmentMean
  }
  
  #plot 2D
  if(plotting==TRUE){
    for(i in 1:dim(boundaryPoints1)[1]){
      plot(rbind(boundaryPoints1[i,],boundaryPoints2[i,]),type="l",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "")
      par(new=TRUE)
    } 
    # plot 2D mesh normals
    for (i in 1:dim(boundaryPoints1)[1]) {
      plot(rbind(allNormalVecTails[i,],allNormalVecTips[i,]),type="l",
           xlim = plotlim,ylim = plotlim,xlab = "",ylab = "")
      par(new=TRUE) 
    }
  }
  
  boundaryPoints2D<-rbind(boundaryPoints1,boundaryPoints2)
  
  # rayTriangleIntersection is in 3D so for 2D we consider the second point of the
  # triangle as c(colMeans(boundaryPoints1),100) which is the vertex of the pyramid
  tipOfRays<-array(NA,dim = dim(boundaryPoints1))
  allTargetNormals2DMesh<-array(NA,dim = dim(boundaryPoints1))
  allTargetNormalsTip<-array(NA,dim = dim(boundaryPoints1))
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
    
    normalTemp<-convertVec2unitVec(c(-dx_dy[2],dx_dy[1]))
    
    normalVecTip<-tipOfRays[i,]+normalTemp
    
    #check the tip of the normal is inside the alpha convex hall
    normalTipIsInHull<-inahull(ahull.obj, normalVecTip)
    
    if(normalTipIsInHull==FALSE){
      normalTemp<-(-normalTemp)
      normalVecTip<-tipOfRays[i,]+normalTemp
    }
    
    allTargetNormals2DMesh[i,]<-normalTemp
    allTargetNormalsTip[i,]<-normalVecTip
  }
  close(pb)
  
  #plot 2D
  if(plotting==TRUE){
    for(i in 1:dim(boundaryPoints1)[1]){
      plot(rbind(boundaryPoints1[i,],boundaryPoints2[i,]),type="l",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "")
      par(new=TRUE)
      plot(rbind(allNormalVecTails[i,],allNormalVecTips[i,]),type="l",
           xlim = plotlim,ylim = plotlim,xlab = "",ylab = "")
      par(new=TRUE) 
      plot(rbind(allNormalVecTails[i,],tipOfRays[i,]),type="l",
           xlim = plotlim,ylim = plotlim,xlab = "",ylab = "")
      par(new=TRUE)
    } 
  }
  
  #plot 2D
  if(plotting==TRUE){
    for (i in 1:dim(boundaryPoints1)[1]) {
      plot(rbind(boundaryPoints1[i,],boundaryPoints2[i,]),type="l",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "")
      par(new=TRUE)
      plot(rbind(tipOfRays[i,],allTargetNormalsTip[i,]),type="l",
           xlim = plotlim,ylim = plotlim,xlab = "",ylab = "")
      par(new=TRUE)
      
    } 
  }
  
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
  
  # length(selectedLineSegments)
  
  #plot 2D
  if(plotting==TRUE){
    for(i in 1:dim(boundaryPoints1)[1]){
      plot(rbind(boundaryPoints1[i,],boundaryPoints2[i,]),type="l",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "")
      par(new=TRUE)
    }
    par(new=TRUE)
    plot(allMiddlepoints2D,pch=20,col="red",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(allMiddlepoints2D[selectedLineSegments,],pch=20,col="blue",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "")
  }
  
  centroid2D<-colMeans(allMiddlepoints2D[selectedLineSegments,])
  centeredPoints2D<-scale(allMiddlepoints2D[selectedLineSegments,],center = TRUE,scale = FALSE)
  pca2D<-prcomp(centeredPoints2D)
  
  
  # We need to control PCA rotation to preserve correspondence
  v1<-c(0,1)%*%(pca2D$rotation)
  v2<-c(1,0)%*%(pca2D$rotation)
  d_g1<-acos(pmin(pmax(sum(c(0,1)*v1) ,-1.0), 1.0))
  d_g2<-acos(pmin(pmax(sum(c(1,0)*v2) ,-1.0), 1.0))
  if(d_g1>(pi/2) & d_g2>(pi/2)){
    rotationMatrix2D<-(pca2D$rotation)*rbind(c(-1,-1),c(-1,-1))
    selectedMiddlePoint2D<-centeredPoints2D%*%rotationMatrix2D
  }else if(d_g1>(pi/2)){
    rotationMatrix2D<-(pca2D$rotation)*rbind(c(1,-1),c(1,-1))
    selectedMiddlePoint2D<-centeredPoints2D%*%rotationMatrix2D
  }else if(d_g2>(pi/2)){
    rotationMatrix2D<-(pca2D$rotation)*rbind(c(-1,1),c(-1,1))
    selectedMiddlePoint2D<-centeredPoints2D%*%rotationMatrix2D
  }else{
    rotationMatrix2D<-pca2D$rotation
    selectedMiddlePoint2D<-pca2D$x
  }
  
  #plot 2D
  if(plotting==TRUE){
    plot(selectedMiddlePoint2D,pch=20,col="blue",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "") 
  }
  
  # selectedMiddlePoint2D<-allMiddlepoints2D[selectedLineSegments,]
  
  # x y z
  x<-selectedMiddlePoint2D[,1]
  y<-selectedMiddlePoint2D[,2]
  
  fit4_2D <- lm(y ~ poly(x, degree = polyDegree2D,raw = TRUE),
                data=as.data.frame(cbind(x,y))) #NB!!! raw=F calculate orthogonal polynomial
  
  xData1<-data.frame(x=x)
  curvePoints1<-predict(fit4_2D, newdata = xData1)
  
  medialPoints_2D<-cbind(x,curvePoints1)
  
  #plot 2D
  if(plotting==TRUE){
    plot(medialPoints_2D,pch=20,col="blue",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "") 
  }
  
  minVal<-min(medialPoints_2D[,1])
  maxVal<-max(medialPoints_2D[,1])
  
  # tail to head
  xValues<-seq(minVal, maxVal, length = numberOfSpanialPoints)
  # #head to tail
  # xValues<-seq(maxVal, minVal, length = numberOfSpanialPoints) #head to tail
  
  xData2<-data.frame(x=xValues)
  curvePoints2<-predict(fit4_2D, newdata = xData2)
  
  medialPointsEqualyDis_2D<-cbind(xValues,curvePoints2)
  
  #rotate boundary
  boundaryPoints1<-(boundaryPoints1-matrix(rep(centroid2D,dim(boundaryPoints1)[1]),
                                           ncol = 2,byrow = TRUE))%*%rotationMatrix2D
  boundaryPoints2<-(boundaryPoints2-matrix(rep(centroid2D,dim(boundaryPoints2)[1]),
                                           ncol = 2,byrow = TRUE))%*%rotationMatrix2D
  
  #plot 2D
  if(plotting==TRUE){
    for(i in 1:dim(boundaryPoints1)[1]){
      plot(rbind(boundaryPoints1[i,],boundaryPoints2[i,]),type="l",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "")
      par(new=TRUE)
    }
    plot(medialPointsEqualyDis_2D,pch=20,col="blue",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "")
  }
  
  # find the average of k closet boundary points
  allBoundaryPoints2D<-rbind(boundaryPoints1,boundaryPoints2)
  medialPoints2D_Dis2Boundary<-array(NA,dim = c(dim(medialPointsEqualyDis_2D)[1],dim(allBoundaryPoints2D)[1]))
  for (i in 1:dim(medialPointsEqualyDis_2D)[1]) {
    for (j in 1:dim(allBoundaryPoints2D)[1]) {
      medialPoints2D_Dis2Boundary[i,j]<-norm(medialPointsEqualyDis_2D[i,]-allBoundaryPoints2D[j,],type = "2")
    }
  }
  
  # choose the number of closest points
  k<-k_closestPoint2D
  averageOf_K_ClosestDistance<-rep(NA,dim(medialPointsEqualyDis_2D)[1])
  for (i in 1:dim(medialPointsEqualyDis_2D)[1]) {
    tempArray<-sort(medialPoints2D_Dis2Boundary[i,])
    averageOf_K_ClosestDistance[i]<-mean(tempArray[1:k])
  }
  
  radiiAverage2D<-averageOf_K_ClosestDistance
  
  #plot 2D
  if(plotting==TRUE){
    for(i in 1:dim(boundaryPoints1)[1]){
      plot(rbind(boundaryPoints1[i,],boundaryPoints2[i,]),type="l",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "")
      par(new=TRUE)
    }
    plot(allBoundaryPoints2D,xlim = plotlim,ylim = plotlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(medialPointsEqualyDis_2D,pch=20,col="blue",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "")
    for (i in 1:dim(medialPointsEqualyDis_2D)[1]) {
      draw.circle(x = medialPointsEqualyDis_2D[i,1],
                  y = medialPointsEqualyDis_2D[i,2],
                  radius = averageOf_K_ClosestDistance[i]) 
    } 
  }
  
  regeneratedInternalPoints<-c()
  for (i in 1:dim(medialPointsEqualyDis_2D)[1]) {
    a<-averageOf_K_ClosestDistance[i]
    b<-averageOf_K_ClosestDistance[i]
    regeneratedInternalPoints<-rbind(regeneratedInternalPoints,
                                     ellipsoidGenerator_2D(center = medialPointsEqualyDis_2D[i,]
                                                           ,a = a,b = b, 10, 1))
  }
  #plot 2D
  if(plotting==TRUE){
    plotshapes(regeneratedInternalPoints) 
  }
  
  # Alpha-convex hull
  ahull.obj <- ahull(regeneratedInternalPoints, alpha = alpha2)
  
  #plot 2D
  if(plotting==TRUE){
    plot(ahull.obj,xlim = plotlim,ylim = plotlim) 
  }
  
  indices<-ahull.obj$ashape.obj$alpha.extremes
  
  #plot 2D
  if(plotting==TRUE){
    plot(regeneratedInternalPoints) # all points
  }
    
  #plot 2D
  if(plotting==TRUE){
    plotshapes(regeneratedInternalPoints[indices,]) #boundary points 
  }
  
  #boundaryPoints1 is connected to boundaryPoints2 to form the edges
  boundaryPoints1<-ahull.obj$ashape.obj$edges[,3:4]
  boundaryPoints2<-ahull.obj$ashape.obj$edges[,5:6]
  
  #plot 2D
  if(plotting==TRUE){
    for(i in 1:dim(boundaryPoints1)[1]){
      plot(rbind(boundaryPoints1[i,],boundaryPoints2[i,]),type="l",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "")
      par(new=TRUE)
    }
    plot(medialPointsEqualyDis_2D,pch=20,col="blue",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "") 
  }
  
  
  #remesh the 2D boundary by increasing the number of points on the boundary
  new2D_mesh<-c()
  for (i in 1:dim(boundaryPoints1)[1]) {
    tempPoints<-generatePointsBetween2Points(boundaryPoints1[i,],
                                             boundaryPoints2[i,],
                                             numberOfPoints = numberOfPoints2Dremesh)
    new2D_mesh<-rbind(new2D_mesh,tempPoints)
  }
  
  #plot 2D
  if(plotting==TRUE){
    plot(new2D_mesh,col="black",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(medialPointsEqualyDis_2D,col="blue",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "")
    
  }

  #Test whether the fold is inside the mesh or no
  new2D_meshTest<-new2D_mesh%*%t(rotationMatrix2D)+matrix(rep(centroid2D,dim(new2D_mesh)[1]),ncol = 2,byrow = TRUE)
  edgeData_2d<-data.frame(x=new2D_meshTest[,1], y=new2D_meshTest[,2])
  surfacePointsFold<-predict(fit4, newdata = edgeData_2d)
  fold_3D<-cbind(new2D_meshTest[,1],new2D_meshTest[,2],surfacePointsFold)
  
  foldPointsAreInMesh<-pip3d(Vertices = t(tmesh$vb)[,1:3], Faces = t(tmesh$it),Queries = fold_3D)
  
  #plot  
  if(plotting==TRUE){
    open3d()
    plot3d(fold_3D,type="p",col = "red",expand = 10,box=FALSE,add = TRUE)
    shade3d(tmesh,col="white",alpha=0.2) 
  }
  
  if(sum(foldPointsAreInMesh==-1)>0){
    stop("Fold is not inside the mesh! \n Please change the parameters.")
  }

  ####################
  # generate spokes 2D
  
  xData4<-data.frame(x=new2D_mesh[,1])
  curvePointsY4<-predict(fit4_2D, newdata = xData4)
  medialPointstest<-cbind(new2D_mesh[,1],curvePointsY4)
  
  #plot 2D
  if(plotting==TRUE){
    plot(new2D_mesh,col="black",pch=20,cex=0.2,xlim = plotlim,ylim = plotlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(medialPointstest,col="green",pch=20,cex=0.2,xlim = plotlim,ylim = plotlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(medialPointsEqualyDis_2D,col="red",pch=20,cex=0.2,xlim = plotlim,ylim = plotlim,xlab = "",ylab = "")
  }

  topPart<-new2D_mesh[new2D_mesh[,2]>curvePointsY4,]
  bottomPart<-new2D_mesh[new2D_mesh[,2]<curvePointsY4,]
  
  #plot 2D
  if(plotting==TRUE){
    plot(topPart,col="blue",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(bottomPart,col="red",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "") 
  }
  
  #choose the number Of closet points
  N<-2*round((dim(topPart)[1]/dim(medialPointsEqualyDis_2D)[1])/2)
  # find spokes
  tipOfTopSpokes<-array(NA,dim = dim(medialPointsEqualyDis_2D))
  n_ClosestPointsTop<-array(NA,dim = c(N,2,dim(medialPointsEqualyDis_2D)[1]))
  for (i in 1:dim(medialPointsEqualyDis_2D)[1]) {
    tempDis<-rep(Inf,dim(topPart)[1])
    for (j in 1:dim(topPart)[1]) {
      tempDis[j]<-norm(medialPointsEqualyDis_2D[i,]-topPart[j,],type = "2")
    }
    n_ClosestPoints<-topPart[which.minn(tempDis,n = N),]
    n_ClosestPointsTop[,,i]<-n_ClosestPoints
    tempDirections<-array(NA,dim = c(N,dim(topPart)[2]))
    for (k in 1:N) {
      tempDirections[k,]<-convertVec2unitVec(n_ClosestPoints[k,]-medialPointsEqualyDis_2D[i,])
    }
    
    tipOfTopSpokes[i,]<-medialPointsEqualyDis_2D[i,]+radiiAverage2D[i]*frechetMean(t(tempDirections))
  }
  
  #choose the number Of closet points
  N<-2*round((dim(bottomPart)[1]/dim(medialPointsEqualyDis_2D)[1])/2)
  # find spokes
  tipOfBottomSpokes<-array(NA,dim = dim(medialPointsEqualyDis_2D))
  n_ClosestPointsBottom<-array(NA,dim = c(N,2,dim(medialPointsEqualyDis_2D)[1]))
  for (i in 1:dim(medialPointsEqualyDis_2D)[1]) {
    tempDis<-rep(Inf,dim(bottomPart)[1])
    for (j in 1:dim(bottomPart)[1]) {
      tempDis[j]<-norm(medialPointsEqualyDis_2D[i,]-bottomPart[j,],type = "2")
    }
    n_ClosestPoints<-bottomPart[which.minn(tempDis,n = N),]
    n_ClosestPointsBottom[,,i]<-n_ClosestPoints
    tempDirections<-array(NA,dim = c(N,dim(bottomPart)[2]))
    for (k in 1:N) {
      tempDirections[k,]<-convertVec2unitVec(n_ClosestPoints[k,]-medialPointsEqualyDis_2D[i,])
    }
    
    tipOfBottomSpokes[i,]<-medialPointsEqualyDis_2D[i,]+radiiAverage2D[i]*frechetMean(t(tempDirections))
  }
  
  #plot 2D
  if(plotSubregions==TRUE){
    plot(topPart,pch=20,col="blue",cex=0.2,xlim = plotlim,ylim = plotlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(bottomPart,pch=20,col="red",cex=0.2,xlim = plotlim,ylim = plotlim,xlab = "",ylab = "")
    par(new=TRUE)
    for (i in 1:dim(medialPointsEqualyDis_2D)[1]) {
      for (j in 1:dim(n_ClosestPointsTop)[1]) {
        plot(rbind(medialPointsEqualyDis_2D[i,],n_ClosestPointsTop[j,,i]),type = "l",col="blue",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "") 
        par(new=TRUE)
      }
      for (k in 1:dim(n_ClosestPointsBottom)[1]) {
        plot(rbind(medialPointsEqualyDis_2D[i,],n_ClosestPointsBottom[k,,i]),type = "l",col="red",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "") 
        par(new=TRUE)
      }
    }
  }
  
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
  
  
  tipOfTopSpokes<-tipOfCuttedTopSpokes2D
  tipOfBottomSpokes<-tipOfCuttedBottomSpokes2D
  tipOfExtraSpokeHead<-tipOfCuttedHeadAndTailSpokes2D[1,]
  tipOfExtraSpokeTail<-tipOfCuttedHeadAndTailSpokes2D[2,]
  
  #plot 2D
  if(plotting==TRUE){
    plot(topPart,pch=20,col="blue",cex=0.2,xlim = plotlim,ylim = plotlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(bottomPart,pch=20,col="red",cex=0.2,xlim = plotlim,ylim = plotlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(medialPointsEqualyDis_2D,pch=20,cex=0.2,col="black",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "")
    for (i in 1:dim(medialPointsEqualyDis_2D)[1]) {
      par(new=TRUE)
      plot(rbind(medialPointsEqualyDis_2D[i,],tipOfTopSpokes[i,]),type = "l",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "") 
      par(new=TRUE)
      plot(rbind(medialPointsEqualyDis_2D[i,],tipOfBottomSpokes[i,]),type = "l",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "") 
    }
    par(new=TRUE)
    plot(rbind(tipOfExtraSpokeHead,medialPointsEqualyDis_2D[1,]),type = "l",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "") 
    par(new=TRUE)
    plot(rbind(tipOfExtraSpokeTail,medialPointsEqualyDis_2D[numberOfSpanialPoints,]),type = "l",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "") 
  }
  
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
  
  
  #plot 2D
  if(plotting==TRUE){
    plot(topSpokesPoints,pch=20,cex=0.2,col="blue",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(bottomSpokesPoints,pch=20,cex=0.2,col="red",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(extraSpokeHeadPoints,pch=20,cex=0.2,col="black",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(extraSpokeTailPoints,pch=20,cex=0.2,col="black",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "")
  }
  
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
  
  numberOfLayers<-2*numberOf2DspokePoints-1
  matrixParentsChildren<-c()
  for (i in 1:numberOfLayers) {
    matrixParentsChildren<-rbind(matrixParentsChildren,
                                 seq(i,dim(skeletalPoints_2D)[1],numberOfLayers))
  }
  
  
  #add two extra head and tail points
  skeletalPoints_2D<-rbind(extraSpokeHeadPoints[-1,],
                           skeletalPoints_2D,
                           extraSpokeTailPoints[-1,])
  
  matrixParentsChildren<-matrixParentsChildren+3
  tempPointsNumber<-max(matrixParentsChildren)+1
  matrixParentsChildren<-cbind(c(1,1:(numberOf2DspokePoints-1),1:(numberOf2DspokePoints-1)),
                               matrixParentsChildren,
                               c(tempPointsNumber,tempPointsNumber:(tempPointsNumber+(numberOf2DspokePoints-2)),
                                 tempPointsNumber:(tempPointsNumber+(numberOf2DspokePoints-2))))
  
  #plot 2D
  if(plotting==TRUE){
    plot(skeletalPoints_2D,type = "p",pch=20,col="blue",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "")   
  }
  
  #transfer back
  skeletalPoints_2D<-skeletalPoints_2D%*%t(rotationMatrix2D)+matrix(rep(centroid2D,dim(skeletalPoints_2D)[1]),ncol = 2,byrow = TRUE)
  
  #plot 2D
  if(plotting==TRUE){
    plot(skeletalPoints_2D,type = "p",pch=20,col="blue",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "") 
    par(new=TRUE)
    for (i in 1:numberOfLayers) {
      plot(skeletalPoints_2D[matrixParentsChildren[i,],],type = "l",pch=20,col="blue",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "") 
      par(new=TRUE)
    }  
  }
  
  newData_2d<-data.frame(x=skeletalPoints_2D[,1], y=skeletalPoints_2D[,2])
  surfacePointsZ2<-predict(fit4, newdata = newData_2d)
  medialPoints3D<-cbind(skeletalPoints_2D[,1],skeletalPoints_2D[,2],surfacePointsZ2)
  
  numberOfFrames<-dim(medialPoints3D)[1]

  #draw object fold
  edgeData_2d<-data.frame(x=new2D_mesh[,1], y=new2D_mesh[,2])
  surfacePointsFold<-predict(fit4, newdata = edgeData_2d)

  #plot
  if(plotting==TRUE){
    open3d()
    plot3d(medialPoints3D,type="s",radius = 0.2,col = "black",expand = 10,box=FALSE,add = TRUE)
    plot3d(medialPoints3D[seq(numberOf2DspokePoints,numberOfFrames,numberOfLayers)[1:numberOfSpanialPoints],],type="l",lwd = 2,col = "blue",expand = 10,box=FALSE,add = TRUE)
    for (i in 1:numberOfLayers) {
      plot3d(medialPoints3D[matrixParentsChildren[i,],],type="l",lwd = 1,col = "blue",expand = 10,box=FALSE,add = TRUE)
    }
    # wire3d(tmesh, col="blue")  #wire mesh
    shade3d(tmesh, col="white",alpha=0.2)
  }
  
  #####################
  # 3D m-rep

  # find the average of k closet boundary points
  medialPoints3D_Dis2Boundary<-array(NA,dim = c(numberOfFrames,dim(t(tmesh$vb[1:3,]))[1]))
  pb <- txtProgressBar(style = 3)
  for (i in 1:numberOfFrames) {
    setTxtProgressBar(pb,i/numberOfFrames)
    
    tempMatrix1<-matrix(rep(medialPoints3D[i,],dim(t(tmesh$vb)[,1:3])[1]),ncol = 3,byrow = T)
    tempMatrix2<-tempMatrix1-t(tmesh$vb)[,1:3]
    tempDis<-apply(tempMatrix2, 1, myNorm)
    
    medialPoints3D_Dis2Boundary[i,]<-tempDis
  }
  close(pb)
  
  # choose the number of closest points
  k<-k_closestPoint3D
  averageOf_K_ClosestDistance3D<-rep(NA,numberOfFrames)
  for (i in 1:numberOfFrames) {
    tempArray<-sort(medialPoints3D_Dis2Boundary[i,])
    averageOf_K_ClosestDistance3D[i]<-mean(tempArray[1:k])
  }

  radiiAverage3D<-averageOf_K_ClosestDistance3D

  #plot
  if(plotting==TRUE){
    open3d()
    plot3d(medialPoints3D,type="s",radius = 0.2,col = "blue",expand = 10,box=FALSE,add = TRUE)
    for (i in 1:numberOfFrames) {
      spheres3d(medialPoints3D[i,], radius = radiiAverage3D[i],col = "lightblue")
    }
    shade3d(tmesh, col="grey",alpha=0.4)
  }

  # #save an stl file. It will save all the objects that are in the open3d() windows
  # filename <- tempfile(fileext = ".stl")
  # writeSTL(filename)
  # filename #see the file location
  # #read the stl file
  # readSTL(filename, col = "red")


  # increase the number of 3D boundary points by reduncing the
  # voxelSize of remeshing if it is necessary
  if(increase3DBoundaryPoint==TRUE){
    tmesh<-vcgUniformRemesh(tmesh,voxelSize = 0.5)
    if(plotting==TRUE){
      open3d()
      wire3d(tmesh, col="blue")
    }
  }

  xData5<-data.frame(x=t(tmesh$vb)[,1],y=t(tmesh$vb)[,2])
  surfacePointsZ5<-predict(fit4, newdata = xData5)

  new3D_mesh<-cbind(t(tmesh$vb)[,1:2],surfacePointsZ5)

  #plot
  if(plotting==TRUE){
    open3d()
    shade3d(tmesh, col="white",alpha=0.2)
    plot3d(cbind(t(tmesh$vb)[,1:2],surfacePointsZ5),type="p",
           col = "green",expand = 10,box=FALSE,add = TRUE)
    }

  topPartMesh<-t(tmesh$vb)[,1:3][t(tmesh$vb)[,3]>surfacePointsZ5,]
  bottomPartMesh<-t(tmesh$vb)[,1:3][t(tmesh$vb)[,3]<surfacePointsZ5,]

  #plot
  if(plotting==TRUE){
    open3d()
    plot3d(topPartMesh,type="p",radius = 0.2,col = "blue",expand = 10,box=FALSE,add = TRUE)
    plot3d(bottomPartMesh,type="p",radius = 0.2,col = "red",expand = 10,box=FALSE,add = TRUE)
    }

  #choose the number Of closet points
  N<-7*round((dim(topPartMesh)[1]/numberOfFrames))
  tipOfUpSpokes3D<-array(NA,dim = dim(medialPoints3D))
  pb <- txtProgressBar(style = 3)
  for (i in 1:numberOfFrames) {
    setTxtProgressBar(pb,i/numberOfFrames)
    tempDis<-rep(Inf,dim(topPartMesh)[1])

    tempMatrix1<-matrix(rep(medialPoints3D[i,],dim(topPartMesh)[1]),ncol = 3,byrow = T)
    tempMatrix2<-tempMatrix1-topPartMesh
    tempDis<-apply(tempMatrix2, 1, myNorm)

    n_ClosestPoints<-topPartMesh[which.minn(tempDis,n = N),]

    tempDirections<-array(NA,dim = c(N,3))
    for (k in 1:N) {
      tempDirections[k,]<-convertVec2unitVec(n_ClosestPoints[k,]-medialPoints3D[i,])
    }
    tipOfUpSpokes3D[i,]<-medialPoints3D[i,]+radiiAverage3D[i]*frechetMean(t(tempDirections))
  }
  close(pb)
  
  #choose the number Of closet points
  N<-7*round((dim(bottomPartMesh)[1]/numberOfFrames))
  tipOfDownSpokes3D<-array(NA,dim = dim(medialPoints3D))
  pb <- txtProgressBar(style = 3)
  for (i in 1:numberOfFrames) {
    setTxtProgressBar(pb,i/numberOfFrames)
    tempDis<-rep(Inf,dim(bottomPartMesh)[1])

    tempMatrix1<-matrix(rep(medialPoints3D[i,],dim(bottomPartMesh)[1]),ncol = 3,byrow = T)
    tempMatrix2<-tempMatrix1-bottomPartMesh
    tempDis<-apply(tempMatrix2, 1, myNorm)


    n_ClosestPoints<-bottomPartMesh[which.minn(tempDis,n = N),]
    tempDirections<-array(NA,dim = c(N,3))
    for (k in 1:N) {
      tempDirections[k,]<-convertVec2unitVec(n_ClosestPoints[k,]-medialPoints3D[i,])
    }

    tipOfDownSpokes3D[i,]<- medialPoints3D[i,]+radiiAverage3D[i]*frechetMean(t(tempDirections))

  }
  close(pb)

  tempMesh4Up<-as.mesh3d(medialPoints3D)
  tempMesh4Up$normals<-rbind(apply(tipOfUpSpokes3D-medialPoints3D,FUN = convertVec2unitVec,MARGIN = 1),
                             rep(1,dim(medialPoints3D)[1]))
  tipOfCuttedUpSpokes<-vert2points(vcgRaySearch(tempMesh4Up,mesh = tmesh))

  tempMesh4Down<-as.mesh3d(medialPoints3D)
  tempMesh4Down$normals<-rbind(apply(tipOfDownSpokes3D-medialPoints3D,FUN = convertVec2unitVec,MARGIN = 1),
                               rep(1,dim(medialPoints3D)[1]))
  tipOfCuttedDownSpokes<-vert2points(vcgRaySearch(tempMesh4Down,mesh = tmesh))

  #plot
  if(plotting==TRUE){
    open3d()
    shade3d(tmesh,col="white",alpha=0.2)
    spheres3d(tipOfCuttedUpSpokes,col="blue",radius = 0.2)
    spheres3d(tipOfCuttedDownSpokes,col="red",radius = 0.2)
  }

  #plot s-rep
  if(plotting==TRUE){
    open3d()
    plot3d(medialPoints3D,type="s",radius = 0.1,col = "blue",expand = 10,box=FALSE,add = TRUE)
    for (i in 1:numberOfFrames) {
      plot3d(rbind(medialPoints3D[i,],tipOfCuttedUpSpokes[i,]),type="l",radius = 0.2,col = "blue",expand = 10,box=FALSE,add = TRUE)
    }
    for (i in 1:numberOfFrames) {
      plot3d(rbind(medialPoints3D[i,],tipOfCuttedDownSpokes[i,]),type="l",radius = 0.2,col = "red",expand = 10,box=FALSE,add = TRUE)
    }
    #plot mesh
    # wire3d(tmesh, col="blue")  #wire mesh
    shade3d(tmesh, col="white",alpha=0.2)
  }
  
  result<-list("medialPoints3D"=medialPoints3D,
               "tipOfCuttedUpSpokes"=tipOfCuttedUpSpokes,
               "tipOfCuttedDownSpokes"=tipOfCuttedDownSpokes)
  
  return(result)
}




