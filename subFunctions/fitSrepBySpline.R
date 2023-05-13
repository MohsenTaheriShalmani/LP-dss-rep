
fitSrepBySpline <- function(tmesh=tmesh,
                            plotting=TRUE,
                            plotSubregions=TRUE,
                            fitByRotationOrPCA="Rotation",
                            numberOf2DspokePoints=numberOf2DspokePoints,
                            numberOfSpanialPoints=numberOfSpanialPoints,
                            alpha1=alpha1,
                            alpha2=alpha2,
                            alpha3dSmoother=alpha3dSmoother,
                            k_closestPoint2D,
                            k_closestPoint3D,
                            threshold4EachLineLength=threshold4EachLineLength,
                            thresholdAngle=thresholdAngle,
                            threshold4EachLineLength2D=threshold4EachLineLength2D,
                            thresholdAngle2D=thresholdAngle2D,
                            polyDegree3D=polyDegree3D,
                            polyDegree2D=polyDegree2D,
                            rotationGap=10,
                            numberOfPoints2Dremesh=numberOfPoints2Dremesh,
                            typeOf3Dsmoother="smoothingByInscribedSpheres",
                            increase3DBoundaryPoint=TRUE,
                            voxelSize=0.5){
  
  
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
  
  allGeodesicDistances<-acos(pmin(pmax(as.numeric(rowSums(t(tmesh$normals[1:3,])*t(intersections$normals[1:3,]))),-1.0),1.0))
  
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
                             allInternalLineLengths/2<threshold4EachLineLength*allClosestDistances2Boundary)
  # length(selectedTriangles)
  
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
  
  indices<-ahull.obj$alpha.extremes
  # ahull.obj$ashape.obj$edges
  
  #plot
  if(plotting==TRUE){
    plot(projectedOnXYplane) # all points
    plotshapes(projectedOnXYplane[indices,]) #boundary points 
  }
  
  #boundaryPoints1 is connected to boundaryPoints2 to form the edges
  boundaryPoints1<-ahull.obj$edges[,3:4]
  boundaryPoints2<-ahull.obj$edges[,5:6]
  
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
                                allInternalLineLengths2D/2<threshold4EachLineLength2D*
                                allClosestDistances2Boundary2D)
  
  #plot 2D
  if(plotting==TRUE){
    for(i in 1:dim(boundaryPoints1)[1]){
      plot(rbind(boundaryPoints1[i,],boundaryPoints2[i,]),type="l",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "")
      par(new=TRUE)
    }
    par(new=TRUE)
    plot(allMiddlepoints2D[selectedLineSegments,],pch=20,col="blue",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "")
  }
  
  
  centroid2D<-colMeans(allMiddlepoints2D[selectedLineSegments,])
  centeredPoints2D<-scale(allMiddlepoints2D[selectedLineSegments,],center = TRUE,scale = FALSE)
  
  if(fitByRotationOrPCA=="PCA"){
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
    
    # x y z
    x<-selectedMiddlePoint2D[,1]
    y<-selectedMiddlePoint2D[,2]
    
    fit4_2D <- lm(y ~ poly(x, degree = polyDegree2D,raw = TRUE),
                  data=as.data.frame(cbind(x,y))) #NB!!! raw=F calculate orthogonal polynomial
    
    xData1<-data.frame(x=x)
    curvePoints1<-predict(fit4_2D, newdata = xData1)
    
    medialPoints_2D<-cbind(x,curvePoints1)
    
  }else if(fitByRotationOrPCA=="Rotation"){
    
    #plot2d
    if(plotting==TRUE){
      boundaryPointsTemp1<-(boundaryPoints1-matrix(rep(centroid2D,dim(boundaryPoints1)[1]),
                                                   ncol = 2,byrow = TRUE))
      boundaryPointsTemp2<-(boundaryPoints2-matrix(rep(centroid2D,dim(boundaryPoints2)[1]),
                                                   ncol = 2,byrow = TRUE))
      for (i in 1:dim(boundaryPointsTemp1)[1]) {
        plot(rbind(boundaryPointsTemp1[i,],boundaryPointsTemp2[i,]),type="l",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "") 
        par(new=TRUE)
      }
      par(new=TRUE)
      plot(centeredPoints2D,pch=20,col="blue",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "")
    }
    
    
    if(plotting==TRUE){
      x4plot<-seq(from = plotlim[1],to = plotlim[2],length.out = 100) 
    }
    sumOfAllResiduals<-c()
    for (k in seq(0,360,rotationGap)) {
      tempTheta<-k*2*pi/360
      # cat("Rotation angle theta:",tempTheta,"\n")
      planeRotationMatrix<-matrix(c(cos(tempTheta),-sin(tempTheta),sin(tempTheta),cos(tempTheta)),ncol = 2,byrow = T)
      centeredPoints2DTemp<-centeredPoints2D%*%planeRotationMatrix
      
      # par(new=T)
      # plot(centeredPoints2DTemp,col="green",xlim = plotlim,ylim = plotlim, xlab = "",ylab = "")
      
      x<-centeredPoints2DTemp[,1]
      y<-centeredPoints2DTemp[,2]
      fit <- lm(y ~ poly(x, degree = polyDegree2D,raw = TRUE),
                data=as.data.frame(cbind(x,y))) #NB!!! raw=F calculate orthogonal polynomial
      
      if(plotting==TRUE){
        boundaryPointsTemp1<-(boundaryPoints1-matrix(rep(centroid2D,dim(boundaryPoints1)[1]),
                                                     ncol = 2,byrow = TRUE))%*%planeRotationMatrix
        boundaryPointsTemp2<-(boundaryPoints2-matrix(rep(centroid2D,dim(boundaryPoints2)[1]),
                                                     ncol = 2,byrow = TRUE))%*%planeRotationMatrix
        for (i in 1:dim(boundaryPointsTemp1)[1]) {
          plot(rbind(boundaryPointsTemp1[i,],boundaryPointsTemp2[i,]),type="l",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "") 
          par(new=TRUE) 
        }
        plot(centeredPoints2DTemp,pch=20,col="blue",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "") 
        par(new=TRUE)
        xDataTemp<-data.frame(x=x4plot)
        y4Plot<-predict(fit, newdata = xDataTemp)
        par(new=TRUE)
        plot(cbind(x4plot,y4Plot),type="l",lwd=2,col="green",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "") 
      }
     
      sumOfAllResiduals<-c(sumOfAllResiduals,sum(fit$residuals^2))
      
    }
    
    k<-which.min(sumOfAllResiduals)
    tempTheta<-seq(0,360,rotationGap)[k]*2*pi/360   
    rotationMatrix2D<-matrix(c(cos(tempTheta),-sin(tempTheta),sin(tempTheta),cos(tempTheta)),ncol = 2,byrow = T)
    selectedMiddlePoint2D<-centeredPoints2D%*%rotationMatrix2D
    # x y z
    x<-selectedMiddlePoint2D[,1]
    y<-selectedMiddlePoint2D[,2]
    
    fit4_2D <- lm(y ~ poly(x, degree = polyDegree2D,raw = TRUE),
                  data=as.data.frame(cbind(x,y))) #NB!!! raw=F calculate orthogonal polynomial
    
    xData1<-data.frame(x=x)
    curvePoints1<-predict(fit4_2D, newdata = xData1)
    
    medialPoints_2D<-cbind(x,curvePoints1)
    
  }else{
    stop("Please choose method of fitting as PCA or Rotation!")
  }
  
  #plot 2D
  if(plotting==TRUE){
    plot(selectedMiddlePoint2D,pch=20,col="blue",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "") 
  }
  
  
  #plot 2D
  if(plotting==TRUE){
    plot(medialPoints_2D,pch=20,col="blue",xlim = plotlim,ylim = plotlim,xlab = "",ylab = "") 
  }
  
  startPointX<-min(medialPoints_2D[,1])
  endPointX<-max(medialPoints_2D[,1])
  
  # tail to head
  allCloseX<-seq(startPointX,endPointX,length.out = 10000)
  
  #function of the polynomial
  f_fit <- function(x) { 
    result<-0
    for (i in polyDegree2D:0) {
      result<-result+fit4_2D$coefficients[i+1]*x^i
    }
    return(result)
  }
  #distance function to calculate the distance between two points on the polynomial
  distanceFuncfit <- function(x) { 
    temp<-0
    for (i in polyDegree2D:1) {
      temp<-temp+i*fit4_2D$coefficients[i+1]*x^(i-1)
    }
    return(sqrt(1+temp^2))
  }
  #distance between two points on a polynomial
  distanceBetween2PointsOnCurve <- function(xStartAndxEnd) {
    return(integrate(f = distanceFuncfit,lower = xStartAndxEnd[1],upper = xStartAndxEnd[2])$value)
  }
  
  totalLength<-integrate(f = distanceFuncfit,lower = startPointX,upper = endPointX)$value
  
  tempMatrixLimit<-cbind(rep(startPointX,length(allCloseX)),allCloseX)
  allDistances2StartingPoint<-apply(tempMatrixLimit, MARGIN = 1,FUN = distanceBetween2PointsOnCurve)
  
  equalLengths<-seq(0,totalLength,length.out = numberOfSpanialPoints)
  
  selectedX<-rep(NA,numberOfSpanialPoints)
  for (i in 1:numberOfSpanialPoints) {
    tempIndex<-which.min((allDistances2StartingPoint-equalLengths[i])^2)
    selectedX[i]<-allCloseX[tempIndex]
  }
  
  xData2<-data.frame(x=selectedX)
  curvePoints2<-predict(fit4_2D, newdata = xData2)
  
  medialPointsEqualyDis_2D<-cbind(selectedX,curvePoints2)
  
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
  N<-2*round((dim(topPart)[1]/numberOfSpanialPoints))
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
  N<-2*round((dim(bottomPart)[1]/numberOfSpanialPoints))
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
  
  allSphereMeshes<-vector(mode = "list",length = numberOfFrames)
  for (i in 1:numberOfFrames) {
    allSphereMeshes[[i]]<-makeSphereMesh(center = medialPoints3D[i,],radius = radiiAverage3D[i])
  }
  
  #plot
  if(plotting==TRUE){
    open3d()
    for (i in 1:numberOfFrames) {
      shade3d(allSphereMeshes[[i]],col = "lightblue")
    }
  }
  
  # #plot
  # if(plotting==TRUE){
  #   open3d()
  #   plot3d(medialPoints3D,type="s",radius = 0.2,col = "blue",expand = 10,box=FALSE,add = TRUE)
  #   for (i in 1:numberOfFrames) {
  #     spheres3d(medialPoints3D[i,], radius = radiiAverage3D[i],col = "lightblue")
  #   }
  #   shade3d(tmesh, col="grey",alpha=0.4)
  # }
  
  # #save an stl file. It will save all the objects that are in the open3d() windows
  # filename <- tempfile(fileext = ".stl")
  # writeSTL(filename)
  # filename #see the file location
  # #read the stl file
  # readSTL(filename, col = "red")

  if(typeOf3Dsmoother=="smoothingByInscribedSpheres"){
    
    # allOutsidePoints<-vector(mode = "list",length = numberOfFrames)
    # pb <- txtProgressBar(style = 3)
    # for (k in 1:numberOfFrames) {
    #   setTxtProgressBar(pb,k/numberOfFrames)
    #   
    #   outsidePoints<-1:dim(t(allSphereMeshes[[k]]$vb))[1]
    #   tempMatrix<-matrix(rep(medialPoints3D[k,],numberOfFrames),ncol = 3,byrow = TRUE)-medialPoints3D
    #   sphereDistances<-as.vector(apply(tempMatrix, MARGIN = 1,FUN = myNorm))
    #   sumOfRadii<-radiiAverage3D[k]+radiiAverage3D
    #   adjacentSpheres<-which(sphereDistances<sumOfRadii)
    #   tempIndices<-adjacentSpheres[!adjacentSpheres %in% k] #without considering its own shape
    #   for (i in tempIndices) {
    #     pointsInMesh<-pip3d(Vertices = t(allSphereMeshes[[i]]$vb)[,1:3],
    #                         Faces = t(allSphereMeshes[[i]]$it),
    #                         Queries = t(allSphereMeshes[[k]]$vb)[,1:3])
    #     outsidePoints<-outsidePoints[!outsidePoints %in% which(pointsInMesh==1)]
    #   }
    #   allOutsidePoints[[k]]<-outsidePoints
    # }
    # close(pb)
    # boundaryPointsSpheres<-c()
    # for (i in 1:numberOfFrames) {
    #   boundaryPointsSpheres<-rbind(boundaryPointsSpheres,
    #                                t(allSphereMeshes[[i]]$vb[1:3,allOutsidePoints[[i]]]))
    # }
    # #plot
    # if(plotting==TRUE){
    #   open3d()
    #   plot3d(boundaryPointsSpheres,type="p", col = "black",expand = 10,box=FALSE,add = TRUE)
    # }
    # # alpha 3D shape
    # tempMesh<-as.mesh3d(ashape3d(boundaryPointsSpheres, alpha=alpha3dSmoother,pert = TRUE))
    # shade3d(tempMesh)
    #set of mesh3d objects
    # tempMesh3dObjects<-vector(mode = "list",length = numberOfFrames)
    # for (i in 1:numberOfFrames) {
    #   tempMesh3dObjects[[i]]<-vcgSphere(subdivision = 3, normals = TRUE)
    # }
    # 
    # for (k in 1:numberOfFrames) {
    #   tempMatrix<-t(allSphereMeshes[[k]]$it)
    #   verticesMatrix<-c()
    #   outsidePoints<-allOutsidePoints[[k]]
    #   for (i in 1:length(allOutsidePoints[[k]])) {
    #     verticesMatrix<-c(verticesMatrix,which(tempMatrix[,1]==outsidePoints[i] | 
    #                                              tempMatrix[,2]==outsidePoints[i] | 
    #                                              tempMatrix[,3]==outsidePoints[i] ))
    #     
    #   }
    #   tempMesh3dObjects[[k]]$vb<-allSphereMeshes[[k]]$vb
    #   tempMesh3dObjects[[k]]$it<-t(tempMatrix[verticesMatrix,])
    # }
    # 
    # #plot
    # if(plotting==TRUE){
    #   open3d()
    #   for (i in 1:numberOfFrames) {
    #     shade3d(tempMesh3dObjects[[i]]) 
    #   }
    # }
    # 
    # mergedShapes<-mergeMeshes(tempMesh3dObjects[[1]],tempMesh3dObjects[[2]])
    # for (i in 3:numberOfFrames) {
    #   mergedShapes<-mergeMeshes(mergedShapes,tempMesh3dObjects[[i]])
    # }
    # 
    # #plot
    # if(plotting==TRUE){
    #   open3d()
    #   shade3d(mergedShapes,alpha=0.2) 
    # }
    # 
    # 
    # mergedShapesCleaned<-vcgClean(mergedShapes,sel = c(0,1,6),tol = 0.1)
    # #plot
    # if(plotting==TRUE){
    #   open3d()
    #   shade3d(mergedShapesCleaned,alpha=0.2)
    # }
    # 
    # mergedShapesCleanedSmoothed<-vcgSmooth(mergedShapesCleaned)
    # #plot
    # if(plotting==TRUE){
    #   open3d()
    #   shade3d(mergedShapesCleanedSmoothed,col="blue",alpha=1)
    # }
    # 
    # mergedShapesRemeshed<-vcgUniformRemesh(mergedShapesCleanedSmoothed,offset = 0)
    # #plot
    # if(plotting==TRUE){
    #   open3d()
    #   shade3d(mergedShapesRemeshed,col="blue",alpha=0.5)
    # }
    # mergedShapesRemeshedSmoothed<-vcgSmooth(mergedShapesRemeshed,mu = 0)
    # #plot
    # if(plotting==TRUE){
    #   open3d()
    #   shade3d(mergedShapesRemeshedSmoothed,col="blue",alpha=0.5)
    # }
    # 
    # mergedShapesRemeshed2<-vcgUniformRemesh(mergedShapesRemeshedSmoothed,voxelSize = 0.5)
    # #plot
    # if(plotting==TRUE){
    #   open3d()
    #   shade3d(mergedShapesRemeshed2,col="blue",alpha=1)
    # }
    # xData5<-data.frame(x=t(mergedShapesRemeshed2$vb)[,1],y=t(mergedShapesRemeshed2$vb)[,2])
    # surfacePointsZ5<-predict(fit4, newdata = xData5)
    # 
    # #plot
    # if(plotting==TRUE){
    #   open3d()
    #   shade3d(mergedShapesRemeshed2, col="white",alpha=0.2)
    #   plot3d(cbind(t(mergedShapesRemeshed2$vb)[,1:2],surfacePointsZ5),type="p",
    #          col = "green",expand = 10,box=FALSE,add = TRUE)
    # }
    #plot
    # if(plotting==TRUE){
    #   open3d()
    #   shade3d(mergedShapesRemeshed2, col="white",alpha=0.2)
    #   plot3d(cbind(t(mergedShapesRemeshed2$vb)[,1:2],surfacePointsZ5),type="p",
    #          col = "green",expand = 10,box=FALSE,add = TRUE)
    # }
    # 
    # topPartMesh<-t(mergedShapesRemeshed2$vb)[,1:3][t(mergedShapesRemeshed2$vb)[,3]>surfacePointsZ5,]
    # bottomPartMesh<-t(mergedShapesRemeshed2$vb)[,1:3][t(mergedShapesRemeshed2$vb)[,3]<surfacePointsZ5,]
    
    
    allSpherePoints<-c()
    for (i in 1:numberOfFrames) {
      allSpherePoints<-rbind(allSpherePoints,
                                   t(allSphereMeshes[[i]]$vb[1:3,]))
    }
    #plot
    if(plotting==TRUE){
      open3d()
      plot3d(allSpherePoints,type="p", col = "black",expand = 10,box=FALSE,add = TRUE)
    }
    
    # alpha 3D shape
    cat("Please wait for a minute.")
    alpha_3DShape<-as.mesh3d(ashape3d(allSpherePoints, alpha=alpha3dSmoother,pert = TRUE))
    
    if(plotting==TRUE){
      open3d()
      shade3d(alpha_3DShape)
    }
    
    #remove internal vertices
    alpha_3DShapeCleaned<-vcgClean(alpha_3DShape,sel = c(0,1))
    
    smoothPointsOfAlpha3DShape<-c()
    pb <- txtProgressBar(style = 3)
    for (i in 1:dim(alpha_3DShapeCleaned$it)[2]) {
      setTxtProgressBar(pb,i/numberOfFrames)
      
      p1<-alpha_3DShapeCleaned$vb[1:3,alpha_3DShapeCleaned$it[1,i]]
      p2<-alpha_3DShapeCleaned$vb[1:3,alpha_3DShapeCleaned$it[2,i]]
      p3<-alpha_3DShapeCleaned$vb[1:3,alpha_3DShapeCleaned$it[3,i]]
      
      trianglepoints<-generatePointsBetween3Points(point1 = p1,
                                                   point2 = p2,
                                                   point3 = p3,
                                                   numberOf2DspokePoints = 3)
      
      smoothPointsOfAlpha3DShape<-rbind(smoothPointsOfAlpha3DShape,trianglepoints)
    }
    smoothPointsOfAlpha3DShape<-rbind(t(alpha_3DShapeCleaned$vb[1:3,]),
                                      smoothPointsOfAlpha3DShape)
    close(pb)
    
    if(plotting==TRUE){
      open3d()
      plot3d(smoothPointsOfAlpha3DShape,type="p", col = "black",expand = 10,box=FALSE,add = TRUE)
    }
    
    #regenerate alpha3d shape
    alpha_3DShape2<-as.mesh3d(ashape3d(smoothPointsOfAlpha3DShape, alpha=alpha3dSmoother,pert = TRUE))
    
    #plot
    if(plotting==TRUE){
      open3d()
      shade3d(alpha_3DShape2)
    }
    
    if(plotting==TRUE){
      open3d()
      plot3d(t(alpha_3DShape2$vb[1:3,]),type="p", col = "black",expand = 10,box=FALSE,add = TRUE)
    }
    
    #remove internal vertices
    alpha3DShapeCleaned<-vcgClean(alpha_3DShape2,sel = c(0,1))
    
    #plot
    if(plotting==TRUE){
      open3d()
      shade3d(alpha3DShapeCleaned,col="white",alpha=0.2)
    }
    
    xData5<-data.frame(x=t(alpha3DShapeCleaned$vb)[,1],y=t(alpha3DShapeCleaned$vb)[,2])
    surfacePointsZ5<-predict(fit4, newdata = xData5)
    
    #plot
    if(plotting==TRUE){
      open3d()
      shade3d(alpha3DShapeCleaned, col="white",alpha=0.2)
      plot3d(cbind(t(alpha3DShapeCleaned$vb)[,1:2],surfacePointsZ5),type="p",
             col = "green",expand = 10,box=FALSE,add = TRUE)
    }
    
    topPartMesh<-t(alpha3DShapeCleaned$vb)[,1:3][t(alpha3DShapeCleaned$vb)[,3]>surfacePointsZ5,]
    bottomPartMesh<-t(alpha3DShapeCleaned$vb)[,1:3][t(alpha3DShapeCleaned$vb)[,3]<surfacePointsZ5,]
    
    
  }else if(typeOf3Dsmoother=="smoothingByTheOriginalMesh"){
    # increase the number of 3D boundary points by reduncing the
    # voxelSize of remeshing if it is necessary
    if(increase3DBoundaryPoint==TRUE){
      tmesh<-vcgUniformRemesh(tmesh,voxelSize = voxelSize)
      if(plotting==TRUE){
        open3d()
        wire3d(tmesh, col="blue")
      }
    }
    
    xData5<-data.frame(x=t(tmesh$vb)[,1],y=t(tmesh$vb)[,2])
    surfacePointsZ5<-predict(fit4, newdata = xData5)
    
    #plot
    if(plotting==TRUE){
      open3d()
      shade3d(tmesh, col="white",alpha=0.2)
      plot3d(cbind(t(tmesh$vb)[,1:2],surfacePointsZ5),type="p",
             col = "green",expand = 10,box=FALSE,add = TRUE)
    }
    
    topPartMesh<-t(tmesh$vb)[,1:3][t(tmesh$vb)[,3]>surfacePointsZ5,]
    bottomPartMesh<-t(tmesh$vb)[,1:3][t(tmesh$vb)[,3]<surfacePointsZ5,]
    
  }else{
    stop("Please specify the type Of 3D smoother as smoothingByInscribedSpheres or smoothingByTheOriginalMesh ! ")
  }
  
  #plot
  if(plotting==TRUE){
    open3d()
    plot3d(topPartMesh,type="p",col = "blue",expand = 10,box=FALSE,add = TRUE)
    plot3d(bottomPartMesh,type="p",col = "red",expand = 10,box=FALSE,add = TRUE)
  }
  
  #choose the number Of closet points
  N<-2*round((dim(topPartMesh)[1]/numberOfFrames))
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
  N<-2*round((dim(bottomPartMesh)[1]/numberOfFrames))
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
  
  #to use as.mesh3d for medial points the number of points must be devidable by 3
  numberOfVirtualPoints<-(3-dim(medialPoints3D)[1]%%3)%%3
  medialPoints3D_2<-rbind(medialPoints3D,
                          matrix(rep(colMeans(medialPoints3D),numberOfVirtualPoints),ncol = 3,byrow = TRUE))
  tipOfUpSpokes3D_2<-rbind(tipOfUpSpokes3D,
                           matrix(rep(tipOfUpSpokes3D[1,],numberOfVirtualPoints),ncol = 3,byrow = TRUE))
  tipOfDownSpokes3D_2<-rbind(tipOfDownSpokes3D,
                           matrix(rep(tipOfDownSpokes3D[1,],numberOfVirtualPoints),ncol = 3,byrow = TRUE))
  
  
  tempMesh4Up<-as.mesh3d(medialPoints3D_2,triangles = TRUE)
  tempMesh4Up$normals<-rbind(apply(tipOfUpSpokes3D_2-medialPoints3D_2,FUN = convertVec2unitVec,MARGIN = 1),
                             rep(1,dim(medialPoints3D_2)[1]))
  tipOfCuttedUpSpokesTemp<-vert2points(vcgRaySearch(tempMesh4Up,mesh = tmesh))
  tipOfCuttedUpSpokes<-tipOfCuttedUpSpokesTemp[1:dim(medialPoints3D)[1],]
  
  tempMesh4Down<-tempMesh4Up
  tempMesh4Down$normals<-rbind(apply(tipOfDownSpokes3D_2-medialPoints3D_2,FUN = convertVec2unitVec,MARGIN = 1),
                               rep(1,dim(medialPoints3D_2)[1]))
  tipOfCuttedDownSpokesTemp<-vert2points(vcgRaySearch(tempMesh4Down,mesh = tmesh))
  tipOfCuttedDownSpokes<-tipOfCuttedDownSpokesTemp[1:dim(medialPoints3D)[1],]
  
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
