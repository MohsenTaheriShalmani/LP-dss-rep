fit_LP_dss_rep<- function(tmesh,
                          plotting=TRUE,
                          numberOfCoverPoints3D=100000,
                          k_Ring3D=10, # for dividing 2D and 3D object
                          lambda3D=0.5, # for dividing 2D and 3D object
                          k_Ring2D=5, # for dividing 2D and 3D object
                          lambda2D=2.5, # for dividing 2D and 3D object
                          sphereResolution=1, # choose 1,2,or 3 for the resolution of the urchin
                          circleResolution=24, # resolution of the 2D urchin
                          circleResolutionHigh=150, # resolution of the 2D urchin
                          urchinRadius=0.5, # resolution of the 2D urchin
                          thresholdDistance2D=0.2,
                          polyDegree3D=4,
                          polyDegree2D=4,
                          alpha1=2, # for alpha convex hull
                          numberOfPoints4alphaHull=5000,
                          numberOf2DspokePoints=4,
                          numberOfSpanialPoints=27,
                          numberOfSpokes4Interpolation=5,
                          rotationGap=10, #to fit best medial curve
                          cross_sections_visualization=FALSE) {
  
  #check mesh integrity
  meshintegrity(tmesh,facecheck = TRUE)
  
  #update normals
  tmesh<-vcgUpdateNormals(tmesh)
  
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
  
  
  numberOfVertices<-dim(tmesh$vb)[2]
  geoDistancesMatrix<-array(NA,dim = c(numberOfVertices,numberOfVertices))
  for (i in 1:numberOfVertices) {
    geoDistancesMatrix[i,]<-vcgDijkstra(tmesh,i)
  }

  #make geoDistancesMatrix symmetric
  geoDistancesMatrix<-forceSymmetric(geoDistancesMatrix)
  # geoDistancesMatrix<-upper.tri(geoDistancesMatrix,diag = TRUE)+t(upper.tri(geoDistancesMatrix,diag = FALSE))

  # mean curvature of the vertices
  vcgInfo<-vcgCurve(tmesh)
  meanCurvatureOfAllVertices<-vcgInfo$meanvb

  #adjancy matrix of neighbors
  allAdjacentVertices<-vcgVertexNeighbors(tmesh,numstep = k_Ring3D)

  unitNormals<-t(tmesh$normals[1:3,])/sqrt(rowSums(t(tmesh$normals[1:3,])^2))

  normalMultiplicationMatrix<-unitNormals%*%t(unitNormals)

  elementWiseMultipleMatrix<-exp(lambda3D*normalMultiplicationMatrix*geoDistancesMatrix)

  #compute affinity matrix W
  W<-elementWiseMultipleMatrix
  numberOfPoints<-dim(W)[2]
  tempVec<-1:numberOfPoints
  pb <- txtProgressBar(style = 3)
  for (i in 1:numberOfPoints) {
    setTxtProgressBar(pb,i/numberOfPoints)
    notNeighbor<-tempVec[!tempVec %in% allAdjacentVertices[[i]]]
    W[i,notNeighbor]<-0
  }
  close(pb)

  d<-(rowSums(W))^(-1/2)
  D<-diag(d)

  L <- D%*%W%*%D

  L <- diag(nrow(L))-L

  eigenValAndVec<-eigen(L)

  smallestEigenVectorIndex<-nrow(L)-1
  secondSmallestEigenValue<-eigenValAndVec$values[smallestEigenVectorIndex]
  secondSmallestEigenVector<-eigenValAndVec$vectors[,smallestEigenVectorIndex]

  upIndices<-which(secondSmallestEigenVector>0)
  downIndices<-which(secondSmallestEigenVector<=0)


  #plot
  if(plotting==TRUE){
    open3d()
    spheres3d(vert2points(tmesh)[upIndices,],radius = 2,col='orange')
    spheres3d(vert2points(tmesh)[downIndices,],radius = 2,col='blue')
  }

  # we consider the top and bottom parts of the object as two meshes by
  # tmesh$it which is the face indices
  #Up
  tempMatrix<-matrix(t(tmesh$it) %in% upIndices,ncol = 3,byrow = FALSE)
  polyIndicesUp<-which(rowSums(tempMatrix)==3)
  trglsUp <- as.matrix(t(t(tmesh$it)[polyIndicesUp,]))
  tmeshUp <- tmesh3d(tmesh$vb, trglsUp)
  #Edge
  polyIndicesEdge<-which(rowSums(tempMatrix)==2 | rowSums(tempMatrix)==1)
  trglsEdge <- as.matrix(t(t(tmesh$it)[polyIndicesEdge,]))
  tmeshEdge <- tmesh3d(tmesh$vb, trglsEdge)
  #Down
  tempMatrix<-matrix(t(tmesh$it) %in% downIndices,ncol = 3,byrow = FALSE)
  polyIndicesDown<-which(rowSums(tempMatrix)==3)
  trglsDown <- as.matrix(t(t(tmesh$it)[polyIndicesDown,]))
  tmeshDown <- tmesh3d(tmesh$vb, trglsDown)


  #remove unreferenced vertices
  tmeshUp<-vcgClean(tmeshUp,sel = 1)
  tmeshEdge<-vcgClean(tmeshEdge,sel = 1)
  tmeshDown<-vcgClean(tmeshDown,sel = 1)

  #plot
  if(plotting==TRUE){
    open3d()
    shade3d(tmeshUp,col="blue")
    shade3d(tmeshEdge,col="yellow")
    shade3d(tmeshDown,col="red")
  }
  
  
  ######################################################################################################
  ######################################################################################################
  # fit spline
  
  minTempX<-min(tmesh$vb[1,])-1
  maxTempX<-max(tmesh$vb[1,])+1
  minTempY<-min(tmesh$vb[2,])-1
  maxTempY<-max(tmesh$vb[2,])+1
  minTempZ<-min(tmesh$vb[3,])-1
  maxTempZ<-max(tmesh$vb[3,])+1
  
  xlim<-c(minTempX,maxTempX)
  ylim<-c(minTempY,maxTempY)
  zlim<-c(minTempZ,maxTempZ)
  
  xs<-runif(n = numberOfCoverPoints3D,min =minTempX ,max =maxTempX)
  ys<-runif(n = numberOfCoverPoints3D,min =minTempY ,max =maxTempY)
  zs<-runif(n = numberOfCoverPoints3D,min =minTempZ ,max =maxTempZ)
  
  coverpoint<-as.matrix(cbind(xs,ys,zs)) #as.matrix is necessary
  
  # #plot
  # if(plotting==TRUE){
  #   open3d()
  #   shade3d(tmesh, col="white",alpha=0.5)  #surface mesh
  #   plot3d(coverpoint,type="p",col = "black",expand = 10,box=FALSE,add = TRUE)
  # }
  
  #make the mesh slightly smaller
  smallMesh<-tmesh
  smallMesh$vb[1:3,]<-smallMesh$vb[1:3,]-0.4*smallMesh$normals[1:3,]
  
  #find internal points
  insidePoints<-pip3d(Vertices = t(smallMesh$vb[1:3,]),
                      Faces = t(smallMesh$it),
                      Queries = coverpoint)
  
  #find internal points
  insidePoints<-pip3d(Vertices = t(tmesh$vb[1:3,]),
                      Faces = t(tmesh$it),
                      Queries = coverpoint)
  innerPoints<-coverpoint[insidePoints==1,]
  
  #plot
  if(plotting==TRUE){
    open3d()
    shade3d(tmesh, col="white",alpha=0.5)  #surface mesh
    plot3d(innerPoints,type="p",col = "black",expand = 10,box=FALSE,add = TRUE)
  }
  
  
  pb <- txtProgressBar(style = 3)
  force<-array(NA,nrow(innerPoints))
  shortestLengths<-array(NA,nrow(innerPoints))
  # angles<-array(NA,nrow(innerPoints))
  for (k in 1:nrow(innerPoints)) {
    setTxtProgressBar(pb,k/nrow(innerPoints))

    center<-innerPoints[k,]

    tempSphere<-makeSphereMesh(center = center,radius = urchinRadius,subdivision = sphereResolution)

    # #plot
    # shade3d(tmesh,col="white",alpha=0.1)
    # shade3d(tempSphere,col="blue")

    #find the intercetions of rays with internal normal directions
    intersectionsUp <- vcgRaySearch(tempSphere,mesh = tmeshUp,mindist = TRUE)
    intersectionsDown <- vcgRaySearch(tempSphere,mesh = tmeshDown,mindist = TRUE)

    tipOfCuttedSpokesUp<-vert2points(intersectionsUp)[intersectionsUp$quality==1,]
    tipOfCuttedSpokesDown<-vert2points(intersectionsDown)[intersectionsDown$quality==1,]

    if(!(is.matrix(tipOfCuttedSpokesUp) & is.matrix(tipOfCuttedSpokesDown))){

      force[k]<-10^4
      shortestLengths[k]<-0

      next
    }

    # #plot
    # if(plotting==TRUE){
    #   open3d()
    #   # spheres3d(vert2points(intersections),col="blue",radius = 0.2) #plot intersections
    #   for (i in 1:dim(tipOfCuttedSpokesUp)[1]) {
    #     plot3d(rbind(center,tipOfCuttedSpokesUp[i,]),type="l",col = "blue",expand = 10,box=FALSE,add = TRUE)
    #   }
    #   for (i in 1:dim(tipOfCuttedSpokesDown)[1]) {
    #     plot3d(rbind(center,tipOfCuttedSpokesDown[i,]),type="l",col = "red",expand = 10,box=FALSE,add = TRUE)
    #   }
    #   #NB!!! we have the information of normals at intersections intersections$normals !!!!
    #   shade3d(tmesh, col="white",alpha=0.2)  #surface mech
    # }

    tipOfCuttedSpokes<-rbind(tipOfCuttedSpokesUp,tipOfCuttedSpokesDown)
    labels<-c(rep(1,nrow(tipOfCuttedSpokesUp)),
              rep(2,nrow(tipOfCuttedSpokesDown)))

    vectorsMatrix<-tipOfCuttedSpokes-matrix(rep(center,nrow(tipOfCuttedSpokes)),ncol = 3,byrow = TRUE)

    # spokeLengths<-apply(X = vectorsMatrix, MARGIN = 1,FUN = myNorm )
    # spokesUnitDirections<-t(apply(X = vectorsMatrix, MARGIN = 1,FUN = convertVec2unitVec ))

    forceTemp<-forceFunctionAngleBased_4MultiObject(vectorsMatrix = vectorsMatrix,
                                                    labels = labels,
                                                    type = 'one')

    # forceTemp<-forceFunctionAngleBased(vectorsMatrix = vectorsMatrix)

    force[k]<-forceTemp$force

    shortestLengths[k]<-forceTemp$shortestLength

    # angles[k]<-urchin_AngleOfShortestSpokes(vectorsMatrix = vectorsMatrix)

    # angles[k]<-urchin_AngleOfShortestSpokes_4MultiObject(vectorsMatrix = vectorsMatrix,
    #                                           labels = labels,
    #                                           type = 'one')


  }
  close(pb)

  forceVectorsMagnitudes<-force
  # forceVectorsMagnitudes<-angles
  # plot(1:length(forceVectorsMagnitudes),forceVectorsMagnitudes,type = 'p')

  roundedForce<-round(forceVectorsMagnitudes,digits = 2)
  uniqueRounded<-sort(unique(roundedForce))
  # table(roundedForce)
  # barplot(uniqueRounded,roundedForce)
  
  distanceThreshold<-quantile(shortestLengths,0.1)
  
  # #plot levels of conditional medial points
  # for (i in length(uniqueRounded):1) {
  #   open3d()
  #   shade3d(tmesh,col="white",alpha=0.1)
  #   plot3d(innerPoints[(roundedForce %in% uniqueRounded[1:i]) &
  #                        shortestLengths>distanceThreshold,],type="p",col = "blue",expand = 10,box=FALSE,add = TRUE)
  # }
  
  #conditional medial points
  selectedMiddlePoint<-innerPoints[(roundedForce %in% uniqueRounded[1:(length(uniqueRounded)-1)]) &
                                     shortestLengths>distanceThreshold,]
  
  # #plot conditional medial points L1
  # if(plotting==TRUE){
  #   open3d()
  #   shade3d(tmesh,col="white",alpha=0.1)
  #   plot3d(selectedMiddlePoint,type="p",col = "blue",expand = 10,box=FALSE,add = TRUE)
  # }
  
  
  ######################################################################################################
  ######################################################################################################
  # medial surface
  
  # i<-6
  # 
  # selectedMiddlePoint<-innerPoints[(roundedForce %in% uniqueRounded[1:i]) &
  #                                    shortestLengths>distanceThreshold,]
  # 
  # plot3d(selectedMiddlePoint,type="p",col = "blue",expand = 10,box=FALSE,add = TRUE)
  # plot3d(edgePoints,type="p",col = "blue",expand = 10,box=FALSE,add = TRUE)
  # shade3d(tmesh,col="white",alpha=0.1)
  
  edgePoints<-c()
  for (i in 1:dim(t(tmeshEdge$it))[1]) {
    edgePoints<-rbind(edgePoints,
                      colMeans((t(tmeshEdge$vb)[t(tmeshEdge$it)[i,],1:3])))
  }
  
  # if(plotting==TRUE){
  #   open3d()
  #   # plot3d(selectedMiddlePoint,type="p",col = "blue",expand = 10,box=FALSE,add = TRUE)
  #   plot3d(edgePoints,type="p",col = "red",expand = 10,box=FALSE,add = TRUE)
  #   shade3d(tmesh,col="white",alpha=0.2)
  # }
  
  
  #plot conditional medial points L1
  if(plotting==TRUE){
    open3d()
    shade3d(tmesh,col="white",alpha=0.1)
    plot3d(selectedMiddlePoint,type="p",col = "blue",expand = 10,box=FALSE,add = TRUE)
    plot3d(edgePoints,type="p",col = "red",expand = 10,box=FALSE,add = TRUE)
  }
  
  allSelectedPoints<-rbind(selectedMiddlePoint,edgePoints)
  
  # x y z
  # allSelectedPoints<-vert2points(tmesh)
  
  x<-allSelectedPoints[,1]
  y<-allSelectedPoints[,2]
  z<-allSelectedPoints[,3]
  
  fit4 <- lm(z ~ poly(x, y, degree = polyDegree3D ,raw = TRUE), data=as.data.frame(cbind(z,x,y)))
  # summary(fit2)
  x2<-coverpoint[,1]
  y2<-coverpoint[,2]
  newDATA<-data.frame(x=x2, y=y2)
  surfacePointsZ<-predict(fit4, newdata = newDATA)
  
  medialPoints<-cbind(x2,y2,surfacePointsZ)
  
  insidePoints2<-pip3d(Vertices = t(tmesh$vb[1:3,]),
                       Faces = t(tmesh$it),
                       Queries = medialPoints)
  medialPoints2<-medialPoints[insidePoints2==1,]
  
  #plot
  if(plotting==TRUE){
    # open3d()
    # plot3d(medialPoints,type="p",col = "red",expand = 10,box=FALSE,add = TRUE)
    # shade3d(tmeshUp, col="blue",alpha=1)
    # shade3d(tmeshEdge, col="yellow",alpha=1) 
    # shade3d(tmeshDown, col="blue",alpha=1)
    
    open3d()
    plot3d(medialPoints2,type="p",col = "red",expand = 10,box=FALSE,add = TRUE)
    shade3d(tmesh, col="blue",alpha=0.2) 
  }
  
  
  ######################################################################################################
  ######################################################################################################
  # collapse mesh to measure and plot curvature
  # 
  # if(plotting==TRUE){
  #   
  #   boundaryPoints<-vert2points(tmeshUp)
  #   
  #   medialMeshPoints_4plot<-array(NA,dim = dim(boundaryPoints))
  #   for (i in 1:nrow(boundaryPoints)) {
  #     tempMatrix<-matrix(rep(boundaryPoints[i,],nrow(medialPoints)),ncol = 3,byrow = TRUE)
  #     tempDistances<-apply(medialPoints-tempMatrix, MARGIN = 1,FUN = myNorm)
  #     medialMeshPoints_4plot[i,]<-medialPoints[which.min(tempDistances),]
  #   }
  #   medialMesh3D_4plot<-tmeshUp
  #   vertsTemp <- rbind(t(as.matrix(medialMeshPoints_4plot)),1)
  #   medialMesh3D_4plot$vb<-vertsTemp
  #   
  #   #plot
  #   open3d()
  #   shade3d(tmesh,col='white',alpha=0.2)
  #   shade3d(medialMesh3D_4plot,col='red')
  #   plotNormals(medialMesh3D_4plot)
  #   
  #   #remesh and clean only for curvature plotting
  #   medialMesh3D_4plot<-vcgUniformRemesh(medialMesh3D_4plot,voxelSize = 0.1)
  #   medialMesh3D_4plot<-vcgClean(medialMesh3D_4plot)
  #   medialMesh3D_4plot<-updateNormals(medialMesh3D_4plot)
  #   
  #   curvatures<-vcgCurve(medialMesh3D_4plot)
  #   hist(curvatures$gaussvb,breaks = 100)
  #   
  #   meshDist(medialMesh3D_4plot,distvec=curvatures$meanvb,from=-0.2,to=0.2,tol=0.01)
  #   shade3d(tmesh,col='white',alpha=0.2)
  #   
  # }
  
  
  ######################################################################################################
  ######################################################################################################
  # projection of the medial points
  
  projectedOnXYplane<-cbind(medialPoints2[,1],medialPoints2[,2])
  
  # reduce number of points to boost Alpha-convex hull
  tempLables<-sample(1:nrow(projectedOnXYplane),size = numberOfPoints4alphaHull,replace = FALSE)
  sampledProjectedOnXYplane<-projectedOnXYplane[tempLables,]
  
  #plot
  if(plotting==TRUE){
    plot(sampledProjectedOnXYplane,pch='.',xlim = xlim,ylim = xlim,xlab = "",ylab = "") 
  }
  
  
  # Alpha-convex hull
  ahull.obj <- ahull(sampledProjectedOnXYplane, alpha = alpha1)
  
  #plot
  if(plotting==TRUE){
    plot(ahull.obj,xlim = xlim,ylim = xlim,xlab = "",ylab = "") 
  }
  
  indices<-ahull.obj$ashape.obj$alpha.extremes
  
  #boundaryPoints1 is connected to boundaryPoints2 to form the edges
  boundaryPoints1<-ahull.obj$ashape.obj$edges[,3:4]
  boundaryPoints2<-ahull.obj$ashape.obj$edges[,5:6]
  
  if(plotting==TRUE){
    for(i in 1:dim(boundaryPoints1)[1]){
      plot(rbind(boundaryPoints1[i,],boundaryPoints2[i,]),type="l",xlim = xlim,ylim = xlim,xlab = "",ylab = "")
      par(new=TRUE)
    } 
  }
  
  
  alphaHullIndices1<-ahull.obj$ashape.obj$edges[,1]
  alphaHullIndices2<-ahull.obj$ashape.obj$edges[,2]
  
  tempIndex<-alphaHullIndices1[1]
  alphaHullIndices_sorted<-c(tempIndex)
  while (length(alphaHullIndices1)>0 & length(alphaHullIndices2)>0) {
    a<-which(alphaHullIndices1==tempIndex)[1]
    b<-which(alphaHullIndices2==tempIndex)[1]
    if(!is.na(a)){
      tempIndex<-alphaHullIndices2[a]
      alphaHullIndices_sorted<-c(alphaHullIndices_sorted,tempIndex)
      alphaHullIndices1<-alphaHullIndices1[-a]
      alphaHullIndices2<-alphaHullIndices2[-a]
    }else if(!is.na(b)){
      tempIndex<-alphaHullIndices1[b]
      alphaHullIndices_sorted<-c(alphaHullIndices_sorted,tempIndex)
      alphaHullIndices1<-alphaHullIndices1[-b]
      alphaHullIndices2<-alphaHullIndices2[-b]
    }else{
      break
    }
  }
  
  mesh2DPoints<-ahull.obj$ashape.obj$x[alphaHullIndices_sorted,]
  
  if(plotting==TRUE){
    plot(mesh2DPoints,type="l",xlim = xlim,ylim = xlim,xlab = "",ylab = "") 
  }
  
  # mesh2D must be open
  mesh2D<-mesh2DPoints[-nrow(mesh2DPoints),]
  
  # #remove duplicated vertices
  # duplicatedIndices<-which(duplicated(mesh2D[,1])+duplicated(mesh2D[,2])==2)
  # mesh2D<-mesh2D[-duplicatedIndices,]
  
  if(plotting==TRUE){
    plot(rbind(mesh2D,mesh2D[1,]),type = "l",xlim = xlim,ylim = xlim,xlab = "",ylab = "") 
  }
  
  #reverse the mesh so that it contains the innerpoints
  innerpointsIndices<-which(pip2d(mesh2D,sampledProjectedOnXYplane)==1)
  if(length(innerpointsIndices)<nrow(sampledProjectedOnXYplane)/2){
    mesh2D<-cbind(rev(mesh2D[,1]),rev(mesh2D[,2]))
    innerpointsIndices<-which(pip2d(mesh2D,sampledProjectedOnXYplane)==1)
  }
  innerPoints2D<-sampledProjectedOnXYplane[innerpointsIndices,]
  
  if(plotting==TRUE){
    plot(rbind(mesh2D,mesh2D[1,]),type = "l",xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(innerPoints2D,pch='.',xlim = xlim,ylim = xlim,xlab = "",ylab = "") 
  }
  
  ######################################################################################################
  ######################################################################################################
  # divide 2D mesh into top and bottom parts
  
  edgeLengths<-rep(NA,(nrow(mesh2D)-1))
  for (i in 1:(nrow(mesh2D)-1)) {
    edgeLengths[i]<-norm(mesh2D[i,]-mesh2D[i+1,],type = '2')
  }
  
  total_Perimeter<-sum(edgeLengths)+norm(mesh2D[nrow(mesh2D),]-mesh2D[1,],type = '2')
  
  #Choose 100 points for the 2D mesh
  # acceptableEdges<-sort(which.maxn(edgeLengths,n = 100))
  # mesh2D<-mesh2D[acceptableEdges,]
  
  if(plotting==TRUE){
    plot(rbind(mesh2D,mesh2D[1,]),type = "l",xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(rbind(mesh2D,mesh2D[1,]),col='blue',xlim = xlim,ylim = xlim,xlab = "",ylab = "")
  }
  
  numberOfVertices<-nrow(mesh2D)
  geoDistancesMatrix<-array(NA,dim = c(numberOfVertices,numberOfVertices))
  for (i in 1:nrow(mesh2D)) {
    for (j in 1:nrow(mesh2D)) {
      geoDistancesMatrix[i,j]<-geodesic2D_of2Vertices(mesh2D = mesh2D,
                                                      vertexIndex1 = i,
                                                      vertexIndex2 = j,
                                                      edgeLengths = edgeLengths,
                                                      total_Perimeter = total_Perimeter)
    }
  }
  
  
  #make geoDistancesMatrix symmetric
  geoDistancesMatrix<-forceSymmetric(geoDistancesMatrix)
  # geoDistancesMatrix<-upper.tri(geoDistancesMatrix,diag = TRUE)+t(upper.tri(geoDistancesMatrix,diag = FALSE))
  
  
  # refVec<-c(1:nrow(mesh2D),1:nrow(mesh2D),1:nrow(mesh2D))
  # 
  # verticesCurvature<-rep(NA,nrow(mesh2D))
  # kRing4Curvature<-5
  # for (i in 1:nrow(mesh2D)) {
  #   verticesCurvature[i]<-curvatureOf3Points(point1 = mesh2D[refVec[nrow(mesh2D)+i-kRing4Curvature],],
  #                                            point2 = mesh2D[i,],
  #                                            point3 = mesh2D[refVec[nrow(mesh2D)+i+kRing4Curvature],],
  #                                            mesh2D=mesh2D)
  # }
  
  
  # if(plotting==TRUE){
  #   plot(rbind(mesh2D,mesh2D[1,]),type = "l",xlim = xlim,ylim = xlim,xlab = "",ylab = "")
  #   par(new=TRUE)
  #   plot(mesh2D[which(verticesCurvature>=0),],col='blue',pch=20,xlim = xlim,ylim = xlim,xlab = "",ylab = "")
  #   par(new=TRUE)
  #   plot(mesh2D[which(verticesCurvature<0),],col='red',pch=20,xlim = xlim,ylim = xlim,xlab = "",ylab = "")
  #   par(new=TRUE)
  #   maxCurveIndx<-which.max(verticesCurvature)
  #   minCurveIndx<-which.min(verticesCurvature)
  #   plot(rbind(mesh2D[maxCurveIndx,],mesh2D[maxCurveIndx,]),col='black',pch=16,xlim = xlim,ylim = xlim,xlab = "",ylab = "")
  #   par(new=TRUE)
  #   plot(rbind(mesh2D[minCurveIndx,],mesh2D[minCurveIndx,]),col='black',pch=16,xlim = xlim,ylim = xlim,xlab = "",ylab = "") 
  # }
  
  unitNormals<-verticesNormals_mesh2D(mesh2D = mesh2D)
  
  if(plotting==TRUE){
    plot(rbind(mesh2D,mesh2D[1,]),col='blue',xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    # par(new=TRUE)
    # plot(rbind(mesh2D,mesh2D[1,]),type = "l",xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    for (i in 1:nrow(mesh2D)) {
      plot(rbind(mesh2D[i,],mesh2D[i,]+unitNormals[i,]),type = "l",xlim = xlim,ylim = xlim,xlab = "",ylab = "")
      par(new=TRUE)
    } 
  }
  
  normalMultiplicationMatrix<-unitNormals%*%t(unitNormals)
  
  normalDistanceMatrix<-array(NA,dim = dim(geoDistancesMatrix))
  for (i in 1:nrow(mesh2D)) {
    for (j in 1:nrow(mesh2D)) {
      normalDistanceMatrix[i,j]<-geodesicDistance(unitNormals[i,],unitNormals[j,])
    }
  }
  
  # curvatureDistanceMatrix<-array(NA,dim = dim(geoDistancesMatrix))
  # for (i in 1:nrow(mesh2D)) {
  #   for (j in 1:nrow(mesh2D)) {
  #     # curvatureDistanceMatrix[i,j]<-abs(abs(verticesCurvature[i])-abs(verticesCurvature[j]))
  #     curvatureDistanceMatrix[i,j]<-abs((verticesCurvature[i])-verticesCurvature[j])
  #   }
  # }
  
  geoDistancesMatrix_pseudo<-array(NA,dim = dim(geoDistancesMatrix))
  total_Perimeter_pseudo<-nrow(mesh2D)
  for (i in 1:dim(geoDistancesMatrix_pseudo)[1]) {
    for (j in 1:dim(geoDistancesMatrix_pseudo)[2]) {
      if(abs(i-j)<total_Perimeter_pseudo/2){
        geoDistancesMatrix_pseudo[i,j]<-abs(i-j)
      }else{
        geoDistancesMatrix_pseudo[i,j]<-abs(total_Perimeter_pseudo-abs(i-j)) 
      }
    }
  }
  
  
  # elementWiseMultipleMatrix<-exp(lambda2D*normalMultiplicationMatrix*geoDistancesMatrix)
  # elementWiseMultipleMatrix<-exp(lambda2D*normalMultiplicationMatrix*geoDistancesMatrix_pseudo)
  elementWiseMultipleMatrix<-exp(lambda2D*cos(normalDistanceMatrix)*geoDistancesMatrix_pseudo)
  # elementWiseMultipleMatrix<-exp(lambda2D*normalMultiplicationMatrix)
  
  
  #compute affinity matrix W
  allAdjacentVertices<-allAdjacentVertices2D(mesh2D=mesh2D,k_Ring=k_Ring2D)
  
  W<-elementWiseMultipleMatrix
  tempVec<-1:nrow(mesh2D)
  for (i in 1:nrow(mesh2D)) {
    notNeighbor<-tempVec[!tempVec %in% allAdjacentVertices[[i]]]
    W[i,notNeighbor]<-0
  }
  diag(W)<-1
  
  d<-(rowSums(W))^(-1/2)
  D<-diag(d)
  
  L <- D%*%W%*%D
  
  L <- diag(nrow(L))-L
  
  eigenValAndVec<-eigen(L)
  
  # plot(1:nrow(mesh2D),eigenValAndVec$values[1:nrow(mesh2D)])
  
  smallestEigenVectorIndex<-nrow(L)-1
  # secondSmallestEigenValue<-eigenValAndVec$values[smallestEigenVectorIndex]
  secondSmallestEigenVector<-eigenValAndVec$vectors[,smallestEigenVectorIndex]
  
  
  # upIndices<-which(secondSmallestEigenVector>0)
  # downIndices<-which(secondSmallestEigenVector<=0)
  # if(plotting==TRUE){
  #   plot(rbind(mesh2D,mesh2D[1,]),type = "l",xlim = xlim,ylim = xlim,xlab = "",ylab = "")
  #   par(new=TRUE)
  #   plot(mesh2D[upIndices,],type = "l",col='blue',xlim = xlim,ylim = xlim,xlab = "",ylab = "")
  #   par(new=TRUE)
  #   plot(mesh2D[downIndices,],type = "l",col='red',xlim = xlim,ylim = xlim,xlab = "",ylab = "")
  # }
  
  edgeSuspects<-which.minn(abs(secondSmallestEigenVector),n = 4*k_Ring2D)
  
  # edgeSuspects<-which(round(abs(secondSmallestEigenVector),diits = 2)==0)
  
  if(plotting==TRUE){
    plot(rbind(mesh2D,mesh2D[1,]),type = "l",xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(mesh2D[edgeSuspects,],pch=16,xlim = xlim,ylim = xlim,xlab = "",ylab = "")
  }
  
  
  #classifying suspicious points into two groups
  dada4Classification<-mesh2D[edgeSuspects,]
  kmeanClustering<-kmeans(x = dada4Classification, centers = 2)
  cluster1_indices<-edgeSuspects[which(kmeanClustering$cluster==1)]
  cluster2_indices<-edgeSuspects[which(kmeanClustering$cluster==2)]
  if(plotting==TRUE){
    plot(rbind(mesh2D,mesh2D[1,]),type = "l",xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(mesh2D[c(cluster1_indices,cluster1_indices),],col='blue',pch=16,xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(mesh2D[c(cluster2_indices,cluster2_indices),],col='red',pch=16,xlim = xlim,ylim = xlim,xlab = "",ylab = "")
  }
  
  if(length(cluster1_indices)>1){
    geoDistancesMatrix_cluster1<-geoDistancesMatrix[cluster1_indices,cluster1_indices] 
    meanIndexCluster1<-cluster1_indices[which.min(rowSums(geoDistancesMatrix_cluster1))]
  }
  if(length(cluster2_indices)>1){
    geoDistancesMatrix_cluster2<-geoDistancesMatrix[cluster2_indices,cluster2_indices]
    meanIndexCluster2<-cluster2_indices[which.min(rowSums(geoDistancesMatrix_cluster2))]
  }
  if(length(cluster1_indices)==1){
    meanIndexCluster1<-cluster1_indices
  }
  if(length(cluster2_indices)==1){
    meanIndexCluster2<-cluster2_indices
  }
  
  if(plotting==TRUE){
    plot(rbind(mesh2D,mesh2D[1,]),type = "l",xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(mesh2D[c(meanIndexCluster1,meanIndexCluster1),],col='blue',pch=16,xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(mesh2D[c(meanIndexCluster2,meanIndexCluster2),],col='red',pch=16,xlim = xlim,ylim = xlim,xlab = "",ylab = "")
  }
  
  startUp<-meanIndexCluster1
  endUp<-meanIndexCluster2
  
  # distanceControlMinimumLimit<-total_Perimeter_pseudo/3
  # 
  # #infinity matrix based on the k-ring distance
  # # geoDistancesMatrix_Kring<-geoDistancesMatrix
  # geoDistancesMatrix_Kring<-geoDistancesMatrix_pseudo
  # #phase 1
  # geoDistancesMatrix_Kring[geoDistancesMatrix_Kring<distanceControlMinimumLimit]<-Inf
  # #phase 2
  # tempVec<-1:nrow(mesh2D)
  # for (i in 1:nrow(mesh2D)) {
  #   geoDistancesMatrix_Kring[i,allAdjacentVertices[[i]]]<-Inf
  # }
  # diag(geoDistancesMatrix_Kring)<-Inf
  # 
  # geoDistancesMatrix_KringCut<-geoDistancesMatrix_Kring[edgeSuspects,edgeSuspects]
  # suspectsWithMaxGeodesicDistance<-which(geoDistancesMatrix_KringCut == min(geoDistancesMatrix_KringCut), arr.ind = TRUE)[1,]
  # # geoDistancesMatrixCut<-geoDistancesMatrix[edgeSuspects,edgeSuspects]
  # # suspectsWithMaxGeodesicDistance<-which(geoDistancesMatrixCut == max(geoDistancesMatrixCut), arr.ind = TRUE)[1,]
  # # normalDistanceMatrixCut<-normalDistanceMatrix[edgeSuspects,edgeSuspects]
  # # suspectsWithMaxGeodesicDistance<-which(normalDistanceMatrixCut == min(normalDistanceMatrixCut), arr.ind = TRUE)[1,]
  # 
  # startUp<-edgeSuspects[suspectsWithMaxGeodesicDistance][1]
  # endUp<-edgeSuspects[suspectsWithMaxGeodesicDistance][2]
  
  
  if(TRUE){
    plot(rbind(mesh2D,mesh2D[1,]),type = "l",xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(mesh2D[c(startUp,endUp),],pch=16,xlim = xlim,ylim = xlim,xlab = "",ylab = "")
  }
  
  readline(prompt="Check the head and tail area and press [enter] to continue")
  
  #find up and down indices
  if(startUp>endUp){
    temp<-startUp
    startUp<-endUp
    endUp<-temp
  }
  
  upIndices<-startUp:endUp
  tempVec<-1:nrow(mesh2D)
  downIndices<-tempVec[!tempVec %in% upIndices]
  
  #fix sorting issue
  upIndices<-sortIndices(upIndices)
  downIndices<-sortIndices(downIndices)
  
  if(plotting==TRUE){
    plot(rbind(mesh2D,mesh2D[1,]),type = "l",xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(mesh2D[upIndices,],type = "l",col='blue',xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(mesh2D[downIndices,],type = "l",col='red',xlim = xlim,ylim = xlim,xlab = "",ylab = "")
  }
  
  ######################################################################################################
  ######################################################################################################
  # divide 2D mesh into top and bottom parts
  
  #we consider the top and bottom parts of the object as two meshes
  # allMeshes2D<-list(mesh2D[upIndices,],mesh2D[downIndices,])
  
  #convert 2D mesh to 3D mesh
  upMesh<-mesh2D[upIndices,]
  downMesh<-mesh2D[downIndices,]
  
  upMeshRev<-cbind(rev(upMesh[,1]),rev(upMesh[,2]))
  downMeshRev<-cbind(rev(downMesh[,1]),rev(downMesh[,2]))
  
  verteciesUP2<-cbind(rbind(upMesh,upMeshRev),rep(c(1,-1),nrow(upMesh)))
  verteciesDown2<-cbind(rbind(downMesh,downMeshRev),rep(c(1,-1),nrow(downMesh)))
  
  
  if(plotting==TRUE){
    open3d()
    plot3d(verteciesUP2,type="l",radius = 0.2,col = "blue",expand = 10,box=FALSE,add = TRUE)
    plot3d(verteciesDown2,type="l",radius = 0.2,col = "red",expand = 10,box=FALSE,add = TRUE)
  }
  
  # readline(prompt="Press [enter] to continue")
  
  
  facesUp<-c()
  k<-1
  for (i in 1:(nrow(verteciesUP2)-2)) {
    facesUp<-rbind(facesUp,c(k,k+1,k+2))
    k<-k+1
  }
  
  facesDown<-c()
  k<-1
  for (i in 1:(nrow(verteciesDown2)-2)) {
    facesDown<-rbind(facesDown,c(k,k+1,k+2))
    k<-k+1
  }
  
  vertsUp2 <- rbind(t(as.matrix(verteciesUP2)),1)
  trglsUp2 <- as.matrix(t(facesUp))
  tmeshUp2 <- tmesh3d(vertsUp2, trglsUp2)
  
  vertsDown2 <- rbind(t(as.matrix(verteciesDown2)),1)
  trglsDown2 <- as.matrix(t(facesDown))
  tmeshDown2 <- tmesh3d(vertsDown2, trglsDown2)
  
  if(plotting==TRUE){
    open3d()
    shade3d(tmeshUp2, col="blue")  #surface mech
    shade3d(tmeshDown2, col="red")  #surface mech 
  }
  
  # v1s<-array(NA,dim(innerPoints2D))
  # v2s<-array(NA,dim(innerPoints2D))
  # forceDirections<-array(NA,dim(innerPoints2D))
  # reduce the number of innerpoints
  # sizeTemp2<-round(nrow(innerPoints2D)/2)
  # innerPoints2D<-innerPoints2D[sample(x = 1:nrow(innerPoints2D),size = sizeTemp2,replace = FALSE),]
  force<-array(NA,nrow(innerPoints2D))
  pb <- txtProgressBar(style = 3)
  for (k in 1:nrow(innerPoints2D)) {
    setTxtProgressBar(pb,k/nrow(innerPoints2D))
    
    center<-innerPoints2D[k,]
    
    # we use asymmetric circles to avoid antipodal vectors
    tempCircle<-sphereGenerator_2D(center = center,r = urchinRadius,
                                   n = circleResolution,asymmetric = TRUE)
    
    circleIn3D<-cbind(tempCircle,rep(0,nrow(tempCircle)))
    centerIn3D<-c(center,0)
    circleMesh3D<-as.mesh3d(circleIn3D)
    
    normalsTemp<-circleIn3D-matrix(rep(centerIn3D,nrow(circleIn3D)),ncol = 3,byrow = TRUE)
    
    circleMesh3D$normals<-rbind(apply(normalsTemp,FUN = convertVec2unitVec,MARGIN = 1),
                                rep(1,nrow(circleIn3D)))
    
    # shade3d(circleMesh3D, col="black")  #surface mech
    
    intersectionsUp<-vcgRaySearch(circleMesh3D,mesh = tmeshUp2)
    intersectionsDown<-vcgRaySearch(circleMesh3D,mesh = tmeshDown2)
    
    if(sum(intersectionsUp$quality)<2 | sum(intersectionsDown$quality)<2 |
       min(c(intersectionsUp$distance,intersectionsDown$distance))<  thresholdDistance2D){
      
      force[k]<-10^4
      
      next
      
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
      #   plot(rbind(center,tipOfCuttedUpSpokes[i,1:2]),type = 'l',xlim = xlim,ylim = xlim,xlab = "",ylab = "")
      #   par(new=TRUE)
      # }
      # for (i in 1:nrow(tipOfCuttedDownSpokes)) {
      #   plot(rbind(center,tipOfCuttedDownSpokes[i,1:2]),type = 'l',xlim = xlim,ylim = xlim,xlab = "",ylab = "")
      #   par(new=TRUE)
      # }
      # plot(vert2points(tmeshUp2)[,1:2],type = 'l',xlim = xlim,ylim = xlim,xlab = "",ylab = "")
      # par(new=TRUE)
      # plot(vert2points(tmeshDown2)[,1:2],type = 'l',xlim = xlim,ylim = xlim,xlab = "",ylab = "")
      
      if(length(tipOfCuttedUpSpokes)==3 & length(tipOfCuttedDownSpokes)!=3){
        tipOfCuttedSpokesWithLabels<-rbind(c(tipOfCuttedUpSpokes[1:2],1),
                                           cbind(tipOfCuttedDownSpokes[,1:2],rep(2,nrow(tipOfCuttedDownSpokes))))
      }else if(length(tipOfCuttedUpSpokes)!=3 & length(tipOfCuttedDownSpokes)==3){
        tipOfCuttedSpokesWithLabels<-rbind(cbind(tipOfCuttedUpSpokes[,1:2],rep(1,nrow(tipOfCuttedUpSpokes))),
                                           c(tipOfCuttedDownSpokes[1:2],2))  
      }else if(length(tipOfCuttedUpSpokes)==3 & length(tipOfCuttedDownSpokes)==3){
        tipOfCuttedSpokesWithLabels<-rbind(c(tipOfCuttedUpSpokes[1:2],1),
                                           c(tipOfCuttedDownSpokes[1:2],2))  
      }else{
        tipOfCuttedSpokesWithLabels<-rbind(cbind(tipOfCuttedUpSpokes[,1:2],rep(1,nrow(tipOfCuttedUpSpokes))),
                                           cbind(tipOfCuttedDownSpokes[,1:2],rep(2,nrow(tipOfCuttedDownSpokes))))  
      }
      
      # plot(allMeshesTemp,pch=20,col='blue',xlim = xlim,ylim = ylim,xlab = "",ylab = "")
      
      tipOfCuttedSpokes<-tipOfCuttedSpokesWithLabels[,1:2]
      
      # plot(rbind(mesh2D,mesh2D[1,]),type = "l",xlim = xlim,ylim = xlim)
      # par(new=TRUE)
      # for (i in 1:dim(tipOfCuttedSpokes)[1]) {
      #   plot(rbind(center,tipOfCuttedSpokes[i,]),type = "l",xlim = xlim,ylim = xlim,xlab = "",ylab = "")
      #   par(new=TRUE)
      # }
      
      tempMatrix1<-matrix(rep(center,nrow(tipOfCuttedSpokes)),ncol = 2,byrow = TRUE)
      vectorsMatrix<-tipOfCuttedSpokes-tempMatrix1
      
      labels<-tipOfCuttedSpokesWithLabels[,3]
      
      forceTemp<-forceFunctionAngleBased_4MultiObject(vectorsMatrix = vectorsMatrix,
                                                      labels = labels,
                                                      type = 'one')
      
      force[k]<-forceTemp$force
      # v1s[k,]<-forceTemp$v1
      # v2s[k,]<-forceTemp$v2
      # forceDirections[k,]<-forceTemp$forceDirection
    }  
  }
  close(pb)
  
  forceVectorsMagnitudes<-force
  
  roundedForce<-round(forceVectorsMagnitudes,digits = 2)
  uniqueRounded<-sort(unique(roundedForce))
  
  # if (plotting==TRUE) {
  #   barplot(uniqueRounded,roundedForce) 
  # }
  
  # #plot levels
  # if(plotting==TRUE){
  #   for (i in length(uniqueRounded):1) {
  #     plot(mesh2D,type = "l",xlim = xlim,ylim = xlim,xlab = "",ylab = "")
  #     cumLevel_ith_Points<-innerPoints2D[roundedForce %in% uniqueRounded[1:i],]
  #     par(new=TRUE)
  #     plot(cumLevel_ith_Points,pch=20,col='blue',xlim = xlim,ylim = xlim,xlab = "",ylab = "")
  #   }    
  # }
  
  # x y z
  i<-length(uniqueRounded)-2
  
  cumLevel_ith_Points<-innerPoints2D[roundedForce %in% uniqueRounded[1:i],]
  
  #add two vertex points i.e., the start and the end of top or bottom part
  cumLevel_ith_Points<-rbind(mesh2D[upIndices[1],],
                             mesh2D[upIndices[1],],
                             mesh2D[upIndices[1],],
                             cumLevel_ith_Points,
                             mesh2D[upIndices[length(upIndices)],],
                             mesh2D[upIndices[length(upIndices)],],
                             mesh2D[upIndices[length(upIndices)],])
  
  #plot 2D
  if(plotting==TRUE){
    plot(rbind(mesh2D,mesh2D[1,]),type = "l",xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(cumLevel_ith_Points,pch=20,col="blue",xlim = xlim,ylim = xlim,xlab = "",ylab = "") 
  }
  
  centroid2D<-colMeans(cumLevel_ith_Points)
  centeredPoints2D<-scale(cumLevel_ith_Points,center = TRUE,scale = FALSE)
  
  boundaryPoints<-mesh2D
  
  
  x4plot<-seq(from = xlim[1],to = xlim[2],length.out = 1000) 
  
  sumOfAllResiduals<-c()
  for (k in seq(0,360,rotationGap)) {
    tempTheta<-k*2*pi/360
    # cat("Rotation angle theta:",tempTheta,"\n")
    planeRotationMatrix<-matrix(c(cos(tempTheta),-sin(tempTheta),sin(tempTheta),cos(tempTheta)),ncol = 2,byrow = T)
    centeredPoints2DTemp<-centeredPoints2D%*%planeRotationMatrix
    
    # par(new=T)
    # plot(centeredPoints2DTemp,col="green",xlim = xlim,ylim = xlim, xlab = "",ylab = "")
    
    x<-centeredPoints2DTemp[,1]
    y<-centeredPoints2DTemp[,2]
    fit <- lm(y ~ poly(x, degree = polyDegree2D,raw = TRUE),
              data=as.data.frame(cbind(x,y))) #NB!!! raw=F calculate orthogonal polynomial
    
    
    if(plotting==TRUE){
      boundaryPointsTemp<-(boundaryPoints-matrix(rep(centroid2D,dim(boundaryPoints)[1]),
                                                 ncol = 2,byrow = TRUE))%*%planeRotationMatrix
      plot(rbind(boundaryPointsTemp,boundaryPointsTemp[1,]),type="l",xlim = xlim,ylim = xlim,xlab = "",ylab = "") 
      par(new=TRUE) 
      plot(centeredPoints2DTemp,pch=20,col="blue",xlim = xlim,ylim = xlim,xlab = "",ylab = "") 
      
      xDataTemp<-data.frame(x=x4plot)
      y4Plot<-predict(fit, newdata = xDataTemp)
      
      par(new=TRUE)
      plot(cbind(x4plot,y4Plot),type="l",lwd=2,col="green",xlim = xlim,ylim = xlim,xlab = "",ylab = "") 
    }
    
    sumOfAllResiduals<-c(sumOfAllResiduals,sum(fit$residuals^2))
    
  }
  
  k<-which.min(sumOfAllResiduals)
  tempTheta<-seq(0,360,rotationGap)[k]*2*pi/360   
  rotationMatrix2D<-matrix(c(cos(tempTheta),-sin(tempTheta),sin(tempTheta),cos(tempTheta)),ncol = 2,byrow = T)
  selectedMiddlePoint2D<-centeredPoints2D%*%rotationMatrix2D
  boundaryPointsTemp<-(boundaryPoints-matrix(rep(centroid2D,dim(boundaryPoints)[1]),
                                             ncol = 2,byrow = TRUE))%*%rotationMatrix2D
  
  if(plotting==TRUE){
    plot(rbind(boundaryPointsTemp,boundaryPointsTemp[1,]),type="l",xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(selectedMiddlePoint2D,pch=20,col="blue",xlim = xlim,ylim = xlim,xlab = "",ylab = "")
  }
  
  # x y z
  x<-selectedMiddlePoint2D[,1]
  y<-selectedMiddlePoint2D[,2]
  fit4_2D <- lm(y ~ poly(x, degree = polyDegree2D,raw = TRUE),
                data=as.data.frame(cbind(x,y))) #NB!!! raw=F calculate orthogonal polynomial
  
  
  x4plot<-seq(from = min(boundaryPointsTemp[,1]),to = max(boundaryPointsTemp[,1]),length.out = 1000)
  xDataTemp<-data.frame(x=x4plot)
  y4plot<-predict(fit4_2D, newdata = xDataTemp)
  
  curvePoints<-cbind(x4plot,y4plot)
  curvePointsInside<-curvePoints[pip2d(boundaryPointsTemp,curvePoints)==1,]
  
  if(plotting==TRUE){
    plot(rbind(boundaryPointsTemp,boundaryPointsTemp[1,]),type="l",xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(curvePointsInside,pch=20,col="blue",xlim = xlim,ylim = xlim,xlab = "",ylab = "") 
  }
  
  medialPointsEqualyDis_2D<-splineFitWithEquidistancePoints(samplePoints = curvePointsInside,
                                                            polyDegree2D = polyDegree2D,
                                                            numberOfEquidistancePoints = numberOfSpanialPoints)
  
  if(plotting==TRUE){
    plot(rbind(boundaryPointsTemp,boundaryPointsTemp[1,]),type="l",xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(medialPointsEqualyDis_2D,pch=20,col="blue",xlim = xlim,ylim = xlim,xlab = "",ylab = "")
  }
  #rotate boundary
  boundaryPoints<-(boundaryPoints-matrix(rep(centroid2D,dim(boundaryPoints)[1]),
                                         ncol = 2,byrow = TRUE))%*%rotationMatrix2D
  
  # plot(rbind(boundaryPoints,boundaryPoints[1,]),type="l",xlim = xlim,ylim = xlim,xlab = "",ylab = "")
  # par(new=TRUE)
  # plot(medialPointsEqualyDis_2D,pch=20,col="blue",xlim = xlim, ylim = xlim,xlab = "",ylab = "")
  
  #transfer back
  boundaryAndMedial<-rbind(boundaryPoints,medialPointsEqualyDis_2D)
  boundaryAndMedial_TransferBack<-boundaryAndMedial%*%t(rotationMatrix2D)+matrix(rep(centroid2D,dim(boundaryAndMedial)[1]),ncol = 2,byrow = TRUE)
  boundaryPoints<-boundaryAndMedial_TransferBack[1:nrow(boundaryPoints),]
  medialPointsEqualyDis_2D<-boundaryAndMedial_TransferBack[(nrow(boundaryPoints)+1):nrow(boundaryAndMedial),]
  
  #sort medial points from left to right
  if(medialPointsEqualyDis_2D[1,1]>medialPointsEqualyDis_2D[nrow(medialPointsEqualyDis_2D),1]){
    medialPointsEqualyDis_2D<-cbind(rev(medialPointsEqualyDis_2D[,1]),rev(medialPointsEqualyDis_2D[,2]))
  }
  
  if(plotting==TRUE){
    
    plot(mesh2D,type="l",xlim = xlim,ylim = xlim,xlab = "",ylab = "") 
    par(new=TRUE) 
    plot(medialPointsEqualyDis_2D,pch=20,col="blue",xlim = xlim,ylim = xlim,xlab = "",ylab = "")  
  }
  
  #local frames 2D on conditional medial locus
  
  skeletalPointsMethod_2D<-medialPointsEqualyDis_2D
  
  normalVectorsMedialLocus<-c()
  for (i in 2:(nrow(skeletalPointsMethod_2D)-1)) {
    p1<-skeletalPointsMethod_2D[i-1,]
    p2<-skeletalPointsMethod_2D[i,]
    p3<-skeletalPointsMethod_2D[i+1,]
    normalVectorsMedialLocus<-rbind(normalVectorsMedialLocus,normalOfaVertex(point1 = p1,vertex = p2,point2 = p3))
  }
  normalVectorsMedialLocus<-rbind(normalVectorsMedialLocus[1,],
                                  normalVectorsMedialLocus,
                                  normalVectorsMedialLocus[nrow(normalVectorsMedialLocus),])
  R1<-rotMat(c(0,1),c(1,0))
  tangentVectorsMedialLocus<-normalVectorsMedialLocus%*%t(R1)
  scaleFactor<-1
  normalVectorsMedialLocusScaled<-normalVectorsMedialLocus*scaleFactor
  tangentVectorsMedialLocusScaled<-tangentVectorsMedialLocus*scaleFactor
  
  
  #plot
  if(plotting==TRUE){
    plot(rbind(boundaryPoints,boundaryPoints[1,]),type = "l",xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    for (i in 2:(nrow(skeletalPointsMethod_2D)-1)) {
      plot(rbind(skeletalPointsMethod_2D[i,],skeletalPointsMethod_2D[i,]+normalVectorsMedialLocusScaled[i,]),col='black',lwd=2,type = 'l',xlim = xlim,ylim = xlim,xlab = "",ylab = "")
      par(new=TRUE)
      plot(rbind(skeletalPointsMethod_2D[i,],skeletalPointsMethod_2D[i,]+tangentVectorsMedialLocusScaled[i,]),col='black',lwd=2,type = 'l',xlim = xlim,ylim = xlim,xlab = "",ylab = "")
      par(new=TRUE)
    }
    plot(skeletalPointsMethod_2D,col='orange',type = 'l',xlim = xlim,ylim = xlim,xlab = "",ylab = "") 
  }
  
  spinalPoints<-skeletalPointsMethod_2D[2:(nrow(skeletalPointsMethod_2D)-1),]
  tipOfExtraSpokeHead<-skeletalPointsMethod_2D[1,]
  tipOfExtraSpokeTail<-skeletalPointsMethod_2D[nrow(skeletalPointsMethod_2D),]
  
  
  for (i in 2:(nrow(skeletalPointsMethod_2D)-1)) {
    tipOfTopSpokes<-cutAndStretchSpokes2D(allSpokesTips = spinalPoints+normalVectorsMedialLocus[2:(nrow(normalVectorsMedialLocus)-1),],
                                          allSpokesTails = spinalPoints,
                                          boundaryPoints1 = boundaryPoints1,
                                          boundaryPoints2 = boundaryPoints2)
    tipOfBottomSpokes<-cutAndStretchSpokes2D(allSpokesTips = spinalPoints-normalVectorsMedialLocus[2:(nrow(normalVectorsMedialLocus)-1),],
                                             allSpokesTails = spinalPoints,
                                             boundaryPoints1 = boundaryPoints1,
                                             boundaryPoints2 = boundaryPoints2)
  }
  
  
  if(plotting==TRUE){
    plot(rbind(boundaryPoints,boundaryPoints[1,]),type = "l", lty=3,xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    for (i in 1:nrow(spinalPoints)) {
      plot(rbind(tipOfTopSpokes[i,],spinalPoints[i,]),col="blue",type = 'l',xlim = xlim,ylim = xlim,xlab = "",ylab = "")
      par(new=TRUE)
    }
    plot(spinalPoints,col="black",pch=20,xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    for (i in 1:nrow(spinalPoints)) {
      plot(rbind(tipOfBottomSpokes[i,],spinalPoints[i,]),col="red",type = 'l',xlim = xlim,ylim = xlim,xlab = "",ylab = "")
      par(new=TRUE)
    }
    par(new=TRUE)
    plot(rbind(tipOfExtraSpokeHead,spinalPoints[1,]),col="black",type = 'l',xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(rbind(tipOfExtraSpokeTail,spinalPoints[nrow(spinalPoints),]),col="black",type = 'l',xlim = xlim,ylim = xlim,xlab = "",ylab = "")
  }
  
  
  dataTemp<-data.frame(x=tipOfTopSpokes[,1], y=tipOfTopSpokes[,2])
  tipOfTopSpokes_Z<-predict(fit4, newdata = dataTemp)
  tipOfTopSpokes3D<-cbind(tipOfTopSpokes,tipOfTopSpokes_Z)
  dataTemp<-data.frame(x=tipOfBottomSpokes[,1], y=tipOfBottomSpokes[,2])
  tipOfBottomSpokes_Z<-predict(fit4, newdata = dataTemp)
  tipOfBottomSpokes3D<-cbind(tipOfBottomSpokes,tipOfBottomSpokes_Z)
  dataTemp<-data.frame(x=spinalPoints[,1], y=spinalPoints[,2])
  spinalPoints_Z<-predict(fit4, newdata = dataTemp)
  spinalPoints3D<-cbind(spinalPoints,spinalPoints_Z)
  dataTemp<-data.frame(x=skeletalPointsMethod_2D[,1], y=skeletalPointsMethod_2D[,2])
  spinalPointsComplete_Z<-predict(fit4, newdata = dataTemp)
  spinalPoints3DComplete<-cbind(skeletalPointsMethod_2D,spinalPointsComplete_Z)
  tipOfExtraSpokes<-rbind(tipOfExtraSpokeHead,tipOfExtraSpokeTail)
  dataTemp<-data.frame(x=tipOfExtraSpokes[,1], y=tipOfExtraSpokes[,2])
  tipOfExtraSpokes_Z<-predict(fit4, newdata = dataTemp)
  tipOfExtraSpokes3D<-cbind(tipOfExtraSpokes,tipOfExtraSpokes_Z)
  tipOfExtraSpokeHead3D<-tipOfExtraSpokes3D[1,]
  tipOfExtraSpokeTail3D<-tipOfExtraSpokes3D[2,]
  
  #plot
  if(plotting==TRUE){
    open3d()
    for (i in 1:nrow(spinalPoints3D)) {
      plot3d(rbind(spinalPoints3D[i,],tipOfTopSpokes3D[i,]),type="l",col = "blue",expand = 10,box=FALSE,add = TRUE)
      plot3d(rbind(spinalPoints3D[i,],tipOfBottomSpokes3D[i,]),type="l",col = "red",expand = 10,box=FALSE,add = TRUE)
    }
    plot3d(rbind(spinalPoints3D[1,],tipOfExtraSpokeHead3D),type="l",col = "black",expand = 10,box=FALSE,add = TRUE)
    plot3d(rbind(spinalPoints3D[nrow(spinalPoints3D),],tipOfExtraSpokeTail3D),type="l",col = "black",expand = 10,box=FALSE,add = TRUE)
    #plot points
    plot3d(spinalPoints3D,type = 's',radius = 0.2,col = "black",expand = 10,box=FALSE,add = TRUE)
    plot3d(tipOfTopSpokes3D,type = 's',radius = 0.2,col = "blue",expand = 10,box=FALSE,add = TRUE)
    plot3d(tipOfBottomSpokes3D,type = 's',radius = 0.2,col = "red",expand = 10,box=FALSE,add = TRUE)
    plot3d(rbind(tipOfExtraSpokeHead3D,tipOfExtraSpokeTail3D),type = 's',radius = 0.2,col = "black",expand = 10,box=FALSE,add = TRUE)
    #plot medial surface
    # plot3d(medialPoints2,type="p",col = "red",expand = 10,box=FALSE,add = TRUE)
    shade3d(tmesh, col="blue",alpha=0.2) 
  }
  
  ######################################################################################################
  ######################################################################################################
  # cross-sections!!!!
  ######################################################################################################
  ######################################################################################################
  # cross-sections!!!!
  
  #plot cross-sections
  if(plotting==TRUE & cross_sections_visualization==TRUE){
    
    normalOfCrossSections<-c()
    for (i in 1:nrow(spinalPoints3D)) {
      normalOfCrossSections<-rbind(normalOfCrossSections,
                                   convertVec2unitVec2(myCrossProduct(tipOfTopSpokes3D[i,]-spinalPoints3D[i,],
                                                                      tipOfBottomSpokes3D[i,]-spinalPoints3D[i,])))
    }
    
    for (i in 1:nrow(spinalPoints3D)) {
      # plot3d(rbind(spinalPoints3D[i,],spinalPoints3D[i,]+normalOfCrossSections[i,]),type="l",col = "black",expand = 10,box=FALSE,add = TRUE)
      vectors3d(spinalPoints3D[i,]+normalOfCrossSections[i,],origin = spinalPoints3D[i,],headlength = 0.2,radius = 1/6, col="blue", lwd=1)
    }
    
    shade3d(tmesh, col="white",alpha=0.2)  
    
    plot3d(spinalPoints3D,type="l",lwd = 4,col = "black",expand = 10,box=FALSE,add = TRUE)
    plot3d(spinalPoints3D,type = 's',radius = 0.2,col = "black",expand = 10,box=FALSE,add = TRUE)
    plot3d(rbind(spinalPoints3D[1,],tipOfExtraSpokeHead3D),type="l",lwd = 2,col = "black",expand = 10,box=FALSE,add = TRUE)
    plot3d(rbind(spinalPoints3D[nrow(spinalPoints3D),],tipOfExtraSpokeTail3D),type="l",lwd = 2,col = "black",expand = 10,box=FALSE,add = TRUE)
    
    plot3d(tipOfTopSpokes3D,type = 's',radius = 0.2,col = "blue",expand = 10,box=FALSE,add = TRUE)
    plot3d(tipOfBottomSpokes3D,type = 's',radius = 0.2,col = "red",expand = 10,box=FALSE,add = TRUE)
    plot3d(rbind(tipOfExtraSpokeHead3D,tipOfExtraSpokeHead3D),type = 's',radius = 0.2,col = "black",expand = 10,box=FALSE,add = TRUE)
    plot3d(rbind(tipOfExtraSpokeTail3D,tipOfExtraSpokeTail3D),type = 's',radius = 0.2,col = "black",expand = 10,box=FALSE,add = TRUE)
    
    Sigma1 <- matrix(c(20, 0, 0, 0, 20, 0, 0, 0, 0.01), 3, 3)
    elipseTemp1<-ellipse3d(Sigma1)
    
    slicingEllipsoids<-list()
    for (i in 1:nrow(spinalPoints3D)) {
      # elipseTemp2<-rotate3d(elipseTemp1,
      #                       matrix = rbind(cbind(rotMat(c(0,0,1),normalOfCrossSections[i,]),c(0,0,0)),c(0,0,0,1)))
      # 
      # elipseTemp4<-translate3d(elipseTemp2,
      #                          x = spinalPoints3D[i,1],
      #                          y = spinalPoints3D[i,2],
      #                          z = spinalPoints3D[i,3])
      # slicingEllipsoids[[i]]<-elipseTemp4
      
      R1<-rotMat(c(0,0,1),normalOfCrossSections[i,])
      elipseTemp2<-elipseTemp1
      elipseTemp2$vb<-t(cbind(vert2points(elipseTemp2)%*%t(R1),rep(1,nrow(vert2points(elipseTemp2)))))
      
      elipseTemp4<-translate3d(elipseTemp2,
                               x = spinalPoints3D[i,1],
                               y = spinalPoints3D[i,2],
                               z = spinalPoints3D[i,3])
      slicingEllipsoids[[i]]<-elipseTemp4
    }
    
    for (i in 1:nrow(spinalPoints3D)) {
      plot3d(slicingEllipsoids[[i]], col = "grey", alpha = 0.2, add = TRUE)
    }
    shade3d(tmesh, col="white",alpha=0.2)
    
    # convert slicingEllipsoids to triangle meshes
    triangleMeshesOfSlicingEllipsoids<-list()
    for (i in 1:nrow(spinalPoints3D)) {
      triangleMeshesOfSlicingEllipsoids[[i]]<-vcgUpdateNormals(as.tmesh3d(slicingEllipsoids[[i]]))
    }
    
    #increase the resolution of the cross-sections and the medial surface
    numberOfCoverPoints3D_2<-5*(10^6)
    xs<-runif(n = numberOfCoverPoints3D_2,min =minTempX ,max =maxTempX)
    ys<-runif(n = numberOfCoverPoints3D_2,min =minTempY ,max =maxTempY)
    zs<-runif(n = numberOfCoverPoints3D_2,min =minTempZ ,max =maxTempZ)
    
    coverpoint<-as.matrix(cbind(xs,ys,zs)) #as.matrix is necessary
    
    
    #find internal points
    insidePoints<-pip3d(Vertices = t(tmesh$vb[1:3,]),
                        Faces = t(tmesh$it),
                        Queries = coverpoint)
    innerPoints<-coverpoint[insidePoints==1,]
    
    x2<-coverpoint[,1]
    y2<-coverpoint[,2]
    newDATA<-data.frame(x=x2, y=y2)
    surfacePointsZ<-predict(fit4, newdata = newDATA)
    
    medialPoints<-cbind(x2,y2,surfacePointsZ)
    
    insidePoints2<-pip3d(Vertices = t(tmesh$vb[1:3,]),
                         Faces = t(tmesh$it),
                         Queries = medialPoints)
    medialPoints2<-medialPoints[insidePoints2==1,]
    
    
    #cross-sections
    plot3d(innerPoints,type = 'p',col = "gray",expand = 10,box=FALSE,add = TRUE)
    open3d()
    plot3d(medialPoints2,type = 'p',col = "red",expand = 10,box=FALSE,add = TRUE)
    shade3d(tmesh,col='white',alpha=0.2)
    
    cross_sections_points<-list()
    cross_sections_Medial_points<-list()
    for (i in 1:nrow(spinalPoints3D)) {
      print(i/nrow(spinalPoints3D))
      choose_cross_sections_points<-pip3d(Vertices = t(slicingEllipsoids[[i]]$vb[1:3,]),
                                          Faces = t(triangleMeshesOfSlicingEllipsoids[[i]]$it),
                                          Queries = innerPoints)
      choose_cross_sections_Medial_points<-pip3d(Vertices = t(slicingEllipsoids[[i]]$vb[1:3,]),
                                                 Faces = t(triangleMeshesOfSlicingEllipsoids[[i]]$it),
                                                 Queries = medialPoints2)
      
      cross_sections_points[[i]]<-innerPoints[choose_cross_sections_points==1,]
      cross_sections_Medial_points[[i]]<-medialPoints2[choose_cross_sections_Medial_points==1,] 
    }
    
    cross_sections_points_2D<-list()
    cross_sections_Medial_points_2D<-list()
    for (i in 1:nrow(spinalPoints3D)) {
      pcaTemp<-prcomp(rbind(cross_sections_points[[i]],cross_sections_Medial_points[[i]]))
      cross_sections_points_2D[[i]]<-pcaTemp$x[1:nrow(cross_sections_points[[i]]),]
      cross_sections_Medial_points_2D[[i]]<-pcaTemp$x[-c(1:nrow(cross_sections_points[[i]])),]
    }
    
    xlim2<-c(-15,15)
    for (i in 1:nrow(spinalPoints3D)) {
      print(i/nrow(spinalPoints3D))
      plot(cross_sections_points_2D[[i]],col='black',pch='.',xlim = xlim2,ylim = xlim2,xlab = "",ylab = "")
      par(new=TRUE)
      plot(cross_sections_Medial_points_2D[[i]],col='red',pch='.',xlim = xlim2,ylim = xlim2,xlab = "",ylab = "")
    }
    
  }
  
  
  
  
  ###################################################################################################################
  ###################################################################################################################
  # 3D skeletal
  
  #store and remove first and the last point
  # tipOfTopSpokes<-chordTipUpSelected[2:nrow(chordTipUpSelected),]
  # tipOfBottomSpokes<-chordTipDownSelected[2:nrow(chordTipDownSelected),]
  if(skeletalPointsMethod_2D[1,1]<skeletalPointsMethod_2D[nrow(skeletalPointsMethod_2D),1]){
    tipOfExtraSpokeHead<-skeletalPointsMethod_2D[1,]
    tipOfExtraSpokeTail<-skeletalPointsMethod_2D[nrow(skeletalPointsMethod_2D),]
  }else{
    tipOfExtraSpokeHead<-skeletalPointsMethod_2D[nrow(skeletalPointsMethod_2D),]
    tipOfExtraSpokeTail<-skeletalPointsMethod_2D[1,]
  }
  spinalPoints<-skeletalPointsMethod_2D[2:(nrow(skeletalPointsMethod_2D)-1),]
  
  # #store and remove first and the last point
  # tipOfTopSpokes<-chordTipUpSelected[2:nrow(chordTipUpSelected),]
  # tipOfBottomSpokes<-chordTipDownSelected[2:nrow(chordTipDownSelected),]
  # tipOfExtraSpokeHead<-middileChordsSelected[1,]
  # tipOfExtraSpokeTail<-middileChordsSelected[nrow(middileChordsSelected),]
  # spinalPoints<-middileChordsSelected[2:(nrow(middileChordsSelected)-1),]
  
  if(plotting==TRUE){
    plot(rbind(boundaryPoints,boundaryPoints[1,]),type = "l", lty=3,xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    for (i in 1:nrow(spinalPoints)) {
      plot(rbind(tipOfTopSpokes[i,],spinalPoints[i,]),col="blue",type = 'l',xlim = xlim,ylim = xlim,xlab = "",ylab = "")
      par(new=TRUE)
    }
    plot(spinalPoints,col="black",pch=20,xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    for (i in 1:nrow(spinalPoints)) {
      plot(rbind(tipOfBottomSpokes[i,],spinalPoints[i,]),col="red",type = 'l',xlim = xlim,ylim = xlim,xlab = "",ylab = "")
      par(new=TRUE)
    }
    par(new=TRUE)
    plot(rbind(tipOfExtraSpokeHead,spinalPoints[1,]),col="black",type = 'l',xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(rbind(tipOfExtraSpokeTail,spinalPoints[nrow(spinalPoints),]),col="black",type = 'l',xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    
  }
  
  
  # spokes interpolation
  topSpokesPoints<-c()
  for (i in 1:nrow(spinalPoints)) {
    tempPoints<-generatePointsBetween2Points(spinalPoints[i,],
                                             tipOfTopSpokes[i,],
                                             numberOfPoints = numberOf2DspokePoints) 
    topSpokesPoints<-rbind(topSpokesPoints,tempPoints)
  }
  bottomSpokesPoints<-c()
  for (i in 1:nrow(spinalPoints)) {
    tempPoints<-generatePointsBetween2Points(spinalPoints[i,],
                                             tipOfBottomSpokes[i,],
                                             numberOfPoints = numberOf2DspokePoints) 
    bottomSpokesPoints<-rbind(bottomSpokesPoints,tempPoints)
  }
  extraSpokeHeadPoints<-tipOfExtraSpokeHead
  extraSpokeTailPoints<-tipOfExtraSpokeTail
  
  #plot 2D
  if(plotting==TRUE){
    plot(rbind(boundaryPoints,boundaryPoints[1,]),type = "l", lty=3,xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(topSpokesPoints,pch=20,cex=0.2,col="blue",xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(bottomSpokesPoints,pch=20,cex=0.2,col="red",xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(rbind(extraSpokeHeadPoints,extraSpokeTailPoints),pch=20,cex=0.2,col="black",xlim = xlim,ylim = xlim,xlab = "",ylab = "")
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
  
  #add two extra head and tail points
  skeletalPoints_2D<-rbind(extraSpokeHeadPoints,
                           skeletalPoints_2D,
                           extraSpokeTailPoints)
  
  #plot 2D
  if(plotting==TRUE){
    plot(skeletalPoints_2D,type = "p",pch=20,col="blue",xlim = xlim,ylim = xlim,xlab = "",ylab = "")   
  }
  
  #transfer back
  # skeletalPoints_2D<-skeletalPoints_2D%*%t(rotationMatrix2D)+matrix(rep(centroid2D,dim(skeletalPoints_2D)[1]),ncol = 2,byrow = TRUE)
  
  # #plot 2D
  # if(plotting==TRUE){
  #   plot(skeletalPoints_2D,type = "l",pch=20,col="blue",xlim = xlim,ylim = xlim,xlab = "",ylab = "") 
  # }
  
  
  newData_2d<-data.frame(x=skeletalPoints_2D[,1], y=skeletalPoints_2D[,2])
  surfacePointsZ2<-predict(fit4, newdata = newData_2d)
  medialPoints3D<-cbind(skeletalPoints_2D[,1],skeletalPoints_2D[,2],surfacePointsZ2)
  
  #plot
  if(plotting==TRUE){
    open3d()
    plot3d(medialPoints3D,type="s",radius = 0.2,col = "black",expand = 10,box=FALSE,add = TRUE)
    # wire3d(tmesh, col="blue")  #wire mesh
    shade3d(tmesh, col="white",alpha=0.2)
  }
  
  ######################################################################################################
  ######################################################################################################
  # local frames 3D
  
  #rearrange medial points
  medialPoints3D<-medialPoints3D[c(2:(nrow(medialPoints3D)-1),1,nrow(medialPoints3D)),]
  
  numberOfMedialPoints3D<-nrow(medialPoints3D)
  
  numberOfLayers<-2*numberOf2DspokePoints-1
  
  groupIndicesOfVeins<-list()
  k<-1
  i<-1
  while (k<(numberOfMedialPoints3D-2)) {
    groupIndicesOfVeins[[i]]<-k:(k+(numberOf2DspokePoints+2))
    k<-k+(2*numberOf2DspokePoints-1)
    i<-i+1
  }
  
  middleCrossSectionIndex<-ceiling(length(groupIndicesOfVeins)/2)
  skeletal_CentroidIndex<-groupIndicesOfVeins[[middleCrossSectionIndex]][1]
  
  pointsIndices<-1:numberOfMedialPoints3D
  parentsIndices<-rep(NA,numberOfMedialPoints3D)
  for (i in 1:(length(pointsIndices)-2)) {
    if(pointsIndices[i]==skeletal_CentroidIndex){
      parentsIndices[i]<-i
    }else if(pointsIndices[i]%%numberOfLayers==1){ #spinal point
      if(pointsIndices[i]<=skeletal_CentroidIndex){
        parentsIndices[i]<-i+numberOfLayers
      }else{
        parentsIndices[i]<-i-numberOfLayers
      }
    }else if(pointsIndices[i]%%numberOfLayers==(numberOf2DspokePoints+1)){
      parentsIndices[i]<-i-numberOf2DspokePoints
    }else{
      parentsIndices[i]<-i-1
    }
  }
  parentsIndices[numberOfMedialPoints3D-1]<-1
  parentsIndices[numberOfMedialPoints3D]<-groupIndicesOfVeins[[length(groupIndicesOfVeins)]][1]
  
  #plot connections
  if(plotting==TRUE){
    for (i in 1:length(pointsIndices)) {
      k1<-pointsIndices[i]
      k2<-parentsIndices[i]
      if(k2==k1){
        next
      }
      vectors3d(medialPoints3D[k1,],
                origin = medialPoints3D[k2,],
                # labels = c("",k1),
                headlength = 0.2,radius = 1/6, col="blue", lwd=1)
    }
  }
  
  
  backPointsIndices<-parentsIndices
  backPointsIndices[skeletal_CentroidIndex]<-groupIndicesOfVeins[[middleCrossSectionIndex+1]][1]
  
  pointsIndices<-1:numberOfMedialPoints3D
  frontPointsIndices<-rep(NA,numberOfMedialPoints3D)
  for (i in 1:(length(pointsIndices)-2)) {
    if(pointsIndices[i]==skeletal_CentroidIndex){
      frontPointsIndices[i]<-groupIndicesOfVeins[[middleCrossSectionIndex-1]][1]
    }else if(pointsIndices[i]%%numberOfLayers==1){ #spinal point
      if(pointsIndices[i]==1){
        frontPointsIndices<-length(pointsIndices)-1
      }else if(pointsIndices[i]==(length(pointsIndices)-numberOfLayers-1)){
        frontPointsIndices[i]<-length(pointsIndices)
      }else if(pointsIndices[i]<=skeletal_CentroidIndex){
        frontPointsIndices[i]<-i-numberOfLayers
      }else{
        frontPointsIndices[i]<-i+numberOfLayers
      }
    }else if(pointsIndices[i]%%numberOfLayers==numberOf2DspokePoints | 
             pointsIndices[i]%%numberOfLayers==0){
      next
    }else{
      frontPointsIndices[i]<-i+1
    }
  }
  frontPointsIndices<-c(frontPointsIndices,NA,NA)
  
  #plot LP_ds_rep
  if(plotting==TRUE){
    for (i in 1:length(pointsIndices)) {
      k1<-pointsIndices[i]
      k2<-backPointsIndices[i]
      if(k2==k1){
        next
      }
      vectors3d(medialPoints3D[k1,],
                origin = medialPoints3D[k2,],
                labels = c("",k1),headlength = 0.2,radius = 1/6, col="red", lwd=1)
    }
    for (i in 1:length(pointsIndices)) {
      k1<-frontPointsIndices[i]
      k2<-pointsIndices[i]
      if(k2==k1 | is.na(k1)){
        next
      }
      vectors3d(medialPoints3D[k1,],
                origin = medialPoints3D[k2,],
                labels = c("",k1),headlength = 0.25,radius = 2/6, col="blue", lwd=1)
    }
  }
  
  
  spineIndices<-c()
  for (i in 1:length(groupIndicesOfVeins)) {
    spineIndices<-c(spineIndices,groupIndicesOfVeins[[i]][1])
  }
  spineIndices<-c(numberOfMedialPoints3D-1,spineIndices,numberOfMedialPoints3D)
  
  framesCenters<-rep(NA,length(pointsIndices))
  for (i in 1:(length(pointsIndices)-2)) {
    if(pointsIndices[i]%%numberOfLayers == 0 |
       pointsIndices[i]%%numberOfLayers == numberOf2DspokePoints){
      next
    }
    framesCenters[i]<-pointsIndices[i]
  }
  
  
  # normal vectors
  # temp<-normalsOfSkeletalSheet(centeredSkel = medialPoints3D)
  # medialNormals<-temp$medialNormals
  
  frameIndices<-which(!is.na(framesCenters))
  
  backPointsIndices4spineNormals<-backPointsIndices
  frontPointsIndices4spineNormals<-frontPointsIndices
  for (i in frameIndices) {
    if(i %in% spineIndices[-c(1,length(spineIndices))]){
      backPointsIndices4spineNormals[i]<-i+1
      frontPointsIndices4spineNormals[i]<-i+numberOf2DspokePoints
    }
  }
  
  medialNormals<-normalsOfSkeletalSheetByConnections(skelPoints = medialPoints3D,
                                                     framesCenters = framesCenters,
                                                     backPointsIndices = backPointsIndices4spineNormals,
                                                     frontPointsIndices = frontPointsIndices4spineNormals)
  

  
  # vectors3d(medialPoints3D[frameIndices,]+medialNormals[frameIndices,],origin = medialPoints3D[frameIndices,],headlength = 0.1,radius = 1/10, col="orange", lwd=2)
  
  
  framesElements<-frameGenerator7(skelPoints = medialPoints3D,
                                  medialNormals=medialNormals,
                                  framesCenters=framesCenters,
                                  framesBackPoints=backPointsIndices,
                                  framesFrontPoints=frontPointsIndices)
  
  framesFirstVectors<-framesElements$framesFirstVec
  framesSecondVectors<-framesElements$framesSecondVec
  framesThirdVectors<-framesElements$framesThirdVec
  
  # Plot LP-ds-rep
  if(plotting==TRUE){
    open3d()
    # plot skeletal points
    plot3d(medialPoints3D,type="s",radius = 0.1,col = "blue",expand = 10,box=FALSE,add = TRUE)
    #plot spine
    # plot3d(medialPoints3D[spineIndices,],type='s',radius = 0.2,col = "black",expand = 10,box=FALSE,add = TRUE)
    # plot connections
    for (i in 1:nrow(medialPoints3D)) {
      k1<-pointsIndices[i]
      k2<-parentsIndices[i]
      if(k2==k1){
        next
      }
      vectors3d(medialPoints3D[k1,],
                origin = medialPoints3D[k2,],
                # labels = c("",k1),
                headlength = 0.2,radius = 1/6, col="blue", lwd=1)
    }
    
    #frames
    for (i in frameIndices) {
      
      vectors3d(medialPoints3D[i,]+framesFirstVectors[i,],origin = medialPoints3D[i,],headlength = 0.1,radius = 1/10, col="orange", lwd=2)
      vectors3d(medialPoints3D[i,]+framesSecondVectors[i,],origin = medialPoints3D[i,],headlength = 0.1,radius = 1/10, col="orange", lwd=2)
      vectors3d(medialPoints3D[i,]+framesThirdVectors[i,],origin = medialPoints3D[i,],headlength = 0.1,radius = 1/10, col="orange", lwd=2)
    }
    
  }
  
  
  #cross-sections
  if(plotting==TRUE & cross_sections_visualization==TRUE){
    open3d()
    shade3d(tmesh, col="white",alpha=0.4)
    plot3d(spinalPoints3DComplete,type="l",lwd = 4,col = "black",expand = 10,box=FALSE,add = TRUE)
    plot3d(spinalPoints3DComplete,type = 's',radius = 0.2,col = "black",expand = 10,box=FALSE,add = TRUE)
    for (i in 1:nrow(spinalPoints3D)) {
      plot3d(slicingEllipsoids[[i]], col = "grey", alpha = 0.05, add = TRUE)
    } 
  }
  
  # ######################################################################################################
  # ######################################################################################################
  # # rotating slicing planes according to the Frenet frames of the spine
  # 
  # groupIndicesOfVeins<-list()
  # k<-4
  # i<-1
  # while (k<(numberOfMedialPoints3D-4)) {
  #   groupIndicesOfVeins[[i]]<-k:(k+6)
  #   k<-k+7
  #   i<-i+1
  # }
  # 
  # 
  # medialPoints3D_Test<-medialPoints3D
  # 
  # for (i in 1:length(groupIndicesOfVeins)) {
  #   plot3d(medialPoints3D_Test[groupIndicesOfVeins[[i]],],type = 's',radius = 0.4,col = sample(1:40,1),expand = 10,box=FALSE,add = TRUE)
  # }
  # 
  # #NB!!!! cross-sections frames are already based on identity matrix as we alighted the object at the first step
  # # thus, we rotate the viens so that they become aligned with the Frenet frames along the spine
  # 
  # framesFirstVectorsSpine<-framesFirstVectors[spineIndices,]
  # framesSecondVectorsSpine<-framesSecondVectors[spineIndices,]
  # framesThirdVectorsSpine<-framesThirdVectors[spineIndices,]
  # 
  # spinalPoints<-centeredSkel[spineIndices,]
  # 
  # #vein frames
  # veinFramesFirstVectors<-matrix(rep(c(0,0,1),length(groupIndicesOfVeins)),ncol = 3,byrow = TRUE)
  # veinFramesThirdVectors<-array(NA,dim = dim(framesThirdVectorsSpine))
  # for (i in 1:length(groupIndicesOfVeins)) {
  #   tempIndices<-groupIndicesOfVeins[[i]]
  #   tempVec<-c((medialPoints3D_Test[tempIndices[1],]-medialPoints3D_Test[tempIndices[2],])[1:2],0)
  #   veinFramesThirdVectors[i,]<-convertVec2unitVec2(tempVec)
  # }
  # 
  # 
  # for (i in 1:nrow(spinalPoints)) {
  #   vectors3d(spinalPoints[i,]+veinFramesFirstVectors[i,],origin = spinalPoints[i,],headlength = 0.1,radius = 1/10, col="orange", lwd=2)
  #   vectors3d(spinalPoints[i,]+veinFramesThirdVectors[i,],origin = spinalPoints[i,],headlength = 0.1,radius = 1/10, col="orange", lwd=2)
  #   # surface frames
  #   vectors3d(spinalPoints[i,]+framesFirstVectorsSpine[i,],origin = spinalPoints[i,],headlength = 0.1,radius = 1/10, col="black", lwd=2)
  #   vectors3d(spinalPoints[i,]+framesThirdVectorsSpine[i,],origin = spinalPoints[i,],headlength = 0.1,radius = 1/10, col="black", lwd=2)
  #   
  # }
  # 
  # 
  # for (i in 1:length(groupIndicesOfVeins)) {
  #   tempPoints<-medialPoints3D_Test[groupIndicesOfVeins[[i]],]
  #   centrPoint<-tempPoints[1,]
  #   #centering based on the spinal point
  #   centeredPoints<-tempPoints-matrix(rep(centrPoint,nrow(tempPoints)),ncol = 3,byrow = TRUE)
  #   R1<-rotMat(veinFramesFirstVectors[i,],framesFirstVectorsSpine[i,])
  #   centeredPointsRotated1<-centeredPoints%*%t(R1)
  #   R2<-rotMat(veinFramesThirdVectors[i,],framesThirdVectorsSpine[i,])
  #   centeredPointsRotated2<-centeredPointsRotated1%*%t(R2)
  #   #tranlate back
  #   newPoints<-centeredPointsRotated2+matrix(rep(centrPoint,nrow(centeredPointsRotated2)),ncol = 3,byrow = TRUE)
  #   medialPoints3D_Test[groupIndicesOfVeins[[i]],]<-newPoints
  #   
  # }
  # 
  # for (i in 1:numberOfMedialPoints3D) {
  #   if(i==numberOf2DspokePoints){
  #     next()
  #   }
  #   k1<-framesCenters[i]
  #   k2<-framesParents[i]
  #   # vectors3d(medialPoints3D_Test[k1,],origin = medialPoints3D_Test[k2,],headlength = 0.2,radius = 1/6, col="blue", lwd=1)
  #   plot3d(rbind(medialPoints3D_Test[k1,],medialPoints3D_Test[k2,]),type="l",lwd = 2,col = "blue",expand = 10,box=FALSE,add = TRUE)
  # }
  # for (i in 1:length(groupIndicesOfVeins)) {
  #   plot3d(medialPoints3D_Test[groupIndicesOfVeins[[i]],],type = 's',radius = 0.4,col = sample(1:40,1),expand = 10,box=FALSE,add = TRUE)
  # }
  # #plot mesh
  # shade3d(tmesh, col="white",alpha=0.1)
  # 
  
  ######################################################################################################
  ######################################################################################################
  # 3D spokes based as normal vectors
  
  
  # cut and stretch normal vectors in both directions
  tipOfCuttedUpSpokes<-array(NA,dim = dim(medialPoints3D))
  tempMesh4Up<-as.mesh3d(medialPoints3D[frameIndices,])
  tempMesh4Up$normals<-t(cbind(framesFirstVectors[frameIndices,],rep(1,nrow(medialPoints3D[frameIndices,]))))
  tipOfCuttedUpSpokes[frameIndices,]<-vert2points(vcgRaySearch(tempMesh4Up,mesh = tmesh))
  
  
  tipOfCuttedDownSpokes<-array(NA,dim = dim(medialPoints3D))
  tempMesh4Down<-as.mesh3d(medialPoints3D[frameIndices,])
  tempMesh4Down$normals<-t(cbind(-framesFirstVectors[frameIndices,],rep(1,nrow(medialPoints3D[frameIndices,]))))
  tipOfCuttedDownSpokes[frameIndices,]<-vert2points(vcgRaySearch(tempMesh4Down,mesh = tmesh))
  
  #plot LP_ds_rep
  if(plotting==TRUE){
    open3d()
    #spokes
    for (i in frameIndices) {
      plot3d(rbind(medialPoints3D[i,],tipOfCuttedUpSpokes[i,]),type="l",radius = 0.2,col = "blue",expand = 10,box=FALSE,add = TRUE)
    }
    for (i in frameIndices) {
      plot3d(rbind(medialPoints3D[i,],tipOfCuttedDownSpokes[i,]),type="l",radius = 0.2,col = "red",expand = 10,box=FALSE,add = TRUE)
    }
    #connections
    for (i in 1:nrow(medialPoints3D)) {
      k1<-pointsIndices[i]
      k2<-parentsIndices[i]
      if(k2==k1){
        next
      }
      vectors3d(medialPoints3D[k1,],
                origin = medialPoints3D[k2,],
                # labels = c("",k1),
                headlength = 0.2,radius = 1/6, col="blue", lwd=1)
      
      
    }
    # frames
    for (i in frameIndices) {
      
      vectors3d(medialPoints3D[i,]+framesFirstVectors[i,],origin = medialPoints3D[i,],headlength = 0.1,radius = 1/10, col="orange", lwd=2)
      vectors3d(medialPoints3D[i,]+framesSecondVectors[i,],origin = medialPoints3D[i,],headlength = 0.1,radius = 1/10, col="orange", lwd=2)
      vectors3d(medialPoints3D[i,]+framesThirdVectors[i,],origin = medialPoints3D[i,],headlength = 0.1,radius = 1/10, col="orange", lwd=2)
    }
    
    #plot mesh
    # wire3d(tmesh, col="blue")  #wire mesh
    shade3d(tmesh, col="white",alpha=0.2)
  }
  
  ######################################################################################################
  ######################################################################################################
  # 3D spokes
  
  result<-list("medialPoints3D"=medialPoints3D,
               "tipOfCuttedUpSpokes"=tipOfCuttedUpSpokes,
               "tipOfCuttedDownSpokes"=tipOfCuttedDownSpokes,
               "framesFirstVectors"=framesFirstVectors,
               "framesSecondVectors"=framesSecondVectors,
               "framesThirdVectors"=framesThirdVectors,
               "framesCenters"=framesCenters,
               "pointsIndices"=pointsIndices,
               "parentsIndices"=parentsIndices)
  
  return(result)
  
  
}
