fit_srep_By_splineAndForceField_FastApproach <-function(tmesh,
                                                        plotting=TRUE,
                                                        numberOfCoverPoints3D=1000000,
                                                        numberOfCoverPoints2D=100000,
                                                        k_Ring3D=10, # for dividing 2D and 3D object
                                                        lambda3D=0.5, # for dividing 2D and 3D object
                                                        k_Ring2D=5, # for dividing 2D and 3D object
                                                        lambda2D=1, # for dividing 2D and 3D object
                                                        sphereResolution=1, # choose 1,2,or 3 for the resolution of the urchin
                                                        circleRsolution=24, # resolution of the 2D urchin
                                                        urchinRadius=0.5, # resolution of the 2D urchin
                                                        thresholdDistance2D=0.2,
                                                        kRing4Normal=2,
                                                        polyDegree3D=4,
                                                        polyDegree2D=5,
                                                        alpha1=2, # for alpha convex hull
                                                        numberOfPoints4alphaHull=5000,
                                                        numberOf2DspokePoints=4,
                                                        numberOfRandomPoints=10, #number of random points for each triangle of the medial surface
                                                        numberOfSpanialPoints=27,
                                                        rotationGap=10, #to fit best medial curve
                                                        increase3DBoundaryPoint=TRUE # to generate 3D spokes
){
  
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

  
  # numberOfVertices<-dim(tmesh$vb)[2]
  # geoDistancesMatrix<-array(NA,dim = c(numberOfVertices,numberOfVertices))
  # for (i in 1:numberOfVertices) {
  #   geoDistancesMatrix[i,]<-vcgDijkstra(tmesh,i)
  # }
  # 
  # #make geoDistancesMatrix symmetric
  # geoDistancesMatrix<-forceSymmetric(geoDistancesMatrix)
  # # geoDistancesMatrix<-upper.tri(geoDistancesMatrix,diag = TRUE)+t(upper.tri(geoDistancesMatrix,diag = FALSE))
  # 
  # # mean curvature of the vertices
  # vcgInfo<-vcgCurve(tmesh)
  # meanCurvatureOfAllVertices<-vcgInfo$meanvb
  # 
  # #adjancy matrix of neighbors
  # allAdjacentVertices<-vcgVertexNeighbors(tmesh,numstep = k_Ring3D)
  # 
  # unitNormals<-t(tmesh$normals[1:3,])/sqrt(rowSums(t(tmesh$normals[1:3,])^2))
  # 
  # normalMultiplicationMatrix<-unitNormals%*%t(unitNormals)
  # 
  # elementWiseMultipleMatrix<-exp(lambda3D*normalMultiplicationMatrix*geoDistancesMatrix)
  # 
  # #compute affinity matrix W
  # W<-elementWiseMultipleMatrix
  # numberOfPoints<-dim(W)[2]
  # tempVec<-1:numberOfPoints
  # pb <- txtProgressBar(style = 3)
  # for (i in 1:numberOfPoints) {
  #   setTxtProgressBar(pb,i/numberOfPoints)
  #   notNeighbor<-tempVec[!tempVec %in% allAdjacentVertices[[i]]]
  #   W[i,notNeighbor]<-0
  # }
  # close(pb)
  # 
  # d<-(rowSums(W))^(-1/2)
  # D<-diag(d)
  # 
  # L <- D%*%W%*%D
  # 
  # L <- diag(nrow(L))-L
  # 
  # eigenValAndVec<-eigen(L)
  # 
  # smallestEigenVectorIndex<-nrow(L)-1
  # secondSmallestEigenValue<-eigenValAndVec$values[smallestEigenVectorIndex]
  # secondSmallestEigenVector<-eigenValAndVec$vectors[,smallestEigenVectorIndex]
  # 
  # upIndices<-which(secondSmallestEigenVector>0)
  # downIndices<-which(secondSmallestEigenVector<=0)
  # 
  # 
  # #plot
  # if(plotting==TRUE){
  #   open3d()
  #   spheres3d(vert2points(tmesh)[upIndices,],radius = 2,col='orange')
  #   spheres3d(vert2points(tmesh)[downIndices,],radius = 2,col='blue')
  # }
  # 
  # # we consider the top and bottom parts of the object as two meshes by
  # # tmesh$it which is the face indices
  # #Up
  # tempMatrix<-matrix(t(tmesh$it) %in% upIndices,ncol = 3,byrow = FALSE)
  # polyIndicesUp<-which(rowSums(tempMatrix)==3)
  # trglsUp <- as.matrix(t(t(tmesh$it)[polyIndicesUp,]))
  # tmeshUp <- tmesh3d(verts, trglsUp)
  # #Edge
  # polyIndicesEdge<-which(rowSums(tempMatrix)==2 | rowSums(tempMatrix)==1)
  # trglsEdge <- as.matrix(t(t(tmesh$it)[polyIndicesEdge,]))
  # tmeshEdge <- tmesh3d(verts, trglsEdge)
  # #Down
  # tempMatrix<-matrix(t(tmesh$it) %in% downIndices,ncol = 3,byrow = FALSE)
  # polyIndicesDown<-which(rowSums(tempMatrix)==3)
  # trglsDown <- as.matrix(t(t(tmesh$it)[polyIndicesDown,]))
  # tmeshDown <- tmesh3d(verts, trglsDown)
  # 
  # 
  # #remove unreferenced vertices
  # tmeshUp<-vcgClean(tmeshUp,sel = 1)
  # tmeshEdge<-vcgClean(tmeshEdge,sel = 1)
  # tmeshDown<-vcgClean(tmeshDown,sel = 1)
  # 
  # #plot
  # if(plotting==TRUE){
  #   open3d()
  #   shade3d(tmeshUp,col="blue")
  #   shade3d(tmeshEdge,col="yellow")
  #   shade3d(tmeshDown,col="red")
  # }
  
  
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
  
  # #make the mesh slightly smaller
  # smallMesh<-tmesh
  # smallMesh$vb[1:3,]<-smallMesh$vb[1:3,]-0.4*smallMesh$normals[1:3,]
  
  # #find internal points
  # insidePoints<-pip3d(Vertices = t(smallMesh$vb[1:3,]),
  #                     Faces = t(smallMesh$it),
  #                     Queries = coverpoint)
  
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
  
  
  # pb <- txtProgressBar(style = 3)
  # force<-array(NA,nrow(innerPoints))
  # shortestLengths<-array(NA,nrow(innerPoints))
  # # angles<-array(NA,nrow(innerPoints))
  # for (k in 1:nrow(innerPoints)) {
  #   setTxtProgressBar(pb,k/nrow(innerPoints))
  #   
  #   center<-innerPoints[k,]
  #   
  #   tempSphere<-makeSphereMesh(center = center,radius = urchinRadius,subdivision = sphereResolution)
  #   
  #   # #plot
  #   # shade3d(tmesh,col="white",alpha=0.1)
  #   # shade3d(tempSphere,col="blue")
  #   
  #   #find the intercetions of rays with internal normal directions
  #   intersectionsUp <- vcgRaySearch(tempSphere,mesh = tmeshUp,mindist = TRUE)
  #   intersectionsDown <- vcgRaySearch(tempSphere,mesh = tmeshDown,mindist = TRUE)
  #   
  #   tipOfCuttedSpokesUp<-vert2points(intersectionsUp)[intersectionsUp$quality==1,]
  #   tipOfCuttedSpokesDown<-vert2points(intersectionsDown)[intersectionsDown$quality==1,]
  #   
  #   if(!(is.matrix(tipOfCuttedSpokesUp) & is.matrix(tipOfCuttedSpokesDown))){
  #     
  #     force[k]<-10^4
  #     shortestLengths[k]<-0
  #     
  #     next
  #   }
  #   
  #   # #plot
  #   # if(plotting==TRUE){
  #   #   open3d()
  #   #   # spheres3d(vert2points(intersections),col="blue",radius = 0.2) #plot intersections
  #   #   for (i in 1:dim(tipOfCuttedSpokesUp)[1]) {
  #   #     plot3d(rbind(center,tipOfCuttedSpokesUp[i,]),type="l",col = "blue",expand = 10,box=FALSE,add = TRUE)
  #   #   }
  #   #   for (i in 1:dim(tipOfCuttedSpokesDown)[1]) {
  #   #     plot3d(rbind(center,tipOfCuttedSpokesDown[i,]),type="l",col = "red",expand = 10,box=FALSE,add = TRUE)
  #   #   }
  #   #   #NB!!! we have the information of normals at intersections intersections$normals !!!!
  #   #   shade3d(tmesh, col="white",alpha=0.2)  #surface mech
  #   # }
  #   
  #   tipOfCuttedSpokes<-rbind(tipOfCuttedSpokesUp,tipOfCuttedSpokesDown)
  #   labels<-c(rep(1,nrow(tipOfCuttedSpokesUp)),
  #             rep(2,nrow(tipOfCuttedSpokesDown)))
  #   
  #   vectorsMatrix<-tipOfCuttedSpokes-matrix(rep(center,nrow(tipOfCuttedSpokes)),ncol = 3,byrow = TRUE)
  #   
  #   # spokeLengths<-apply(X = vectorsMatrix, MARGIN = 1,FUN = myNorm )
  #   # spokesUnitDirections<-t(apply(X = vectorsMatrix, MARGIN = 1,FUN = convertVec2unitVec ))
  #   
  #   forceTemp<-forceFunctionAngleBased_4MultiObject(vectorsMatrix = vectorsMatrix,
  #                                                   labels = labels,
  #                                                   type = 'one')
  #   
  #   # forceTemp<-forceFunctionAngleBased(vectorsMatrix = vectorsMatrix)
  #   
  #   force[k]<-forceTemp$force
  #   
  #   shortestLengths[k]<-forceTemp$shortestLength
  #   
  #   # angles[k]<-urchin_AngleOfShortestSpokes(vectorsMatrix = vectorsMatrix)
  #   
  #   # angles[k]<-urchin_AngleOfShortestSpokes_4MultiObject(vectorsMatrix = vectorsMatrix,
  #   #                                           labels = labels,
  #   #                                           type = 'one')
  #   
  #   
  # }
  # close(pb)
  # 
  # forceVectorsMagnitudes<-force
  # # forceVectorsMagnitudes<-angles
  # # plot(1:length(forceVectorsMagnitudes),forceVectorsMagnitudes,type = 'p')
  # 
  # roundedForce<-round(forceVectorsMagnitudes,digits = 2)
  # uniqueRounded<-sort(unique(roundedForce))
  # table(roundedForce)
  # barplot(uniqueRounded,roundedForce)
  # distanceThreshold<-2*quantile(shortestLengths,0.1)
  # for (i in length(uniqueRounded):1) {
  #   open3d()
  #   shade3d(tmesh,col="white",alpha=0.1)
  #   plot3d(innerPoints[(roundedForce %in% uniqueRounded[1:i]) &
  #                        shortestLengths>distanceThreshold,],type="p",col = "blue",expand = 10,box=FALSE,add = TRUE)
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
  
  # edgePoints<-c()
  # for (i in 1:dim(t(tmeshEdge$it))[1]) {
  #   edgePoints<-rbind(edgePoints,
  #                     colMeans((t(tmeshEdge$vb)[t(tmeshEdge$it)[i,],1:3])))
  # }
  
  # if(plotting==TRUE){
  #   open3d()
  #   # plot3d(selectedMiddlePoint,type="p",col = "blue",expand = 10,box=FALSE,add = TRUE)
  #   plot3d(edgePoints,type="p",col = "red",expand = 10,box=FALSE,add = TRUE)
  #   shade3d(tmesh,col="white",alpha=0.2)
  # }
  
  # allSelectedPoints<-rbind(selectedMiddlePoint,edgePoints)
  
  # x y z
  allSelectedPoints<-vert2points(tmesh)
  
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
  
  a<-1
  tempPoint1<-boundaryPoints1[a,]
  tempPoint2<-boundaryPoints2[a,]
  mesh2DPoints<-rbind(tempPoint1,tempPoint2)
  boundaryPoints1<-boundaryPoints1[-a,]
  boundaryPoints2<-boundaryPoints2[-a,]
  while (is.matrix(boundaryPoints1) & is.matrix(boundaryPoints2)) {
    if(sum(tempPoint1 %in% boundaryPoints1)){
      a<-which(boundaryPoints1[,1]==tempPoint1[1])[1]
      tempPoint1<-boundaryPoints1[a,]
      tempPoint2<-boundaryPoints2[a,]
      mesh2DPoints<-rbind(mesh2DPoints,tempPoint1,tempPoint2)
      boundaryPoints1<-boundaryPoints1[-a,]
      boundaryPoints2<-boundaryPoints2[-a,]
    }else if(sum(tempPoint1 %in% boundaryPoints2)){
      a<-which(boundaryPoints2[,1]==tempPoint1[1])[1]
      tempPoint1<-boundaryPoints1[a,]
      tempPoint2<-boundaryPoints2[a,]
      mesh2DPoints<-rbind(mesh2DPoints,tempPoint2,tempPoint1)
      boundaryPoints1<-boundaryPoints1[-a,]
      boundaryPoints2<-boundaryPoints2[-a,]
    }else if(sum(tempPoint2 %in% boundaryPoints1)){
      a<-which(boundaryPoints1[,1]==tempPoint2[1])[1]
      tempPoint1<-boundaryPoints1[a,]
      tempPoint2<-boundaryPoints2[a,]
      mesh2DPoints<-rbind(mesh2DPoints,tempPoint1,tempPoint2)
      boundaryPoints1<-boundaryPoints1[-a,]
      boundaryPoints2<-boundaryPoints2[-a,]
    }else if(sum(tempPoint2 %in% boundaryPoints2)){
      a<-which(boundaryPoints2[,1]==tempPoint2[1])[1]
      tempPoint1<-boundaryPoints1[a,]
      tempPoint2<-boundaryPoints2[a,]
      mesh2DPoints<-rbind(mesh2DPoints,tempPoint1,tempPoint2)
      boundaryPoints1<-boundaryPoints1[-a,]
      boundaryPoints2<-boundaryPoints2[-a,]
    }else{
      break
    }
  }
  mesh2DPoints<-rbind(mesh2DPoints,boundaryPoints1,boundaryPoints2)
  
  if(plotting==TRUE){
    plot(mesh2DPoints,type="l",xlim = xlim,ylim = xlim,xlab = "",ylab = "") 
  }
  
  # mesh2D must be open
  mesh2D<-mesh2DPoints[-nrow(mesh2DPoints),]
  
  #remove duplicated vertices
  duplicatedIndices<-which(duplicated(mesh2D[,1])+duplicated(mesh2D[,2])==2)
  mesh2D<-mesh2D[-duplicatedIndices,]
  
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
  acceptableEdges<-sort(which.maxn(edgeLengths,n = 100))
  mesh2D<-mesh2D[acceptableEdges,]
  
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
  
  unitNormals<-verticesNormals_mesh2D(mesh2D = mesh2D,kRing4Normal=kRing4Normal)
  
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
  
  # edgeSuspects<-which(round(abs(secondSmallestEigenVector),digits = 2)==0)
  
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

  
  if(plotting==TRUE){
    plot(rbind(mesh2D,mesh2D[1,]),type = "l",xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(mesh2D[c(startUp,endUp),],pch=16,xlim = xlim,ylim = xlim,xlab = "",ylab = "")
  }
  
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
  
  
  if(TRUE){
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
    tempCircle<-sphereGenerator_2D(center = center,r = urchinRadius,n = circleRsolution,asymmetric = TRUE)
    
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
      
      tipOfCuttedSpokesWithLabels<-rbind(cbind(tipOfCuttedUpSpokes[,1:2],rep(1,nrow(tipOfCuttedUpSpokes))),
                                         cbind(tipOfCuttedDownSpokes[,1:2],rep(2,nrow(tipOfCuttedDownSpokes))))
      
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
  
  if (plotting==TRUE) {
    barplot(uniqueRounded,roundedForce) 
  }
  
  #plot levels
  if(plotting==TRUE){
    for (i in length(uniqueRounded):1) {
      plot(mesh2D,type = "l",xlim = xlim,ylim = xlim,xlab = "",ylab = "")
      cumLevel_ith_Points<-innerPoints2D[roundedForce %in% uniqueRounded[1:i],]
      par(new=TRUE)
      plot(cumLevel_ith_Points,pch=20,col='blue',xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    }    
  }

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
  if(TRUE){
    plot(rbind(mesh2D,mesh2D[1,]),type = "l",xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(cumLevel_ith_Points,pch=20,col="blue",xlim = xlim,ylim = xlim,xlab = "",ylab = "") 
  }
  
  centroid2D<-colMeans(cumLevel_ith_Points)
  centeredPoints2D<-scale(cumLevel_ith_Points,center = TRUE,scale = FALSE)
  
  boundaryPoints<-mesh2D
  
  
  x4plot<-seq(from = xlim[1],to = xlim[2],length.out = 100) 
  
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
      par(new=TRUE)
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
  # x y z
  x<-selectedMiddlePoint2D[,1]
  y<-selectedMiddlePoint2D[,2]
  
  fit4_2D <- lm(y ~ poly(x, degree = polyDegree2D,raw = TRUE),
                data=as.data.frame(cbind(x,y))) #NB!!! raw=F calculate orthogonal polynomial
  
  
  x<-seq(from=min(selectedMiddlePoint2D[,1]),to=max(selectedMiddlePoint2D[,1]),length.out = 1000)
  xData1<-data.frame(x=x)
  curvePoints1<-predict(fit4_2D, newdata = xData1)
  
  medialPoints_2D<-cbind(x,curvePoints1)
  
  # plot(rbind(boundaryPointsTemp,boundaryPointsTemp[1,]),type="l",xlim = xlim,ylim = xlim,xlab = "",ylab = "") 
  # par(new=TRUE) 
  # plot(medialPoints_2D,pch=20,col="blue",xlim = xlim,ylim = xlim,xlab = "",ylab = "") 
  
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
  
  # plot(rbind(boundaryPointsTemp,boundaryPointsTemp[1,]),type="l",xlim = xlim,ylim = xlim,xlab = "",ylab = "") 
  # par(new=TRUE) 
  # plot(medialPointsEqualyDis_2D,pch=20,col="blue",xlim = xlim,ylim = xlim,xlab = "",ylab = "") 
  
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
  
  if(plotting==TRUE){
    plot(rbind(boundaryPoints,boundaryPoints[1,]),type="l",xlim = xlim,ylim = xlim,xlab = "",ylab = "") 
    par(new=TRUE) 
    plot(medialPointsEqualyDis_2D,pch=20,col="blue",xlim = xlim,ylim = xlim,xlab = "",ylab = "")  
  }
  
  #store and remove first and the last point
  firstPoint<-medialPointsEqualyDis_2D[1,]
  lastPoint<-medialPointsEqualyDis_2D[nrow(medialPointsEqualyDis_2D),]
  medialPointsEqualyDis_2D<-medialPointsEqualyDis_2D[-c(1,nrow(medialPointsEqualyDis_2D)),]
  
  if(plotting==TRUE){
    plot(rbind(boundaryPoints,boundaryPoints[1,]),type="l",xlim = xlim,ylim = xlim,xlab = "",ylab = "") 
    par(new=TRUE) 
    plot(medialPointsEqualyDis_2D,pch=20,col="blue",xlim = xlim,ylim = xlim,xlab = "",ylab = "") 
    par(new=TRUE) 
    plot(rbind(firstPoint,lastPoint),pch=20,col="red",xlim = xlim,ylim = xlim,xlab = "",ylab = "") 
  }

  spokesLengths<-rep(NA,nrow(medialPointsEqualyDis_2D))
  tipOfUpSpokes<-array(NA,dim=c(nrow(medialPointsEqualyDis_2D),2))
  tipOfDownSpokes<-array(NA,dim=c(nrow(medialPointsEqualyDis_2D),2))
  pb <- txtProgressBar(style = 3)
  for (k in 1:nrow(medialPointsEqualyDis_2D)) {
    setTxtProgressBar(pb,k/nrow(medialPointsEqualyDis_2D))
    
    center<-medialPointsEqualyDis_2D[k,]
    
    # we use asymmetric circles to avoid antipodal vectors
    tempCircle<-sphereGenerator_2D(center = center,r = 1,n = 50,asymmetric = TRUE)
    
    # #plot
    # plot(rbind(mesh2D,mesh2D[1,]),type = "l",xlim = xlim,ylim = ylim)
    # par(new=TRUE)
    # plot(tempCircle,pch=20,xlim = xlim,ylim = ylim)
    
    # plot(rbind(tempCircle,tempCircle[1,]),type = "l",xlim = center+c(-1.5,1.5),ylim = center+c(-1.5,1.5))
    # for (i in 1:nrow(tempCircle)) {
    #   par(new=TRUE)
    #   plot(rbind(tempCircle[i,],center),type = "l",xlim = center+c(-1.5,1.5),ylim = center+c(-1.5,1.5))
    # }
    
    allMeshes2D<-list(rbind(boundaryPoints[upIndices,],boundaryPoints[downIndices[1],]),
                      rbind(boundaryPoints[downIndices,],boundaryPoints[upIndices[1],]))
    
    
    tipOfCuttedSpokesWithLabels<-cutAndStretchSpokes_4Bivalves_2D(allSpokesTips = tempCircle,
                                                                  allSpokesTails = matrix(rep(center,dim(tempCircle)[1]),ncol = 2,byrow = TRUE),
                                                                  allMeshes2D)
    
    # plot(allMeshesTemp,pch=20,col='blue',xlim = xlim,ylim = ylim,xlab = "",ylab = "")
    
    tipOfCuttedSpokes<-tipOfCuttedSpokesWithLabels[,1:2]
    
    # 
    # plot(rbind(boundaryPoints,boundaryPoints[1,]),type = "l",xlim = xlim,ylim = ylim,xlab = "",ylab = "")
    # par(new=TRUE)
    # for (i in 1:dim(tipOfCuttedSpokes)[1]) {
    #   plot(rbind(center,tipOfCuttedSpokes[i,]),type = "l",xlim = xlim,ylim = ylim,xlab = "",ylab = "")
    #   par(new=TRUE)
    # }
    
    tempMatrix1<-matrix(rep(center,nrow(tipOfCuttedSpokes)),ncol = 2,byrow = TRUE)
    vectorsMatrix<-tipOfCuttedSpokes-tempMatrix1
    
    labels<-tipOfCuttedSpokesWithLabels[,3]
    
    vectorsMatrixUp<-vectorsMatrix[as.vector(which(labels==1)),]
    vectorsMatrixDown<-vectorsMatrix[as.vector(which(labels==2)),]
    
    vectorsLengthsUp<-apply(X = vectorsMatrixUp, MARGIN = 1,FUN = myNorm)
    vectorsLengthsDown<-apply(X = vectorsMatrixDown, MARGIN = 1,FUN = myNorm)
    
    shortestVectorUpIndex<-which.min(vectorsLengthsUp)
    shortestVectorDownIndex<-which.min(vectorsLengthsDown)
    
    tipOfCuttedSpokesUp<-tipOfCuttedSpokes[which(labels==1),]
    tipOfCuttedSpokesDown<-tipOfCuttedSpokes[which(labels==2),]
    
    
    shortestSpokeUp<-tipOfCuttedSpokesUp[shortestVectorUpIndex,]
    shortestSpokeDown<-tipOfCuttedSpokesDown[shortestVectorDownIndex,]
    
    #average length of top and bottom spokes
    spokesLengths[k]<-(min(vectorsLengthsUp)+min(vectorsLengthsDown))/2
    
    tipOfUpSpokes[k,]<-shortestSpokeUp
    tipOfDownSpokes[k,]<-shortestSpokeDown
  }
  close(pb)
  
  #plot
  if(plotting==TRUE){
    plot(rbind(boundaryPoints,boundaryPoints[1,]),type = "l",xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    for (i in 1:nrow(medialPointsEqualyDis_2D)) {
      par(new=TRUE)
      plot(rbind(medialPointsEqualyDis_2D[i,],tipOfUpSpokes[i,]),col='blue',lwd = 2,type = "l",xlim = xlim,ylim = xlim,xlab = "",ylab = "")
      par(new=TRUE)
      plot(rbind(medialPointsEqualyDis_2D[i,],tipOfDownSpokes[i,]),col='red',lwd = 2,type = "l",xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    } 
  }
  
  #N vector (i.e., normal) of medial point
  N_Vec<-c()
  k<-2
  for (i in 2:(nrow(medialPointsEqualyDis_2D)-1)) {
    N_Vec<-rbind(N_Vec,normalOfaVertex(medialPointsEqualyDis_2D[i+1,],
                                       medialPointsEqualyDis_2D[i,],
                                       medialPointsEqualyDis_2D[i-1,]))
    k<-k+1
  }
  N_Vec<-rbind(N_Vec[1,],N_Vec,N_Vec[nrow(N_Vec),])
  
  # T vector (i.e., tangent) of medial points
  rotation4Tangent<-rotMat(c(0,1),c(-1,0))
  T_Vec<-array(NA,dim = dim(N_Vec))
  for (i in 1:nrow(N_Vec)) {
    T_Vec[i,]<-N_Vec[i,]%*%t(rotation4Tangent)
  }
  
  if(plotting==TRUE){
    plot(rbind(boundaryPoints,boundaryPoints[1,]),type = "l",xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    for (i in 1:nrow(medialPointsEqualyDis_2D)) {
      par(new=TRUE)
      plot(rbind(medialPointsEqualyDis_2D[i,],medialPointsEqualyDis_2D[i,]+N_Vec[i,]),col='blue',lwd = 2,type = "l",xlim = xlim,ylim = xlim,xlab = "",ylab = "")
      par(new=TRUE)
      plot(rbind(medialPointsEqualyDis_2D[i,],medialPointsEqualyDis_2D[i,]+T_Vec[i,]),col='red',lwd = 2,type = "l",xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    } 
  }
  
  allPoints<-c()
  allCircles<-list()
  for (i in 1:nrow(medialPointsEqualyDis_2D)) {
    centerTemp<-medialPointsEqualyDis_2D[i,]
    radiusTemp<-spokesLengths[i]
    tempCircle<-sphereGenerator_2D(center = centerTemp,r = radiusTemp,n = circleRsolution,asymmetric = TRUE)
    allCircles[[i]]<-tempCircle
    allPoints<-rbind(allPoints,tempCircle)
  }
  
  if(plotting==TRUE){
    plot(rbind(boundaryPoints,boundaryPoints[1,]),type = "l",lty=3,xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(medialPointsEqualyDis_2D,col='blue',pch=20,xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(allPoints,pch='.',xlim = xlim,ylim = xlim,xlab = "",ylab = "") 
  }
  
  matrixInOut<-array(1,dim=c(nrow(allPoints),nrow(medialPointsEqualyDis_2D)))
  for (j in 1:nrow(medialPointsEqualyDis_2D)) {
    matrixInOut[,j]<-pip2d(allCircles[[j]],allPoints)
  }
  labelsFinal<-rep(NA,nrow(allPoints))
  for (i in 1:nrow(matrixInOut)) {
    labelsFinal[i]<- 1 %in% matrixInOut[i,]
  }
  
  if(plotting==TRUE){
    plot(rbind(boundaryPoints,boundaryPoints[1,]),type = "l",lty=3,xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(medialPointsEqualyDis_2D,col='blue',pch=20,xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(allPoints[!labelsFinal,],pch='.',xlim = xlim,ylim = xlim,xlab = "",ylab = "") 
  }
  
  # #remove head and tail circles
  # matrixInOut<-array(NA,dim=c(nrow(allPoints),nrow(medialPointsEqualyDis_2D)))
  # for (j in 1:nrow(medialPointsEqualyDis_2D)) {
  #   matrixInOut[,j]<-pip2d(allCircles[[j]],allPoints)
  # }
  # labelsFinal<-rep(NA,nrow(allPoints))
  # for (i in 1:nrow(matrixInOut)) {
  #   labelsFinal[i]<- 1 %in% matrixInOut[i,]
  # }
  # firstCircleBoundary<-allCircles[[1]][!labelsFinal[1:nrow(allCircles[[1]])],]
  # tempVec2<-labelsFinal[(length(labelsFinal) - nrow(allCircles[[1]])+1):length(labelsFinal)]
  # lastCircleBoundary<-allCircles[[nrow(medialPointsEqualyDis_2D)]][!tempVec2,]
  # 
  # labelsFinal[1:nrow(allCircles[[1]])]<-1
  # labelsFinal[(length(labelsFinal) - nrow(allCircles[[1]])):length(labelsFinal)]<-1
  
  # plot(rbind(boundaryPoints,boundaryPoints[1,]),type = "l",lty=3,xlim = xlim,ylim = xlim,xlab = "",ylab = "")
  # par(new=TRUE)
  # plot(allPoints[!labelsFinal,],pch='.',xlim = xlim,ylim = xlim,xlab = "",ylab = "")
  # par(new=TRUE)
  # plot(medialPointsEqualyDis_2D,col='blue',pch = 20,xlim = xlim,ylim = xlim,xlab = "",ylab = "")
  # par(new=TRUE)
  # plot(firstCircleBoundary,col='green',xlim = xlim,ylim = xlim,xlab = "",ylab = "")
  # par(new=TRUE)
  # plot(lastCircleBoundary,col='green',xlim = xlim,ylim = xlim,xlab = "",ylab = "")
  
  
  # radiusFirstCircle<-norm(medialPointsEqualyDis_2D[1,]-firstCircleBoundary[1,],type = '2')
  # extraSpokeSpineHead<-radiusFirstCircle*convertVec2unitVec2(colMeans(firstCircleBoundary-
  # matrix(rep(medialPointsEqualyDis_2D[1,],nrow(firstCircleBoundary)),ncol = 2,byrow = TRUE)))
  
  extraSpokeSpineHead<-medialPointsEqualyDis_2D[1,]
  extraSpokeSpineHead_Tip<-firstPoint
  
  # radiusLastCircle<-norm(medialPointsEqualyDis_2D[nrow(medialPointsEqualyDis_2D),]-lastCircleBoundary[1,],type = '2')
  # extraSpokeSpineTail<-radiusLastCircle*convertVec2unitVec2(colMeans(lastCircleBoundary-
  # matrix(rep(medialPointsEqualyDis_2D[nrow(medialPointsEqualyDis_2D),],nrow(lastCircleBoundary)),ncol = 2,byrow = TRUE)))
  
  extraSpokeSpineTail<-medialPointsEqualyDis_2D[nrow(medialPointsEqualyDis_2D),]
  extraSpokeSpineTail_Tip<-lastPoint
  
  
  # x<-c(extraSpokeSpineHead_Tip[1],medialPointsEqualyDis_2D[,1],extraSpokeSpineTail_Tip[1])
  # y<-c(extraSpokeSpineHead_Tip[2],medialPointsEqualyDis_2D[,2],extraSpokeSpineTail_Tip[2])
  # plot(x,y,pch=20,xlim = xlim,ylim = ylim,xlab = "",ylab = "")
  # fit5_2D <- lm(y ~ poly(x , degree = polyDegree2D,raw = TRUE),
  #               data=as.data.frame(cbind(x,y))) #NB!!! raw=F calculate orthogonal polynomial
  # 
  # xData2<-data.frame(x=x)
  # curvePoints2<-predict(fit5_2D, newdata = xData2)
  # 
  # plot(rbind(boundaryPoints,boundaryPoints[1,]),type = "l",xlim = xlim,ylim = ylim,xlab = "",ylab = "")
  # par(new=TRUE)
  # plot(allPoints[!labelsFinal,],xlim = xlim,ylim = ylim,xlab = "",ylab = "")
  # par(new=TRUE)
  # plot(medialPointsEqualyDis_2D,col='blue',pch = 20,xlim = xlim,ylim = ylim,xlab = "",ylab = "")
  # par(new=TRUE)
  # plot(x,as.vector(curvePoints2),pch=20,col='orange',xlim = xlim,ylim = ylim,xlab = "",ylab = "")
  # 
  # new2D_mesh<-allPoints[!labelsFinal,]
  # xData4<-data.frame(x=new2D_mesh[,1])
  # curvePointsY4<-predict(fit5_2D, newdata = xData4)
  # medialPointstest<-cbind(new2D_mesh[,1],curvePointsY4)
  # 
  # #plot 2D
  # plot(new2D_mesh,col="black",pch=20,cex=0.2,xlim = xlim,ylim = xlim,xlab = "",ylab = "")
  # par(new=TRUE)
  # plot(medialPointstest,col="green",pch=20,cex=0.2,xlim = xlim,ylim = xlim,xlab = "",ylab = "")
  # par(new=TRUE)
  # plot(medialPointsEqualyDis_2D,col="red",pch=20,cex=0.2,xlim = xlim,ylim = xlim,xlab = "",ylab = "")
  # 
  # 
  # topPart<-new2D_mesh[new2D_mesh[,2]>curvePointsY4,]
  # bottomPart<-new2D_mesh[new2D_mesh[,2]<curvePointsY4,]
  
  new2D_mesh<-allPoints[!labelsFinal,]
  
  topPart<-c()
  bottomPart<-c()
  for (i in 1:nrow(new2D_mesh)) {
    tempBoundaryPoint<-new2D_mesh[i,]
    tempMatrix<-matrix(rep(tempBoundaryPoint,nrow(medialPointsEqualyDis_2D)),ncol = 2,byrow = TRUE)
    distancesTemp<-apply(X = medialPointsEqualyDis_2D-tempMatrix, MARGIN = 1,FUN = myNorm)
    indexTemp<-which.min(distancesTemp)
    closestMedialPoint<-medialPointsEqualyDis_2D[indexTemp,]
    unitDirectionTemp<-convertVec2unitVec2(tempBoundaryPoint-closestMedialPoint)
    geoDis<-geodesicDistance(unitDirectionTemp,N_Vec[indexTemp,])
    if(geoDis<=pi/2){
      topPart<-rbind(topPart,tempBoundaryPoint)
    }else{
      bottomPart<-rbind(bottomPart,tempBoundaryPoint)
    }
  }
  
  
  closestMedialLabelTop<-rep(NA,nrow(topPart))
  for (i in 1:nrow(topPart)) {
    tempMatrix<-matrix(rep(topPart[i,],nrow(medialPointsEqualyDis_2D)),ncol = 2,byrow = TRUE)
    distancesTemp<-apply(X = medialPointsEqualyDis_2D-tempMatrix, MARGIN = 1,FUN = myNorm)
    closestMedialLabelTop[i]<-which.min(distancesTemp)
  }
  closestMedialLabelBottom<-rep(NA,nrow(topPart))
  for (i in 1:nrow(bottomPart)) {
    tempMatrix<-matrix(rep(bottomPart[i,],nrow(medialPointsEqualyDis_2D)),ncol = 2,byrow = TRUE)
    distancesTemp<-apply(X = medialPointsEqualyDis_2D-tempMatrix, MARGIN = 1,FUN = myNorm)
    closestMedialLabelBottom[i]<-which.min(distancesTemp)
  }
  
  if(plotting==TRUE){
    plot(rbind(boundaryPoints,boundaryPoints[1,]),type = "l",lty=3,xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    for (i in 1:nrow(topPart)) {
      plot(rbind(topPart[i,],medialPointsEqualyDis_2D[closestMedialLabelTop[i],]),col="blue",type = 'l',xlim = xlim,ylim = xlim,xlab = "",ylab = "")
      par(new=TRUE)
    }
    plot(medialPointsEqualyDis_2D,col="black",pch=20,xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    for (i in 1:nrow(bottomPart)) {
      plot(rbind(bottomPart[i,],medialPointsEqualyDis_2D[closestMedialLabelBottom[i],]),col="red",type = 'l',xlim = xlim,ylim = xlim,xlab = "",ylab = "")
      par(new=TRUE)
    }
    for (i in 1:nrow(topPart)) {
      plot(rbind(topPart[i,],medialPointsEqualyDis_2D[closestMedialLabelTop[i],]),col="blue",type = 'l',xlim = xlim,ylim = xlim,xlab = "",ylab = "")
      par(new=TRUE)
    }
  }
  
  upSpokesMeans<-array(NA,dim = c(nrow(medialPointsEqualyDis_2D),2))
  for (i in 1:nrow(medialPointsEqualyDis_2D)) {
    if(sum(topPart[which(closestMedialLabelTop==i),])==0){
      next
    }
    tempBoundary<-topPart[which(closestMedialLabelTop==i),]
    tempSpokes<-tempBoundary-matrix(rep(medialPointsEqualyDis_2D[i,],length(tempBoundary)/2),ncol = 2,byrow = TRUE)
    upSpokesMeans[i,]<-colMeans(tempSpokes)
  }
  downSpokesMeans<-array(NA,dim = c(nrow(medialPointsEqualyDis_2D),2))
  for (i in 1:nrow(medialPointsEqualyDis_2D)) {
    if(sum(bottomPart[which(closestMedialLabelBottom==i),])==0){
      next
    }
    tempBoundary<-bottomPart[which(closestMedialLabelBottom==i),]
    tempSpokes<-tempBoundary-matrix(rep(medialPointsEqualyDis_2D[i,],length(tempBoundary)/2),ncol = 2,byrow = TRUE)
    downSpokesMeans[i,]<-colMeans(tempSpokes)
  }
  
  if(plotting==TRUE){
    for (i in 1:nrow(medialPointsEqualyDis_2D)) {
      plot(rbind(upSpokesMeans[i,]+medialPointsEqualyDis_2D[i,],
                 medialPointsEqualyDis_2D[i,]),col="blue",type = 'l',xlim = xlim,ylim = xlim,xlab = "",ylab = "")
      par(new=TRUE)
    }
    plot(rbind(boundaryPoints,boundaryPoints[1,]),type = "l",lty=3,xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(medialPointsEqualyDis_2D,col="black",pch=20,xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    for (i in 1:nrow(medialPointsEqualyDis_2D)) {
      plot(rbind(downSpokesMeans[i,]+medialPointsEqualyDis_2D[i,],
                 medialPointsEqualyDis_2D[i,]),col="red",type = 'l',xlim = xlim,ylim = xlim,xlab = "",ylab = "")
      par(new=TRUE)
    } 
  }
  
  #Make sure the first point has up and down spokes
  if(sum(is.na(upSpokesMeans[1,]))>0){
    tempMatrix<-matrix(rep(medialPointsEqualyDis_2D[1,],nrow(topPart)),ncol = 2,byrow = TRUE)
    distancesTemp<-apply(X = topPart-tempMatrix, MARGIN = 1,FUN = myNorm)
    upSpokesMeans[1,]<-topPart[which.min(distancesTemp),]-medialPointsEqualyDis_2D[1,]
  }
  if(sum(is.na(downSpokesMeans[1,]))>0){
    tempMatrix<-matrix(rep(medialPointsEqualyDis_2D[1,],nrow(bottomPart)),ncol = 2,byrow = TRUE)
    distancesTemp<-apply(X = bottomPart-tempMatrix, MARGIN = 1,FUN = myNorm)
    downSpokesMeans[1,]<-bottomPart[which.min(distancesTemp),]-medialPointsEqualyDis_2D[1,]
  }
  
  # fill the gaps
  for (i in 1:nrow(medialPointsEqualyDis_2D)) {
    if(sum(is.na(upSpokesMeans[i,]))>0){
      upSpokesMeans[i,]<-upSpokesMeans[i-1,]+medialPointsEqualyDis_2D[i-1,]-medialPointsEqualyDis_2D[i,]
    }
  }
  for (i in 1:nrow(medialPointsEqualyDis_2D)) {
    if(sum(is.na(downSpokesMeans[i,]))>0){
      downSpokesMeans[i,]<-downSpokesMeans[i-1,]+medialPointsEqualyDis_2D[i-1,]-medialPointsEqualyDis_2D[i,]
    }
  }
  
  
  # # symmetric up and down
  # upSpokesMeans<-upSpokesMeans
  # for (i in 1:nrow(medialPointsEqualyDis_2D)) {
  #   if(sum(is.na(upSpokesMeans[i,]))>0 & sum(is.na(downSpokesMeans[i,]))==0){
  #     unitVecTemp<-convertVec2unitVec2(downSpokesMeans[i,])
  #     phi<-geodesicDistance(T_Vec[i,],unitVecTemp)
  #     r<-norm(downSpokesMeans[i,],type = '2')
  #     upSpokesMeans[i,]<-r*cos(phi)*T_Vec[i,]+r*sin(phi)*N_Vec[i,]
  #   }
  # }
  # downSpokesMeans<-downSpokesMeans
  # for (i in 1:nrow(medialPointsEqualyDis_2D)) {
  #   if(sum(is.na(downSpokesMeans[i,]))>0 & sum(is.na(upSpokesMeans[i,]))==0){
  #     unitVecTemp<-convertVec2unitVec2(upSpokesMeans[i,])
  #     phi<-geodesicDistance(T_Vec[i,],unitVecTemp)
  #     r<-norm(upSpokesMeans[i,],type = '2')
  #     downSpokesMeans[i,]<-r*cos(phi)*T_Vec[i,]-r*sin(phi)*N_Vec[i,]
  #   }
  # }
  
  if(plotting==TRUE){
    plot(rbind(boundaryPoints,boundaryPoints[1,]),type = "l", lty=3,xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    for (i in 1:nrow(medialPointsEqualyDis_2D)) {
      plot(rbind(upSpokesMeans[i,]+medialPointsEqualyDis_2D[i,],
                 medialPointsEqualyDis_2D[i,]),col="blue",type = 'l',xlim = xlim,ylim = xlim,xlab = "",ylab = "")
      par(new=TRUE)
    }
    plot(medialPointsEqualyDis_2D,col="black",pch=20,xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    for (i in 1:nrow(medialPointsEqualyDis_2D)) {
      plot(rbind(downSpokesMeans[i,]+medialPointsEqualyDis_2D[i,],
                 medialPointsEqualyDis_2D[i,]),col="red",type = 'l',xlim = xlim,ylim = xlim,xlab = "",ylab = "")
      par(new=TRUE)
    }
    par(new=TRUE)
    plot(rbind(extraSpokeSpineHead,
               extraSpokeSpineHead_Tip),col="black",type = 'l',xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(rbind(extraSpokeSpineTail,
               extraSpokeSpineTail_Tip),col="black",type = 'l',xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    
  }
  
  ######################################################################################################
  ######################################################################################################
  # cut and stretch spokes
  
  allSpokesTipsTemp<-rbind(upSpokesMeans+medialPointsEqualyDis_2D,
                           downSpokesMeans+medialPointsEqualyDis_2D,
                           extraSpokeSpineHead_Tip,
                           extraSpokeSpineTail_Tip)
  allSpokesTailsTemp<-rbind(medialPointsEqualyDis_2D,
                            medialPointsEqualyDis_2D,
                            medialPointsEqualyDis_2D[1,],
                            medialPointsEqualyDis_2D[nrow(medialPointsEqualyDis_2D),])
  
  cutAndStretchSpokes<-cutAndStretchSpokes_4Bivalves_2D(allSpokesTips = allSpokesTipsTemp,
                                                        allSpokesTails = allSpokesTailsTemp,
                                                        allMeshes2D)
  cutAndStretchSpokes<-cutAndStretchSpokes[,1:2]
  
  tipOfTopSpokes<-cutAndStretchSpokes[1:nrow(upSpokesMeans),]
  tipOfBottomSpokes<-cutAndStretchSpokes[(1+nrow(upSpokesMeans)):(nrow(upSpokesMeans)+nrow(downSpokesMeans)),]
  tipOfExtraSpokeHead<-cutAndStretchSpokes[(nrow(cutAndStretchSpokes)-1),]
  tipOfExtraSpokeTail<-cutAndStretchSpokes[nrow(cutAndStretchSpokes),]
  
  if(plotting==TRUE){
    plot(rbind(boundaryPoints,boundaryPoints[1,]),type = "l", lty=3,xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    for (i in 1:nrow(medialPointsEqualyDis_2D)) {
      plot(rbind(tipOfTopSpokes[i,],
                 medialPointsEqualyDis_2D[i,]),col="blue",type = 'l',xlim = xlim,ylim = xlim,xlab = "",ylab = "")
      par(new=TRUE)
    }
    plot(medialPointsEqualyDis_2D,col="black",pch=20,xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    for (i in 1:nrow(medialPointsEqualyDis_2D)) {
      plot(rbind(tipOfBottomSpokes[i,],
                 medialPointsEqualyDis_2D[i,]),col="red",type = 'l',xlim = xlim,ylim = xlim,xlab = "",ylab = "")
      par(new=TRUE)
    }
    par(new=TRUE)
    plot(rbind(tipOfExtraSpokeHead,
               medialPointsEqualyDis_2D[1,]),col="black",type = 'l',xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(rbind(tipOfExtraSpokeTail,
               medialPointsEqualyDis_2D[nrow(medialPointsEqualyDis_2D),]),col="black",type = 'l',xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    
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
  extraSpokeTailPoints<-generatePointsBetween2Points(medialPointsEqualyDis_2D[nrow(medialPointsEqualyDis_2D),],
                                                     tipOfExtraSpokeTail,
                                                     numberOfPoints = numberOf2DspokePoints)
  
  #plot 2D
  if(plotting==TRUE){
    plot(rbind(boundaryPoints,boundaryPoints[1,]),type = "l", lty=3,xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(topSpokesPoints,pch=20,cex=0.2,col="blue",xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(bottomSpokesPoints,pch=20,cex=0.2,col="red",xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(extraSpokeHeadPoints,pch=20,cex=0.2,col="black",xlim = xlim,ylim = xlim,xlab = "",ylab = "")
    par(new=TRUE)
    plot(extraSpokeTailPoints,pch=20,cex=0.2,col="black",xlim = xlim,ylim = xlim,xlab = "",ylab = "")
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
    plot(skeletalPoints_2D,type = "p",pch=20,col="blue",xlim = xlim,ylim = xlim,xlab = "",ylab = "")   
  }
  
  #transfer back
  # skeletalPoints_2D<-skeletalPoints_2D%*%t(rotationMatrix2D)+matrix(rep(centroid2D,dim(skeletalPoints_2D)[1]),ncol = 2,byrow = TRUE)
  
  #plot 2D
  if(plotting==TRUE){
    plot(skeletalPoints_2D,type = "p",pch=20,col="blue",xlim = xlim,ylim = xlim,xlab = "",ylab = "") 
    par(new=TRUE)
    for (i in 1:numberOfLayers) {
      plot(skeletalPoints_2D[matrixParentsChildren[i,],],type = "l",pch=20,col="blue",xlim = xlim,ylim = xlim,xlab = "",ylab = "") 
      par(new=TRUE)
    }  
  }
  
  crestMidlineIndices<-unique(c(matrixParentsChildren[numberOf2DspokePoints,],
                                matrixParentsChildren[2*numberOf2DspokePoints-1,]))
  
  if(plotting==TRUE){
    plot(skeletalPoints_2D[-crestMidlineIndices,],pch=20,col="blue",xlim = xlim,ylim = xlim,xlab = "",ylab = "") 
    par(new=TRUE)
    plot(skeletalPoints_2D[crestMidlineIndices,],pch=20,col="red",xlim = xlim,ylim = xlim,xlab = "",ylab = "") 
  }
  
  newData_2d<-data.frame(x=skeletalPoints_2D[,1], y=skeletalPoints_2D[,2])
  surfacePointsZ2<-predict(fit4, newdata = newData_2d)
  medialPoints3D<-cbind(skeletalPoints_2D[,1],skeletalPoints_2D[,2],surfacePointsZ2)
  
  numberOfFrames<-dim(medialPoints3D)[1]
  
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
  
  ######################################################################################################
  ######################################################################################################
  # 3D m-rep
  
  # increase the number of 3D boundary points by reduncing the
  # voxelSize of remeshing if it is necessary
  if(increase3DBoundaryPoint==TRUE){
    tmesh<-vcgUniformRemesh(tmesh,voxelSize = 0.5)
    if(plotting==TRUE){
      open3d()
      wire3d(tmesh, col="blue")
    }
  }
  
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
  k_closestPoint3D<-1
  averageOf_K_ClosestDistance3D<-rep(NA,numberOfFrames)
  for (i in 1:numberOfFrames) {
    tempArray<-sort(medialPoints3D_Dis2Boundary[i,])
    averageOf_K_ClosestDistance3D[i]<-mean(tempArray[1:k_closestPoint3D])
  }
  
  radiiAverage3D<-averageOf_K_ClosestDistance3D
  
  # #plot
  # if(plotting==TRUE){
  #   open3d()
  #   plot3d(medialPoints3D,type="s",radius = 0.2,col = "blue",expand = 10,box=FALSE,add = TRUE)
  #   for (i in 1:nrow(medialPoints3D)) {
  #     spheres3d(medialPoints3D[i,], radius = radiiAverage3D[i],col = "lightblue")
  #   }
  #   shade3d(tmesh, col="grey",alpha=0.4)
  # }
  
  allSpheresPoints<-c()
  allSpheres<-list()
  for (i in 1:nrow(medialPoints3D)) {
    tempSphere<-makeSphereMesh(center = medialPoints3D[i,],radius = radiiAverage3D[i],subdivision = 3)
    allSpheres[[i]]<-tempSphere
    allSpheresPoints<-rbind(allSpheresPoints,t(tempSphere$vb)[,1:3])
  }
  
  matrixInOut3D<-array(1,dim=c(nrow(allSpheresPoints),nrow(medialPoints3D)))
  for (j in 1:nrow(medialPoints3D)) {
    matrixInOut3D[,j]<-pip3d(Vertices = t(allSpheres[[j]]$vb[1:3,]),
                             Faces = t(allSpheres[[j]]$it),
                             Queries = allSpheresPoints)
  }
  labelsFinal3D<-rep(NA,nrow(allSpheresPoints))
  for (i in 1:nrow(matrixInOut3D)) {
    labelsFinal3D[i]<- 1 %in% matrixInOut3D[i,]
  }
  
  impliedBoundaryPoints<-allSpheresPoints[!labelsFinal3D,]
  
  if(plotting==TRUE){
    open3d()
    plot3d(impliedBoundaryPoints,type="p",col = "blue",expand = 10,box=FALSE,add = TRUE)
  }
  
  xData5<-data.frame(x=impliedBoundaryPoints[,1],y=impliedBoundaryPoints[,2])
  # xData5<-data.frame(x=t(tmesh$vb)[,1],y=t(tmesh$vb)[,2])
  surfacePointsZ5<-predict(fit4, newdata = xData5)
  
  # new3D_mesh<-cbind(t(tmesh$vb)[,1:2],surfacePointsZ5)
  new3D_mesh<-cbind(impliedBoundaryPoints[,1:2],surfacePointsZ5)
  
  #plot
  if(plotting==TRUE){
    open3d()
    shade3d(tmesh, col="white",alpha=0.2)
    plot3d(cbind(impliedBoundaryPoints[,1:2],surfacePointsZ5),type="p",
           col = "green",expand = 10,box=FALSE,add = TRUE)
  }
  
  # topPartMesh<-t(tmesh$vb)[,1:3][t(tmesh$vb)[,3]>surfacePointsZ5,]
  # bottomPartMesh<-t(tmesh$vb)[,1:3][t(tmesh$vb)[,3]<surfacePointsZ5,]
  
  topPartMesh<-impliedBoundaryPoints[impliedBoundaryPoints[,3]>surfacePointsZ5,]
  bottomPartMesh<-impliedBoundaryPoints[impliedBoundaryPoints[,3]<surfacePointsZ5,]
  
  
  #plot
  if(plotting==TRUE){
    open3d()
    plot3d(topPartMesh,type="p",radius = 0.2,col = "blue",expand = 10,box=FALSE,add = TRUE)
    plot3d(bottomPartMesh,type="p",radius = 0.2,col = "red",expand = 10,box=FALSE,add = TRUE)
  }
  
  # calculate 3D spokes 
  closestMedialLabelTop3D<-rep(NA,nrow(topPartMesh))
  pb <- txtProgressBar(style = 3)
  for (i in 1:nrow(topPartMesh)) {
    setTxtProgressBar(pb,i/nrow(topPartMesh))
    tempMatrix<-matrix(rep(topPartMesh[i,],nrow(medialPoints3D)),ncol = 3,byrow = TRUE)
    distancesTemp<-apply(X = medialPoints3D-tempMatrix, MARGIN = 1,FUN = myNorm)
    closestMedialLabelTop3D[i]<-which.min(distancesTemp)
  }
  close(pb)
  closestMedialLabelBottom3D<-rep(NA,nrow(bottomPartMesh))
  pb <- txtProgressBar(style = 3)
  for (i in 1:nrow(bottomPartMesh)) {
    setTxtProgressBar(pb,i/nrow(bottomPartMesh))
    tempMatrix<-matrix(rep(bottomPartMesh[i,],nrow(medialPoints3D)),ncol = 3,byrow = TRUE)
    distancesTemp<-apply(X = medialPoints3D-tempMatrix, MARGIN = 1,FUN = myNorm)
    closestMedialLabelBottom3D[i]<-which.min(distancesTemp)
  }
  
  # if(plotting==TRUE){
  #   open3d()
  #   for (i in 1:nrow(topPartMesh)) {
  #     plot3d(rbind(medialPoints3D[closestMedialLabelTop3D[i],],topPartMesh[i,]),type="l",col = "blue",expand = 10,box=FALSE,add = TRUE)
  #   }
  #   for (i in 1:nrow(bottomPartMesh)) {
  #     plot3d(rbind(medialPoints3D[closestMedialLabelBottom3D[i],],bottomPartMesh[i,]),type="l",col = "red",expand = 10,box=FALSE,add = TRUE)
  #   }
  #   shade3d(tmesh, col="white",alpha=0.2)  #surface mesh
  # }
  
  
  upSpokesMeans3D<-array(NA,dim = c(nrow(medialPoints3D),3))
  for (i in 1:nrow(medialPoints3D)) {
    if(sum(topPartMesh[which(closestMedialLabelTop3D==i),])==0){
      next
    }
    tempBoundary<-topPartMesh[which(closestMedialLabelTop3D==i),]
    tempSpokes<-tempBoundary-matrix(rep(medialPoints3D[i,],length(tempBoundary)/3),ncol = 3,byrow = TRUE)
    unitDirections<-apply(X = tempSpokes, MARGIN = 1, FUN = convertVec2unitVec2)
    frechetMeanTemp<-frechetMean(unitDirections)
    spokeLengths<-apply(X = tempSpokes, MARGIN = 1, FUN = myNorm)
    upSpokesMeans3D[i,]<-frechetMeanTemp*mean(spokeLengths)
  }
  downSpokesMeans3D<-array(NA,dim = c(nrow(medialPoints3D),3))
  for (i in 1:nrow(medialPoints3D)) {
    if(sum(bottomPartMesh[which(closestMedialLabelBottom3D==i),])==0){
      next
    }
    tempBoundary<-bottomPartMesh[which(closestMedialLabelBottom3D==i),]
    tempSpokes<-tempBoundary-matrix(rep(medialPoints3D[i,],length(tempBoundary)/3),ncol = 3,byrow = TRUE)
    unitDirections<-apply(X = tempSpokes, MARGIN = 1, FUN = convertVec2unitVec2)
    frechetMeanTemp<-frechetMean(unitDirections)
    spokeLengths<-apply(X = tempSpokes, MARGIN = 1, FUN = myNorm)
    downSpokesMeans3D[i,]<-frechetMeanTemp*mean(spokeLengths)
  }
  
  if(plotting==TRUE){
    open3d()
    for (i in 1:nrow(upSpokesMeans3D)) {
      plot3d(rbind(medialPoints3D[i,],medialPoints3D[i,]+upSpokesMeans3D[i,]),type="l",col = "blue",expand = 10,box=FALSE,add = TRUE)
    }
    for (i in 1:nrow(downSpokesMeans3D)) {
      plot3d(rbind(medialPoints3D[i,],medialPoints3D[i,]+downSpokesMeans3D[i,]),type="l",col = "red",expand = 10,box=FALSE,add = TRUE)
    }
    shade3d(tmesh, col="white",alpha=0.2)  #surface mesh
  }
  
  # cut and stretch spokes
  tipOfUpSpokes3D<-medialPoints3D+upSpokesMeans3D
  tipOfDownSpokes3D<-medialPoints3D+downSpokesMeans3D
  
  tempMesh4Up<-as.mesh3d(medialPoints3D)
  tempMesh4Up$normals<-rbind(apply(tipOfUpSpokes3D-medialPoints3D,FUN = convertVec2unitVec,MARGIN = 1),
                             rep(1,dim(medialPoints3D)[1]))
  tipOfCuttedUpSpokes<-vert2points(vcgRaySearch(tempMesh4Up,mesh = tmesh))
  
  tempMesh4Down<-as.mesh3d(medialPoints3D)
  tempMesh4Down$normals<-rbind(apply(tipOfDownSpokes3D-medialPoints3D,FUN = convertVec2unitVec,MARGIN = 1),
                               rep(1,dim(medialPoints3D)[1]))
  tipOfCuttedDownSpokes<-vert2points(vcgRaySearch(tempMesh4Down,mesh = tmesh))
  
  
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
  
  
  if(TRUE){
    # Plot LP-ds-rep
    open3d()
    # plot skeletal points
    plot3d(medialPoints3D,type="s",radius = 0.1,col = "blue",expand = 10,box=FALSE,add = TRUE)
    # plot radial flows
    for (i in 1:numberOfLayers) {
      plot3d(medialPoints3D[matrixParentsChildren[i,],],type="l",lwd = 2,col = "blue",expand = 10,box=FALSE,add = TRUE)
    }
    #plot mesh
    shade3d(tmesh, col="white",alpha=0.2) 
    
  }
  
  
  result<-list("medialPoints3D"=medialPoints3D,
               "tipOfCuttedUpSpokes"=tipOfCuttedUpSpokes,
               "tipOfCuttedDownSpokes"=tipOfCuttedDownSpokes)
  
  return(result)
  
}
