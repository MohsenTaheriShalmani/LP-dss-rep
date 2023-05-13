library(RiemBase)
library(pracma)

# Math functions

# function to calculate Mahalanobis distance (Hotelling Metric)
MahalanobisDistance<-function(X,Y){
  k<-dim(X)[2]
  nx<-dim(X)[1]
  ny<-dim(Y)[1]
  Sx<-cov(X)
  Sy<-cov(Y)
  meanX<-colMeans(X)
  meanY<-colMeans(Y)
  n<-nx+ny-1
  S<-((nx-1)*Sx+(ny-1)*Sy)/(nx+ny-2) #S=pooled covariance matrix
  # T2<-t(meanX-meanY)%*%solve(S*(1/nx+1/ny))%*%(meanX-meanY)
  T2<-t(meanX-meanY)%*%Ginv(S*(1/nx+1/ny))%*%(meanX-meanY) # Ginv is the Moore-Penrose generalized inverse
  return(as.double(T2))
}

MahalanobisDistance_1D<-function(X,Y){
  nx<-length(X)
  ny<-length(Y)
  Sx<-var(X)
  Sy<-var(Y)
  meanX<-mean(X)
  meanY<-mean(Y)
  n<-nx+ny-1
  S<-((nx-1)*Sx+(ny-1)*Sy)/(nx+ny-2) #S=pooled variance
  T2<-abs(meanX-meanY)*(1/(S*sqrt(1/nx+1/ny)))
  return(T2)
}

# Hotelling T2 test
HotellingT2<-function(X,Y){
  if(dim(X)[2]!=dim(Y)[2]){
    stop("Dimention Error!\n")
  }
  k<-dim(X)[2]
  nx<-dim(X)[1]
  ny<-dim(Y)[1]
  Sx<-cov(X)
  Sy<-cov(Y)
  meanX<-colMeans(X)
  meanY<-colMeans(Y)
  n<-nx+ny-1
  S<-((nx-1)*Sx+(ny-1)*Sy)/(nx+ny-2) #S=pooled covariance matrix
  # T2<-t(meanX-meanY)%*%solve(S*(1/nx+1/ny))%*%(meanX-meanY)
  T2<-t(meanX-meanY)%*%Ginv(S*(1/nx+1/ny))%*%(meanX-meanY) # Ginv is the Moore-Penrose generalized inverse
  F_value<-((n-k)/(k*(n-1)))*T2
  df1<-k
  df2<-n-k
  p_value<-1-pf(F_value,df1,df2)
  return(p_value)
}

permutation1D <- function(A,B,nPerm=10000) {
  if(!is.null(dim(A)) | !is.null(dim(B))){
    stop("Dimention of A & B must be equal 1!")
  }
  T0<-MahalanobisDistance_1D(A,B) #observed test statistics
  nSamplesG1<-length(A)
  nSamplesG2<-length(B)
  testStatistics<-rep(NA,nPerm)
  pooledGroup<-c(A,B)
  for (j in 1:nPerm) {
    permNumbers<-sample(1:(nSamplesG1+nSamplesG2))
    groupA_elementNumbers<-permNumbers[c(1:nSamplesG1)]
    groupB_elementNumbers<-permNumbers[c((nSamplesG1+1):(nSamplesG1+nSamplesG2))]
    
    tempX<-pooledGroup[groupA_elementNumbers]
    tempY<-pooledGroup[groupB_elementNumbers]
    
    testStatistics[j]<-MahalanobisDistance_1D(tempX,tempY)
  }
  pVal<-sum(testStatistics>T0)/nPerm
  return(pVal)
}

permutationMultivariate <- function(A,B,nPerm=10000) {
  if(dim(A)[2]!=dim(B)[2] | is.null(dim(A)) | is.null(dim(B))){
    stop("Dimention of A & B must be equal and greater than 1!")
  }
  T0<-MahalanobisDistance(A,B) #observed test statistics
  nSamplesG1<-dim(A)[1]
  nSamplesG2<-dim(B)[1]
  testStatistics<-rep(NA,nPerm)
  pooledGroup<-rbind(A,B)
  for (j in 1:nPerm) {
    permNumbers<-sample(1:(nSamplesG1+nSamplesG2))
    groupA_elementNumbers<-permNumbers[c(1:nSamplesG1)]
    groupB_elementNumbers<-permNumbers[c((nSamplesG1+1):(nSamplesG1+nSamplesG2))]
    
    tempX<-pooledGroup[groupA_elementNumbers,]
    tempY<-pooledGroup[groupB_elementNumbers,]
    
    testStatistics[j]<-MahalanobisDistance(tempX,tempY)
  }
  pVal<-sum(testStatistics>T0)/nPerm
  return(pVal)
}


meanDifferenceTest1D <- function(A,B,type="Parametric",nPerm=10000) {
  if(!is.null(dim(A)) | !is.null(dim(B))){
    stop("Dimention of A & B must be equal 1!")
  }
  if(type=="Parametric"){
    pVal<-t.test(A,B)$p.value
    return(pVal)
  }else if(type=="Permutation"){
    
    T0<-MahalanobisDistance_1D(A,B) #observed test statistics
    nSamplesG1<-length(A)
    nSamplesG2<-length(B)
    testStatistics<-rep(NA,nPerm)
    pooledGroup<-c(A,B)
    for (j in 1:nPerm) {
      
      permNumbers<-sample(1:(nSamplesG1+nSamplesG2))
      groupA_elementNumbers<-permNumbers[c(1:nSamplesG1)]
      groupB_elementNumbers<-permNumbers[c((nSamplesG1+1):(nSamplesG1+nSamplesG2))]
      
      tempX<-pooledGroup[groupA_elementNumbers]
      tempY<-pooledGroup[groupB_elementNumbers]
      
      testStatistics[j]<-MahalanobisDistance_1D(tempX,tempY)
    }
    pVal<-(0.5+sum(testStatistics>=T0))/(nPerm+1)
    return(pVal)
  }else{
    stop("Please verify the type as Parametric or Permutation!")
  }
}

meanDifferenceTestMultivariate <- function(A,B,type="Parametric",nPerm=10000) {
  if(dim(A)[2]!=dim(B)[2] | is.null(dim(A)) | is.null(dim(B))){
    stop("Dimention of A & B must be equal and greater than 1!")
  }
  if(type=="Parametric"){
    pVal<-HotellingT2(A,B)
    return(pVal)
  }else if(type=="Permutation"){
    
    T0<-MahalanobisDistance(A,B) #observed test statistics
    nSamplesG1<-dim(A)[1]
    nSamplesG2<-dim(B)[1]
    testStatistics<-rep(NA,nPerm)
    pooledGroup<-rbind(A,B)
    for (j in 1:nPerm) {
      
      permNumbers<-sample(1:(nSamplesG1+nSamplesG2))
      groupA_elementNumbers<-permNumbers[c(1:nSamplesG1)]
      groupB_elementNumbers<-permNumbers[c((nSamplesG1+1):(nSamplesG1+nSamplesG2))]
      
      tempX<-pooledGroup[groupA_elementNumbers,]
      tempY<-pooledGroup[groupB_elementNumbers,]
      
      testStatistics[j]<-MahalanobisDistance(tempX,tempY)
    }
    pVal<-(0.5+sum(testStatistics>=T0))/(nPerm+1)
    return(pVal)
  }else{
    stop("Please verify the type as Parametric or Permutation!")
  }
}



# library(BisRNA)
# The fisher.method function is also available in BisRNA library
fisher.method<-function (pvalues)
{
  df <- 2 * length(pvalues)
  global_pValue<-pchisq(-2 * sum(log(pvalues), na.rm = TRUE), df, lower.tail = FALSE)
  # global_pValue<-1-pchisq(-2 * sum(log(pvalues), na.rm = TRUE), df, lower.tail = TRUE)
  return(global_pValue)
}

# convert vectors to unit vectors
convertVec2unitVec <- function(vec) {
  if(norm(vec,type = "2")==0){
    stop("vector is zero!")
  }
  return(vec/norm(vec,type = "2"))
}

# convert vectors to unit vector or c(1,0,0)
convertVec2unitVec2 <- function(vec) {
  if(norm(vec,type = "2")==0){
    return(c(1,rep(0,length(vec)-1)))
  }
  return(vec/norm(vec,type = "2"))
}

# cross product of 2 vectors
myCrossProduct <- function(v,u) {
  return(c(v[2]*u[3]-v[3]*u[2],v[3]*u[1]-v[1]*u[3],v[1]*u[2]-v[2]*u[1]))
}



# frechet mean

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


clockwiseAngle <- function(v1,v2) {
  u1<-convertVec2unitVec(v1)
  u2<-convertVec2unitVec(v2)
  
  x1<-u1[1]
  y1<-u1[2]
  x2<-u2[1]
  y2<-u2[2]
  
  dot <- x1*x2 + y1*y2      # dot product between [x1, y1] and [x2, y2]
  det <- x1*y2 - y1*x2      # determinant
  angle <- atan2(det, dot)  # atan2(y, x) or atan2(sin, cos)
  
}


# calculate unit normal vector of triangle mesh
unitNormalOfTriangle <- function(point1,point2,point3) {
  a<-point2-point1
  b<-point3-point1
  
  normalVec<-c((a[2]*b[3]-a[3]*b[2]),-(a[1]*b[3]-a[3]*b[1]),(a[1]*b[2]-a[2]*b[1]))
  unitNormal<-normalVec/sqrt(sum(normalVec^2))
  
  return(unitNormal)
  
}


#this function is exactly like rgl::rotate3d(vec, angle, x, y, z)
rotationAboutXYZaxis <- function(vector,angle,axis=1) {
  if(axis==1){
    rotationMatrixTemp<-matrix(c(1,0,0,0,cos(angle),-sin(angle),0,sin(angle),cos(angle)),nrow = 3,byrow =T)
  }else if(axis==2){
    rotationMatrixTemp<-matrix(c(cos(angle),0,sin(angle),0,1,0,-sin(angle),0,cos(angle)),nrow = 3,byrow =T)  
  }else if(axis==3){
    rotationMatrixTemp<-matrix(c(cos(angle),-sin(angle),0,sin(angle),cos(angle),0,0,0,1),nrow = 3,byrow =T)  
  }else{
    stop("Axis is not acceptable! It must be 1 for X 2 for Y and 3 for Z")
  }
  return(as.vector(vector%*%rotationMatrixTemp))
}

# rotate3d(c(2,2,2), pi/4, 1, 0, 0)
# rotationAboutXYZaxis(c(2,2,2),angle = pi/4,axis = 1)


# generate random von Mises distribution on circle in radian
# converted code of Sungkyu Jung 2013, and Byungwon Kim 2017, MATLAB randvonMises.m
# mean is in[0,2pi] and kappa>0
randVonMises <- function(mean, kappa, n) {
  tau<-1+sqrt(1+4*kappa^2)
  rho<-(tau-sqrt(2*tau))/(2*kappa)
  r<-(1+rho^2)/(2*rho)
  
  u1<-runif(n)
  z<-cos(pi*u1)
  f<-(1+r*z)/(r+z)
  c<-kappa*(r-f)
  u2<-runif(n)
  acceptanceid<-(c*(2-c)-u2>0) | (log(c/u2)+1-c>=0)
  u3<-runif(sum(acceptanceid))
  theta<-sign(u3-0.5)*acos(f[acceptanceid])
  nnow<-length(theta)
  
  while (n>nnow) {
    n_more<-ceiling(n/nnow*(n-nnow))
    u1<-runif(n_more)
    z<-cos(pi*u1)
    f<-(1+r*z)/(r+z)
    c<-kappa*(r-f)
    u2<-runif(n_more)
    acceptanceid<-(c*(2-c)-u2>0) | (log(c/u2)+1-c>=0)
    u3<-runif(sum(acceptanceid))
    thetamore<-sign(u3-0.5)*acos(f[acceptanceid])
    
    theta<-c(theta, thetamore)
    nnow<-length(theta)
  }
  
  theta<-theta[1:n] + mean
  
  return(theta)
}

# randVonMisesSamples<-randVonMises(mean = pi/4,kappa = 10,n = 2000)
# hist(randVonMisesSamples)
# plotshapes(cbind(cos(randVonMisesSamples),sin(randVonMisesSamples)))


# converted code of Sungkyu Jung MATLAB randS2.m
# generate random sample of small sphere on S2 of the second kind 
# mu0, mu1 are directions and kappa0>1, kappa1>1
randS2 <- function(mu0,mu1,kappa0,kappa1,n) {
  
  mu0<-mu0/norm(mu0,type = "2")
  mu1<-mu1/norm(mu1,type = "2")
  nu<-sum(mu0*mu1)
  
  #generate Bingham-Mardia random vectors by the north pole
  
  x<-rnorm(n = n,mean = nu,sd = 1/sqrt(2*kappa0))
  x<-x[x<1 & x>-1]
  nnow<-length(x)
  while(n>nnow){
    n_more<-ceiling(n/nnow*(n-nnow))
    xmore<-rnorm(n = n_more,mean = nu,sd = 1/sqrt(2*kappa0))
    xmore<-xmore[xmore<1 & xmore>-1]
    x<-c(x,xmore)
    nnow<-length(x)
  }
  z<-x[1:n]
  
  #generate von Mises for longitude that c=mu0-nu*mu1 is parallel to x-axis
  theta<-randVonMises(mean = 0, kappa = kappa1, n = n)
  X_axis_northpole<-cbind(sqrt(1-z^2)*cos(theta),sqrt(1-z^2)*sin(theta), z)
  
  cx<-(mu1-nu*mu0)/sqrt(1-nu^2)
  cy<-cross(mu0,cx)
  cz<-mu0
  
  #rotate
  X<-X_axis_northpole%*%rbind(cx,cy,cz)
  
  return(X)
}

# draw circle on unit sphere S2 by center of small circle and r
# converted code of Sungkyu Jung MATLAB drawCircleS2.m
drawCircleS2 <- function(center,theta) {
  # NB!!! theta is the angle from center
  if(theta==pi/2){
    t<-cbind(cos(seq(0,2*pi,length.out = 50)),sin(seq(0,2*pi,length.out = 50)),rep(0,50))
    sCirc<-t%*%rotMat(center,c(0,0,1))
  }else{
    t<-cbind(sin(theta)*cos(seq(0,2*pi,length.out = 50)),sin(theta)*sin(seq(0,2*pi,length.out = 50)),cos(theta)*rep(1,50))
    sCirc<-t%*%rotMat(center,c(0,0,1))
  }
  spheres3d(x = 0, y = 0, z = 0, radius = 1,col = "lightblue", alpha=0.1)
  plot3d(sCirc,type="l",col = "black",expand = 10,box=TRUE,add = TRUE)
}


#generate ellipsoid PDM
ellipsoidGenerator_2D <- function(center,a,b,n,n2){
  phi<-seq(0, 2*pi, length.out = 4*n+1)
  points<-center
  for (i in phi[1:(length(phi)-1)]) {
    for (j in 1:n2) {
      x<-a*cos(i)/n2*j
      y<-b*sin(i)/n2*j
      points<-rbind(points,center+c(x,y))
    }
  }
  return(points)
}
# plot(ellipsoidGenerator_2D(center = c(2,3),a = 1,b = 1,10,1),xlim = c(0,5),ylim = c(0,5),xlab = "",ylab = "")

#generate ellipsoid PDM without center
ellipsoidGenerator_2D_2 <- function(center,a,b,n,n2){
  phi<-seq(0, 2*pi, length.out = 4*n+1)
  # points<-center
  points<-c()
  for (i in phi[1:(length(phi)-1)]) {
    for (j in 1:n2) {
      x<-a*cos(i)/n2*j
      y<-b*sin(i)/n2*j
      points<-rbind(points,center+c(x,y))
    }
  }
  return(points)
}

#generate ellipsoid PDM without center
sphereGenerator_2D <- function(center,r=1,n=16,asymmetric=TRUE){
  if(asymmetric==TRUE & (n %% 2) != 0){
    n<-n+1
  }
  phi<-seq(0, 2*pi, length.out = n)
  points<-array(NA,dim = c(n-1,2))
  for (i in 1:(n-1)) {
    x<-r*cos(phi[i])
    y<-r*sin(phi[i])
    points[i,]<-center+c(x,y)
  }
  return(points)
}

#generate a triangle
triangleGenerator <- function(bisectorVectorVertex,posVertex,sideLength,theta){
  u<-convertVec2unitVec(bisectorVectorVertex)
  v1<-c(cos(theta),sin(theta))
  v2<-c(cos(theta),-sin(theta))
  R<-rotMat(c(1,0),u)
  sidesEndPoints<-(rbind(v1,v2)%*%t(R))*(matrix(rep(1,4),nrow = 2)*sideLength)
  triangle<-rbind(posVertex,rbind(posVertex,posVertex)+sidesEndPoints)
  return(triangle)
}

#generate a pyramid
PyramidGenerator <- function(bisectorVectorVertex,posVertex,sideLength,theta){
  u<-convertVec2unitVec(bisectorVectorVertex)
  v1<-c(cos(theta),sin(theta),0)
  v2<-c(cos(theta),0,sin(theta))
  v3<-c(cos(theta),-sin(theta),0)
  v4<-c(cos(theta),0,-sin(theta))
  R<-rotMat(c(1,0,0),u)
  sidesEndPoints<-(rbind(v1,v2,v3,v4)%*%t(R))*(matrix(rep(1,12),ncol = 3)*sideLength)
  Pyramid<-rbind(posVertex,matrix(rep(posVertex,4),ncol = 3,byrow = TRUE)+sidesEndPoints)
  meshInfo<-rbind(c(1,2,3),c(1,3,4),c(1,4,5),c(1,5,2),c(2,3,4),c(4,5,2))
  out<-list(Pyramid=Pyramid,meshInfo=meshInfo)
  return(out)
}



# calculate triangle area
triangleArea <- function(p1,p2,p3) {
  myCrossProduct((p3-p2),(p2-p1))
}


#interpolate points between two points
generatePointsBetween2Points <- function(point1, point2, numberOfPoints) {
  
  dimention<-length(point1)
  
  totalDis<-norm(point1-point2,type = "2")
  smallDis<-seq(0, totalDis, length.out = numberOfPoints)
  
  direction<-convertVec2unitVec(point2-point1)
  
  middlePoints<-array(NA, dim = c(numberOfPoints,dimention))
  for (i in 1:length(smallDis)) {
    tempPoint<-point1+smallDis[i]*direction
    middlePoints[i,]<-tempPoint
  }
  
  return(middlePoints)
}
# testPoints<-generatePointsBetween2Points(point1 = c(1,1),point2 = c(5,2),numberOfPoints = 5)
# plot(testPoints,xlim = c(0,5),ylim = c(0,5),xlab = "",ylab = "")


#generate spheres between 2 spheres in 2D
generateSpheresBetween2Points_2D <- function(center1,
                                             center2,
                                             radius1,
                                             radius2,
                                             numberOfSpheres=10,
                                             resolution=20) {
  
  s1<-sphereGenerator_2D(center = center1,r = radius1,n = resolution,asymmetric = FALSE)
  s2<-sphereGenerator_2D(center = center2,r = radius2,n = resolution,asymmetric = FALSE)
  
  centers<-generatePointsBetween2Points(point1 = center1, point2 = center2,
                                        numberOfPoints = numberOfSpheres)
  
  radii<-seq(from=radius1,to=radius2,length.out = numberOfSpheres)
  
  numberOfBoundaryPoints<-nrow(s1)
  spheres<-array(NA,dim = c(numberOfSpheres*numberOfBoundaryPoints,2))
  
  for (i in 1:numberOfSpheres) {
    tempSphere<-sphereGenerator_2D(center = centers[i,],r = radii[i],n = resolution,asymmetric = FALSE)
    spheres[((i-1)*numberOfBoundaryPoints+1):(i*numberOfBoundaryPoints),]<-tempSphere
  }
  return(spheres)
}

#points whithout tips ans tails
generatePointsBetween2Points2 <- function(point1, point2, numberOfPoints) {
  
  dimention<-length(point1)
  
  totalDis<-norm(point1-point2,type = "2")
  smallDis<-seq(0, totalDis, length.out = numberOfPoints)
  
  direction<-convertVec2unitVec(point2-point1)
  
  middlePoints<-c()
  for (i in 2:(length(smallDis)-1)) {
    tempPoint<-point1+smallDis[i]*direction
    middlePoints<-rbind(middlePoints,tempPoint)
  }
  
  return(middlePoints)
}


#smooth the meshPDM
generatePointsBetween3Points <- function(point1,point2,point3,numberOf2DspokePoints) {
  
  centroid<-colMeans(rbind(point1,point2,point3))
  
  points<-c()
  
  points<-rbind(points,generatePointsBetween2Points2(centroid,point1,numberOfPoints = numberOf2DspokePoints))
  points<-rbind(points,generatePointsBetween2Points2(centroid,point2,numberOfPoints = numberOf2DspokePoints))
  points<-rbind(points,generatePointsBetween2Points2(centroid,point3,numberOfPoints = numberOf2DspokePoints))
  
  return(points)
}

#Euclidean norm
myNorm <- function(vec) {
  return(sqrt(sum(vec^2)))
}

#check a point is between 2 points on a line segment
isOnLineSegment <- function(point1,middlePoint,point2) {
  d1<-norm(point1-middlePoint,type = "2")
  d2<-norm(point2-middlePoint,type = "2")
  d<-norm(point1-point2,type = "2")
  if(abs(d1+d2-d)<1e-6){
    return(TRUE)
  }else{
    return(FALSE)
  }
}


#create sphere mesh with normals
makeSphereMesh <- function(center,radius,subdivision = 3) {
  tempSphere<-vcgSphere(subdivision = subdivision, normals = TRUE)
  tempSphere$vb[1:3,]<-radius*tempSphere$vb[1:3,]
  tempSphere$vb[1:3,]<-tempSphere$vb[1:3,]+
    t(matrix(rep(center,dim(t(tempSphere$vb)[,1:3])[1]),ncol = 3,byrow = TRUE))
  return(tempSphere)
  
}


# geodesicDistance between 2 unit vectors
geodesicDistance <- function(u1,u2) {
  acos(pmin(pmax(sum(u1*u2),-1.0),1.0))
}


#generate random points inside a triangle in 3D
randomPointOnTriangle <- function(n,point1,point2,point3) {
  v1<-point2-point1
  v2<-point3-point1
  triangleArea<-norm(myCrossProduct(v1,v2),type = "2")/2
  
  randPoints<-array(NA,dim=c(n,length(point1)))
  k<-1
  while (k<=n) {
    randPointTemp<-point1+runif(1)*v1+runif(1)*v2
    
    areaTriangle1<-norm(myCrossProduct((randPointTemp-point1),(randPointTemp-point2)),type = "2")/(2*triangleArea)
    areaTriangle2<-norm(myCrossProduct((randPointTemp-point1),(randPointTemp-point3)),type = "2")/(2*triangleArea)
    areaTriangle3<-norm(myCrossProduct((randPointTemp-point2),(randPointTemp-point3)),type = "2")/(2*triangleArea)
    
    if(areaTriangle1+areaTriangle2+areaTriangle3-1<1e-10){
      randPoints[k,]<-randPointTemp
      k<-k+1
    } 
  }
  return(randPoints)
}


geodesic2D_of2Vertices <- function(mesh2D,
                                   vertexIndex1,
                                   vertexIndex2,
                                   edgeLengths,
                                   total_Perimeter) {
  if(vertexIndex1==vertexIndex2){
    return(0)
  }else if(vertexIndex1>vertexIndex2){
    temp<-vertexIndex1
    vertexIndex1<-vertexIndex2
    vertexIndex2<-temp
  }
  
  geodesicDis1<-sum(edgeLengths[vertexIndex1:(vertexIndex2-1)])
  
  # geodesicDis1<-0
  # for (i in vertexIndex1:(vertexIndex2-1)) {
  #   geodesicDis1<-geodesicDis1+norm(mesh2D[i,]-mesh2D[i+1,],type = '2')
  # }
  
  geodesicDis2<-total_Perimeter-geodesicDis1
  
  geodesicDis<-min(geodesicDis1,geodesicDis2)
  
  return(geodesicDis)
  
}

curvatureOf3Points <- function(point1, point2, point3, mesh2D) {
  x1<-point1[1]
  y1<-point1[2]
  x2<-point2[1]
  y2<-point2[2]
  x3<-point3[1]
  y3<-point3[2]
  tempVal1<-2*abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1))
  tempVal2<-sqrt(((x2-x1)^2+(y2-y1)^2)*((x3-x1)^2+(y3-y1)^2)*((x3-x2)^2+(y3-y2)^2))
  if(tempVal1==0 | tempVal2==0){
    K<-0
  }else{
    K = tempVal1/tempVal2
    
    #The curvature is positive if the triangle centroid is inside the object 
    #otherwise it is negative
    triangleCentroid<-colMeans(rbind(point1,point2,point3))
    if(pip2d(mesh2D,rbind(triangleCentroid,triangleCentroid))[1]==-1){
      K<-(-K)
    }
  }
  
  return(K)
}



normalOfaVertex <- function(point1,vertex,point2) {
  
  vecTemp<-convertVec2unitVec2(point2-point1)
  
  R1<-rotMat(c(1,0),c(0,1))
  
  vertexUnitNormal<-vecTemp%*%t(R1)
  
  return(vertexUnitNormal)
  
}

normalOfaVertex2 <- function(point1,vertex,point2) {
  
  dx_line1 = vertex[1] - point1[1]
  dy_line1 = vertex[2] - point1[2]
  
  unitNormal_line1<-c(-dy_line1, dx_line1)/sqrt(sum(c(-dy_line1, dx_line1)^2))
  
  dx_line2 = point2[1] - vertex[1]
  dy_line2 = point2[2] - vertex[2]
  
  unitNormal_line2<-c(-dy_line2, dx_line2)/sqrt(sum(c(-dy_line2, dx_line2)^2))
  
  sumVec<-unitNormal_line1+unitNormal_line2
  vertexUnitNormal<-sumVec/norm(sumVec,type = '2')
  
  return(vertexUnitNormal)
  
}


verticesNormals_mesh2D <- function(mesh2D) {
  
  vetricesNormals<-array(NA,dim = c(nrow(mesh2D),2))
  tempMesh<-rbind(mesh2D[nrow(mesh2D),],mesh2D,mesh2D[1,])
  k<-1
  for (i in 2:(nrow(tempMesh)-1)) {
    vetricesNormals[k,]<-normalOfaVertex(tempMesh[i-1,],tempMesh[i,],tempMesh[i+1,])
    k<-k+1
  }
  
  if(1 %in% pip2d(mesh2D,mesh2D+vetricesNormals*0.1)){
    vetricesNormals<-(-vetricesNormals)
  }
  
}


linesegmentsIntersection <- function(P0,P1,P2,P3) {
  
  #line segment1 P0,P1
  #line segment1 P2,P3
  
  p0_x<-P0[1]
  p0_y<-P0[2]
  p1_x<-P1[1]
  p1_y<-P1[2] 
  p2_x<-P2[1]
  p2_y<-P2[2]
  p3_x<-P3[1]
  p3_y<-P3[2]
  
  s1_x = p1_x - p0_x
  s1_y = p1_y - p0_y
  s2_x = p3_x - p2_x
  s2_y = p3_y - p2_y
  
  s = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) / (-s2_x * s1_y + s1_x * s2_y)
  t = ( s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) / (-s2_x * s1_y + s1_x * s2_y)
  
  if (s >= 0 && s <= 1 && t >= 0 && t <= 1){
    x = p0_x + (t * s1_x)
    y = p0_y + (t * s1_y)
    return(c(x,y))
  }else{
    return(c(NA,NA))
  }
  
}

allAdjacentVertices2D <- function(mesh2D,k_Ring2D) {
  
  allAdjacentVertices<-list()
  refVec<-c(c(1:nrow(mesh2D)),c(1:nrow(mesh2D)),c(1:nrow(mesh2D)))
  for (i in 1:nrow(mesh2D)) {
    forwardIndices<-refVec[(nrow(mesh2D)+i+1):(nrow(mesh2D)+i+k_Ring2D)]
    backwardIndices<-refVec[(nrow(mesh2D)+i-1):(nrow(mesh2D)+i-k_Ring2D)]
    
    
    allAdjacentVertices[[i]]<-c(backwardIndices,forwardIndices)
  }
  
  return(allAdjacentVertices)
}

sortIndices <- function(indices) {
  group1<-indices[1]
  for (i in 2:length(indices)) {
    if((indices[i-1]+1)==indices[i]){
      group1<-c(group1,indices[i])
    }else{
      break
    }
  }
  group2<-indices[!indices %in% group1]
  if(length(group2)==1){
    newIndices<-group1
  }else if(group1[length(group1)]<group2[1] & length(group2)>0){
    newIndices<-c(group2,group1) 
  }else{
    newIndices<-indices
  }
  return(newIndices)
}


#cutted urchin
urchinInterpolation_2D <- function(numberOfMiddleSpokes,degree,coneVector, includeStartEnd=FALSE){
  
  if(degree<0 | degree>pi){
    stop("degree must be between 0 and pi")
  }
  
  phi<-seq(-degree/2, degree/2, length.out = numberOfMiddleSpokes+2) 
  
  if(includeStartEnd==FALSE){
    rangeTemp<-phi[2:(length(phi)-1)]
  }else{
    rangeTemp<-phi[1:length(phi)]
  }
  
  points<-c()
  for (i in rangeTemp) {
    x<-cos(i)
    y<-sin(i)
    points<-rbind(points,c(x,y))
  }
  eyePupil<-convertVec2unitVec(coneVector)
  
  R<-rotMat(c(1,0),eyePupil)
  eyeVectors<-points%*%t(R)
  return(eyeVectors)
}
# urchinCutted<-urchinInterpolation_2D(numberOfMiddleSpokes = 10,degree = pi/2,coneVector = c(1,1),includeStartEnd=FALSE)
# plot(urchinCutted,xlim = c(-1,1),ylim = c(-1,1),xlab = '',ylab = '')
# par(new=TRUE)
# plot(rbind(c(0,1),c(0,0),c(1,0)),type = 'l',xlim = c(-1,1),ylim = c(-1,1),xlab = '',ylab = '')



distanceOfaPoint2aLine_2D <- function(point,linePoint1,linePoint2) {
  
  x0<-point[1]
  y0<-point[2]
  x1<-linePoint1[1]
  y1<-linePoint1[2]
  x2<-linePoint2[1]
  y2<-linePoint2[2]
  
  result<-abs((x2-x1)*(y1-y0)-(x1-x0)*(y2-y1))/norm(linePoint1-linePoint2,type = '2')
  return(result)
  
}
# distanceOfaPoint2aLine_2D(point = c(1,0),linePoint1 = c(1,1),linePoint2 = c(2,1))


splineFitWithEquidistancePoints <- function(samplePoints,
                                            polyDegree2D,
                                            numberOfEquidistancePoints=200,
                                            lengthOut=10000) {
  # x y z
  x<-samplePoints[,1]
  y<-samplePoints[,2]
  
  fit4_2D <- lm(y ~ poly(x, degree = polyDegree2D,raw = TRUE),
                data=as.data.frame(cbind(x,y))) #NB!!! raw=F calculate orthogonal polynomial
  
  # fit4_2D <- ols(y ~ rcs(x,4))
  
  xData1<-data.frame(x=x)
  curvePoints1<-predict(fit4_2D, newdata = xData1)
  
  smoothedPoints<-cbind(x,curvePoints1)
  
  startPointX<-min(smoothedPoints[,1])
  endPointX<-max(smoothedPoints[,1])
  
  # tail to head
  allCloseX<-seq(startPointX,endPointX,length.out = lengthOut)
  
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
  
  
  equalLengths<-seq(0,totalLength,length.out = numberOfEquidistancePoints)
  
  selectedX<-rep(NA,numberOfEquidistancePoints)
  for (i in 1:numberOfEquidistancePoints) {
    tempIndex<-which.min((allDistances2StartingPoint-equalLengths[i])^2)
    selectedX[i]<-allCloseX[tempIndex]
  }
  
  xData2<-data.frame(x=selectedX)
  curvePoints2<-predict(fit4_2D, newdata = xData2)
  
  result<-cbind(selectedX,curvePoints2)
  
  return(result)
  
}
