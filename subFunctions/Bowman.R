
# up and down arrow generator
catEye <- function(numberOfArrows, linkVector, frame) {
  phi<-seq(0, pi/2, length.out = numberOfArrows)
  arrowsUp<-c()
  arrowsDown<-c()
  for (i in phi) {
    x<-cos(i)
    y<-0
    z<-sin(i)
    arrowsUp<-rbind(arrowsUp,c(x,y,z))
    arrowsDown<-rbind(arrowsDown,c(x,y,-z))
  }
  
  normalVec<- convertVec2unitVec(frame[1,])  #normal of a plane that target lives on it
  v<-convertVec2unitVec(linkVector)
  
  #projection of tempPoint on tangent space
  #sum(v*normalVec) is dot product that is the point distance from the plane
  projected_point <- convertVec2unitVec(v-(sum(v*normalVec)/norm(normalVec,type = "2"))*normalVec)
  
  initTarget<-c(1,0,0)
  R1<-rotMat(c(0,0,1),normalVec)
  targetUpdated1<-initTarget%*%t(R1)
  R2<-rotMat(targetUpdated1,projected_point)
  targetUpdated2<-targetUpdated1%*%t(R2)
  R3<-rotMat(projected_point,v)
  
  arrowsUpRotated<-arrowsUp%*%t(R3%*%R2%*%R1)
  arrowsDownRotated<-arrowsDown%*%t(R3%*%R2%*%R1)
  
  result<-list(arrowsUp=arrowsUpRotated, arrowsDown=arrowsDownRotated)
  return(result)
}


# frame<-diag(3)
# linkVector <- convertVec2unitVec(c(1,0,0))
# testOut<-catEye(numberOfArrows = 20,linkVector = linkVector,frame = frame)
# 
# open3d()
# spheres3d(x = 0, y = 0, z = 0, radius = 1,col = "lightblue", alpha=0.1)
# vectors3d(diag(3)*1.5,color="black", lwd=1)
# vectors3d(linkVector, color="black", lwd=1)
# vectors3d(frame[1,], color="blue", lwd=1)
# vectors3d(frame[2,], color="red", lwd=1)
# vectors3d(frame[3,], color="green", lwd=1)
# vectors3d(testOut$arrowsUp, color="purple", lwd=1)
# vectors3d(testOut$arrowsDown, color="pink", lwd=1)
# 
# 
# 
# frame<-framesInGlobalCoordinate[,,34]
# linkVector <- convertVec2unitVec(colMeans(frame))
# 
# phi<-seq(0, pi/2, length.out = 20)
# arrowsUp<-c()
# arrowsDown<-c()
# for (i in phi) {
#   x<-cos(i)
#   y<-0
#   z<-sin(i)
#   arrowsUp<-rbind(arrowsUp,c(x,y,z))
#   arrowsDown<-rbind(arrowsDown,c(x,y,-z))
# }

# normalVec<- convertVec2unitVec(frame[3,])  #normal of a plane that target lives on it
# linkVector<-c(1,0,0)
# v<-convertVec2unitVec(linkVector)
# 
# #projection of tempPoint on tangent space
# #sum(v*normalVec) is dot product that is the point distance from the plane
# projected_point <- convertVec2unitVec(v-(sum(v*normalVec)/norm(normalVec,type = "2"))*normalVec)
# 
# initTarget<-c(1,0,0)
# R1<-rotMat(c(0,0,1),normalVec)
# targetUpdated1<-initTarget%*%t(R1)
# R2<-rotMat(targetUpdated1,projected_point)
# targetUpdated2<-targetUpdated1%*%t(R2)
# R3<-rotMat(projected_point,v)
# 
# arrowsUpRotated<-arrowsUp%*%t(R3%*%R2%*%R1)
# arrowsDownRotated<-arrowsDown%*%t(R3%*%R2%*%R1)
# 
# open3d()
# spheres3d(x = 0, y = 0, z = 0, radius = 1,col = "lightblue", alpha=0.1)
# vectors3d(diag(3)*0.5,color="orange", lwd=1)
# vectors3d(linkVector*1.5, color="Blue", lwd=1)
# # vectors3d(frame[1,], color="black", lwd=1)
# # vectors3d(frame[2,], color="black", lwd=1)
# # vectors3d(frame[3,], color="black", lwd=1)
# vectors3d(arrowsUpRotated, color="purple", lwd=1)
# vectors3d(arrowsDownRotated, color="pink", lwd=1)
# 


eye <- function(numberOflayers,pointsInEachLayer,eyeOpenParameter,eyePupilVector){
  n<-trunc(pointsInEachLayer/4)
  n2<-numberOflayers
  phi<-seq(0, 2*pi, length.out = 4*n+1)
  a<-eyeOpenParameter 
  if(eyeOpenParameter<0 | eyeOpenParameter>pi){
    stop("eyeOpenParameter must be between 0 and pi")
  }
  points<-c(0,0)
  for (i in phi[1:(length(phi)-1)]) {
    for (j in 1:n2) {
      x<-a*cos(i)/n2*j
      y<-a*sin(i)/n2*j
      points<-rbind(points,c(x,y))
    }
  }
  eyeInit<-ExpNPd(t(points))
  eyePupil<-convertVec2unitVec(eyePupilVector)
  
  R<-rotMat(c(0,0,1),eyePupil)
  eyeVectors<-t(eyeInit)%*%t(R)
  return(eyeVectors)
}

eye1<-eye(numberOflayers = 3,pointsInEachLayer = 20,eyeOpenParameter = 0.5,eyePupilVector = c(1,1,1))

#plot three points of a frame
open3d()
spheres3d(x = 0, y = 0, z = 0, radius = 1,col = "lightblue", alpha=0.1)
vectors3d(convertVec2unitVec(c(1,1,1)), color="black", lwd=1,labels = "")
vectors3d(eye1, color="blue", lwd=1,labels = "")
vectors3d(1.5*diag(3), color="black", lwd=1,labels = "")


eyeIntersectionsWithBoundary <- function(meshPoints, meshPolygons, eyeOriginGlobal, eyeVectors) {
  
  numberOfPoints<-dim(meshPoints)[1]
  tipOfEyeVectors<-array(NA, dim = dim(eyeVectors))
  for (i in 1:dim(eyeVectors)[1]) {
    eyeVecTemp<-eyeVectors[i,]
    intersections<-array(NA,dim = c(dim(meshPolygons)[1],3))
    for (j in 1:dim(meshPolygons)[1]) {
      p1<-polyMatrix[j,1]
      p2<-polyMatrix[j,2]
      p3<-polyMatrix[j,3]
      point1<-meshPoints[p1,]
      point2<-meshPoints[p2,]
      point3<-meshPoints[p3,]
      tempIntersection<-rayTriangleIntersection(rayOrigin = eyeOriginGlobal,
                                                rayDirection = eyeVecTemp,
                                                triangleVertex1 = point1,
                                                triangleVertex2 = point2,
                                                triangleVertex3 = point3)
      intersections[j,]<-tempIntersection
    }
    distances<-rep(Inf,dim(meshPolygons)[1])
    for (k in 1:dim(meshPolygons)[1]) {
      if(!is.na(intersections[k,1])){
        distances[k]<-norm(intersections[k,]-eyeOriginGlobal,type = "2")
      }
    }
    
    tipOfEyeVectors[i,]<-intersections[which.min(distances),]
  }
  return(tipOfEyeVectors)
}


