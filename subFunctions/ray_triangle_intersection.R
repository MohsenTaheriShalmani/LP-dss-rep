
# return the intersection of a ray with a triangle if the intersection exist
rayTriangleIntersection <- function(rayOrigin,rayDirection,triangleVertex1,triangleVertex2,triangleVertex3) {
  
  O<-rayOrigin #origin of the ray
  D<-rayDirection #direction of the ray 
  A<-triangleVertex1 #triangle vertices
  B<-triangleVertex2
  C<-triangleVertex3
  
  E1<-B-A
  E2<-C-A
  N<-myCrossProduct(E1,E2)
  
  det<-(-sum(D*N))
  invdet <- 1/det
  AO<-O-A
  DAO<-myCrossProduct(AO,D)
  u<-sum(E2*DAO)*invdet
  v<-(-sum(E1*DAO)*invdet)
  t<-sum(AO*N)*invdet
  if (abs(det) >= 1e-6 & t >= 0 & u >= 0 & v >= 0 & (u+v) <= 1){
    intersection<-O + t * D
  }else{
    intersection<-c(NA,NA,NA)
  }
  return(intersection)
}


rayTriangleIntersection2D <- function(rayOrigin,
                                      rayDirection,
                                      point1,
                                      point2) {
  v1 = rayOrigin - point1
  v2 = point2 - point1
  v3 = c(-rayDirection[2], rayDirection[1])
  
  
  dotProduct = sum(v2 * v3)
  if (abs(dotProduct) < 1e-6){
    return(c(NA,NA)) 
  }else{
    t1 = (v1[2]*v2[1]-v1[1]*v2[2]) / dotProduct
    
    t2 = (v1 * v3) / sum(v2*v3)
    
    if (t1 >= 0.0 && (t2 >= 0.0 && t2 <= 1.0)){
      return(rayOrigin+t1*rayDirection)
    }else{
      return(c(NA,NA)) 
    }
  }
}


