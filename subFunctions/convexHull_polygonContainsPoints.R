# Algorithm to find convex hull polygon of a set of points

orientation <- function(p,q,r){ #insert three points
  
  val <- (q[2] - p[2]) * (r[1] - q[1]) - (q[1] - p[1]) * (r[2] - q[2])
  
  if (val == 0){ 
    return(0)
  }else if(val < 0){ # clock or counterclock wise
    return(1) 
  }else{
    return(2)
  }
}


convexHullFunction <- function(points) {
  
  n<-dim(points)[1]
  
  hullPoints<-c()
  
  # The leftmost point
  l <-which.min(points[,1])
  
  p <- l

  repeat {
    # Add current point
    
    hullPoints<-rbind(hullPoints,points[p,])
    
    q <- (p+1)%%n+1
    
    for (i in 1:n){
      if (orientation(points[p,], points[i,], points[q,]) == 2){
        q <-i 
      }
    }
    
    p <- q
    
    if (p == l){break}
    
  }  
  
    return(hullPoints)  
}

# points<-rbind(c(0,3),c(2,2),c(1,1),c(2,1),c(3,0),c(0,0),c(3,3))
# plotshapes(points)
# convexHullPoints<-convexHullFunction(points)
# plotshapes(convexHullPoints)
