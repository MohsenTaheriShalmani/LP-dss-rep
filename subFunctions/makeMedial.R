# functon to project a point on a line (to improve up and down spokes' lengths)
# a and b are on the line
projectPointOnLine<-function(a,b,p){
  ap = p-a
  ab = b-a
  projection = a + sum(ap*ab)/sum(ab*ab) * ab
  return(projection)
}



makeMedial <- function(centerPoint,upPoint,downPoint) {
  
  midPoint<-(upPoint+downPoint)/2
  
  secondPointOnTheLine<-centerPoint+convertVec2unitVec(downPoint-upPoint)
  
  updatedCenter<-projectPointOnLine(centerPoint,secondPointOnTheLine,midPoint)
  
  return(updatedCenter)
  
}
