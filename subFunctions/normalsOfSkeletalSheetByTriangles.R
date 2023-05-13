source("subFunctions/MathFunctions.R")

normalsOfSkeletalSheetByTriangles <- function(skeletalPDM) {
  
  skel_points<-skeletalPDM
  
  framesCenters <- c(16,13,10,7 ,4 ,1 ,2 ,3 ,19,22,25,28,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,29,30,26,27,23,24,20,21,17,18,14,15,11,12,8 ,9 ,5 ,6 ,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71)
  coreAndNeighboursList<-c(16,13,17,19,42,13,10,14,16,44,10,7 ,11,13,46,7 ,4 ,8 ,10,48,4 ,1 ,5 ,7 ,50,1 ,2 ,5 ,4 ,50,2 ,3 ,5 ,1 ,50,3 ,52,6 ,2 ,51,  19,22,40,16,20,  22,25,38,19,23,  25,28,36,22,26,  28,31,34,25,29,  31,32,34,28,29,  32,33,34,31,29,  33,62,35,32,30,  34,35,36,28,32,  35,63,37,34,33,  36,37,38,25,34,  37,64,39,36,35,  38,39,40,22,36,  39,65,41,38,37,  40,41,42,19,38,  41,66,43,40,39,  42,43,44,16,40,  43,67,45,42,41,  44,45,46,13,42,  45,68,47,44,43,  46,47,48,10,44,  47,69,49,46,45,  48,49,50,7,46,  49,70,51,48,47,  50,51,2,4,48,  51,71,3,50,49,  29,30,32,28,26,  30,61,33,29,27,  26,27,29,25,23,  27,60,30,26,24,  23,24,26,22,20,  24,59,27,23,21,  20,21,23,19,17,  21,58,24,20,18,  17,18,20,16,14,  18,57,21,17,15,  14,15,17,13,11,  15,56,18,14,12,  11,12,14,10,8,  12,55,15,11,9,  8 ,9,11,7,5,  9 ,54,12,8,6,  6 ,53,9,5,3,  5,6,8,4,2)
  crestFrames<-c(52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71)
  crestParents<-c(3,6,9,12,15,18,21,24,27,30,33,35,37,39,41,43,45,47,49,51)
  
  coreAndNeighbours<-matrix(coreAndNeighboursList,ncol = 5,byrow = T) #first column is core
  
  skelCorePointNo<-numberOfFrames-crestSpoksNumber
  normalsCore<-array(NA,dim=c(skelCorePointNo,3))
  for (i in 1:skelCorePointNo) {
    
    coreIndx<-coreAndNeighbours[i,1]
    neighbour1<-coreAndNeighbours[i,2]
    neighbour2<-coreAndNeighbours[i,3]
    neighbour3<-coreAndNeighbours[i,4]
    neighbour4<-coreAndNeighbours[i,5]
    
    core<-skel_points[coreIndx,]
    p1<-skel_points[neighbour1,]
    p2<-skel_points[neighbour2,]
    p3<-skel_points[neighbour3,]
    p4<-skel_points[neighbour4,]
    
    v1<-p1-core
    v2<-p2-core
    v3<-p3-core
    v4<-p4-core
    normal2<-convertVec2unitVec(myCrossProduct(v2,v3))
    normal3<-convertVec2unitVec(myCrossProduct(v3,v4))
    if(sum(v1==c(0,0,0))!=3){
      normal1<-convertVec2unitVec(myCrossProduct(v1,v2)) 
      normal4<-convertVec2unitVec(myCrossProduct(v4,v1))
      normalCore<-convertVec2unitVec(colMeans(rbind(normal1,normal2,normal3,normal4)))
    }else{
      normalCore<-convertVec2unitVec(colMeans(rbind(normal2,normal3)))
    }
    
    normalsCore[coreIndx,]<-(-normalCore)
  }
  
  normalsCrests<-array(NA,dim = c(crestSpoksNumber,3))
  for (i in 1:crestSpoksNumber) {
    normalsCrests[i,]<-normalsCore[crestParents[i],]
  }
  normals<-rbind(normalsCore,normalsCrests)
  
  return(normals)
  
}

# sampleNo<-1
# skel_points<-SkeletalPDMG1[skelRange,,sampleNo]
# normals<-normalsOfSkeletalPDMByTriangles(skel_points)
# open3d()
# for (i in 2:numberOfFrames) {
#   vectors3d(skel_points[framesCenters[i],],origin = skel_points[framesParents[i],],
#             headlength = 0.1,radius = 1/6, col="blue", lwd=1)
# }
# for (i in 1:numberOfFrames) {
#   vectors3d(skel_points[i,]+normals[i,],origin = skel_points[i,],headlength = 0.1,radius = 1/10, col="orange", lwd=1)
# }

