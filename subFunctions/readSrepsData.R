readSrepsData <- function(srepsData=srepsDataG1) {
  
  upSpoeksNumber<-max(srepsData$SpokesNumber[which(srepsData$srepNumber==1 & srepsData$Spoke=='up')])
  downSpoeksNumber<-max(srepsData$SpokesNumber[which(srepsData$srepNumber==1 & srepsData$Spoke=='down')])
  crestSpoeksNumber<-max(srepsData$SpokesNumber[which(srepsData$srepNumber==1 & srepsData$Spoke=='crest')])
  nTotalRadii <- upSpoeksNumber + downSpoeksNumber + crestSpoeksNumber
  skelPointNo <- nTotalRadii-downSpoeksNumber
  nSamples<-max(srepsData$srepNumber)
  skelRange<-c(1:downSpoeksNumber,(2*downSpoeksNumber+1):nTotalRadii)
  nSamples<-max(srepsData$srepNumber)
  
  SkeletalPDM<-array(NA,dim=c(nTotalRadii,3,nSamples))
  BoundaryPDM<-array(NA,dim=c(nTotalRadii,3,nSamples))
  boundaryPlusSkeletal<-array(NA,dim=c((2*nTotalRadii),3,nSamples))
  for (k in 1:nSamples){
    xPositions<-srepsData$xPos[which(srepsData$srepNumber==k )]
    yPositions<-srepsData$yPos[which(srepsData$srepNumber==k )]
    zPositions<-srepsData$zPos[which(srepsData$srepNumber==k )]
    Radii<-srepsData$radii[which(srepsData$srepNumber==k )]
    UxDirection<-srepsData$Ux[which(srepsData$srepNumber==k )]
    UyDirection<-srepsData$Uy[which(srepsData$srepNumber==k )]
    UzDirection<-srepsData$Uz[which(srepsData$srepNumber==k )]
    boundaryPointx<-xPositions+Radii*UxDirection
    boundaryPointy<-yPositions+Radii*UyDirection
    boundaryPointz<-zPositions+Radii*UzDirection
    
    skeletalPDMTemp<-cbind(xPositions,yPositions,zPositions)
    boundaryPDMTemp<-cbind(boundaryPointx,boundaryPointy,boundaryPointz)
    boundaryPlusSkeletalTemp<-rbind(skeletalPDMTemp,boundaryPDMTemp)
    
    SkeletalPDM[,,k]<-skeletalPDMTemp
    BoundaryPDM[,,k]<-boundaryPDMTemp
    boundaryPlusSkeletal[,,k]<-boundaryPlusSkeletalTemp
  }
  results<-list(SkeletalPDM=SkeletalPDM, BoundaryPDM=BoundaryPDM, boundaryPlusSkeletal=boundaryPlusSkeletal)
  
  return(results)
}
