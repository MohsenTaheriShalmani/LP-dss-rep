library(data.table)
source("subFunctions/readSrepsData.R")

fixFilipedSreps <- function(srepsData = srepsDataG1,
                            reference_u0 = reference_u0,
                            reference_u1 = reference_u1,
                            reference_u2 = reference_u2,
                            centerPoint = centerPoint, p2 = p2, p3 = p3,
                            gridSize="5_7") {
  #No. samples and spokes
  upSpoeksNumber<-max(srepsData$SpokesNumber[which(srepsData$srepNumber==1 & srepsData$Spoke=='up')])
  downSpoeksNumber<-max(srepsData$SpokesNumber[which(srepsData$srepNumber==1 & srepsData$Spoke=='down')])
  crestSpoksNumber<-max(srepsData$SpokesNumber[which(srepsData$srepNumber==1 & srepsData$Spoke=='crest')])
  nTotalRadii <- upSpoeksNumber + downSpoeksNumber + crestSpoksNumber
  nSamples<-max(srepsData$srepNumber)
  
  temp<-readSrepsData(srepsData = srepsData)
  SkeletalPDM<-temp$SkeletalPDM
  BoundaryPDM<-temp$BoundaryPDM
  
  if(gridSize=="5*7"){
    
    senario1_upSpokesNumber<-1:upSpoeksNumber
    senario1_downSpokesNumber<-senario1_upSpokesNumber
    senario1_crestSpokesNumber<-1:crestSpoksNumber
    senario1_Numbers<-c(senario1_upSpokesNumber,senario1_downSpokesNumber,senario1_crestSpokesNumber)
    senario2_upSpokesNumber<-c(1,2,3,4,50,51,7,48,49,10,46,47,13,44,45,16,42,43,19,40,41,22,38,39,25,36,
                               37,28,34,35,31,32,33,29,30,26,27,23,24,20,21,17,18,14,15,11,12,8,9,5,6)
    senario2_downSpokesNumber<-senario2_upSpokesNumber
    senario2_crestSpokesNumber<-c(1,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2)
    senario2_Numbers<-c(senario2_upSpokesNumber,senario2_downSpokesNumber,senario2_crestSpokesNumber)
    senario3_upSpokesNumber<-c(31,32,33,28,29,30,25,26,27,22,23,24,19,20,21,16,17,18,13,14,15,10,11,12,7,
                               8,9,4,5,6,1,2,3,50,51,48,49,46,47,44,45,42,43,40,41,38,39,36,37,34,35)
    senario3_downSpokesNumber<-senario3_upSpokesNumber
    senario3_crestSpokesNumber<-c(11,10,9,8,7,6,5,4,3,2,1,20,19,18,17,16,15,14,13,12)
    senario3_Numbers<-c(senario3_upSpokesNumber,senario3_downSpokesNumber,senario3_crestSpokesNumber)
    senario4_upSpokesNumber<-c(31,32,33,28,34,35,25,36,37,22,38,39,19,40,41,16,42,43,13,44,45,10,46,47,7,
                               48,49,4,50,51,1,2,3,5,6,8,9,11,12,14,15,17,18,20,21,23,24,26,27,29,30)
    senario4_downSpokesNumber<-senario4_upSpokesNumber
    senario4_crestSpokesNumber<-c(11,12,13,14,15,16,17,18,19,20,1,2,3,4,5,6,7,8,9,10)
    senario4_Numbers<-c(senario4_upSpokesNumber,senario4_downSpokesNumber,senario4_crestSpokesNumber)
    
    u0<-array(NA, dim = c(nSamples,3))
    u1<-array(NA, dim = c(nSamples,3))
    u2<-array(NA, dim = c(nSamples,3))
    for (k in 1:nSamples) {
      
      v0<-BoundaryPDM[centerPoint,,k]-SkeletalPDM[centerPoint,,k]
      u0[k,]<-v0/sqrt(sum(v0^2))
      v1<-SkeletalPDM[centerPoint,,k]-SkeletalPDM[p2,,k]
      u1[k,]<-v1/sqrt(sum(v1^2))
      v2<-SkeletalPDM[centerPoint,,k]-SkeletalPDM[p3,,k]
      u2[k,]<-v2/sqrt(sum(v2^2))
      
    }
    
    senarioTypes<-rep(NA,nSamples)
    upSpokeIsOK<-rep(NA,nSamples) 
    for (i in 1:nSamples) {
      
      geoDisU0<-acos(pmin(pmax(sum(reference_u0*u0[i,]),-1.0),1.0)) 
      geoDisU1<-acos(pmin(pmax(sum(reference_u1*u1[i,]),-1.0),1.0)) 
      geoDisU2<-acos(pmin(pmax(sum(reference_u2*u2[i,]),-1.0),1.0))
      
      if(geoDisU0<pi/2){
        upSpokeIsOK[i]<-TRUE
      }else{
        upSpokeIsOK[i]<-FALSE
      }
      
      if(geoDisU1<pi/2 & geoDisU2<pi/2){
        senario<-1
      }else if(geoDisU1<pi/2 & geoDisU2>pi/2){
        senario<-2
      }else if(geoDisU1>pi/2 & geoDisU2<pi/2){
        senario<-3
      }else{
        senario<-4
      }
      senarioTypes[i]<-senario
    }
  }else{
    stop("The initial grid size must be 5*7.")
  }
  
  senario1_Spoke<-c(rep("up",upSpoeksNumber),rep("down",downSpoeksNumber),rep("crest",crestSpoksNumber))
  senario2_Spoke<-c(rep("down",downSpoeksNumber),rep("up",upSpoeksNumber),rep("crest",crestSpoksNumber))
  
  newSpokeNumbers<-c()
  newSpoke<-c()
  for(i in 1:nSamples){
    if(senarioTypes[i]==1){
      newSpokeNumbers<-c(newSpokeNumbers,senario1_Numbers)
    }else if(senarioTypes[i]==2){
      newSpokeNumbers<-c(newSpokeNumbers,senario2_Numbers)
    }else if(senarioTypes[i]==3){
      newSpokeNumbers<-c(newSpokeNumbers,senario3_Numbers)
    }else if(senarioTypes[i]==4){
      newSpokeNumbers<-c(newSpokeNumbers,senario4_Numbers)
    }
    
    if(upSpokeIsOK[i]==TRUE){
      newSpoke<-c(newSpoke,senario1_Spoke) 
    }else{
      newSpoke<-c(newSpoke,senario2_Spoke) 
    }
  }
  # cat("Types of fliping:",senarioTypes,"\n")
  # cat("up down spokes was OK:",upSpokeIsOK,"\n")
  
  
  srepDataAdjusted<-copy(srepsData)
  srepDataAdjusted$SpokesNumber<-newSpokeNumbers
  srepDataAdjusted$Spoke<-newSpoke
  
  sorted_srepDataAdjusted<-setorder(srepDataAdjusted,srepNumber,
                                      -Spoke, #Note that minus is for the decending order i.e., up, down, crest
                                      SpokesNumber)

  # # or use the below code for sorting
  # sorted_srepDataAdjusted<-srepDataAdjusted[order(srepDataAdjusted$srepNumber,
  #                                                       srepDataAdjusted$Spoke,
  #                                                       srepDataAdjusted$SpokesNumber,
  #                                                     decreasing = c(FALSE,TRUE,FALSE)),]

  
  cat("Fliping issue is fixed. Please proceed!")
  
  return(sorted_srepDataAdjusted)
}

