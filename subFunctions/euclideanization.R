
euclideanization <- function(directions1, directions2, type="tangent space") {
  
  nSamplesG1<-dim(directions1)[2]
  nSamplesG2<-dim(directions2)[2]
  d<-dim(directions1)[1]
  
  pooledDirection<-cbind(directions1,directions2)
  
  # For extremely concentrated data we use Mardia mean direction
  pcaData1<-prcomp(t(pooledDirection))
  if(pcaData1$sdev[1]<1e-02 | pcaData1$sdev[2]<1e-02){
    mu_g<-convertVec2unitVec(colMeans(t(pooledDirection)))
    
    R <- rotMat(mu_g,c(rep(0,d-1),1))
    shiftedG1<-R%*%directions1
    shiftedG2<-R%*%directions2
    
    #log transfer to the tangent space
    logG1<-t(LogNPd(shiftedG1))
    logG2<-t(LogNPd(shiftedG2))
    
    result<-list(euclideanG1=logG1, euclideanG2=logG2)
    
  }else if(type=="tangent space"){
    
    # mean by Fre'chet mean
    allDirTemp<-t(pooledDirection)
    data1 <- list()
    for (j in 1:dim(allDirTemp)[1]){
      data1[[j]] <-allDirTemp[j,]
    }
    data2 <- riemfactory(data1, name="sphere")
    # Fre'chet Mean
    out1<- rbase.mean(data2)
    mu_g<-as.vector(out1$x)
    
    #rotate data to the north pole 
    R <- rotMat(mu_g,c(rep(0,d-1),1))
    shiftedG1<-R%*%directions1
    shiftedG2<-R%*%directions2
    
    #log transfer to the tangent space
    logG1<-t(LogNPd(shiftedG1))
    logG2<-t(LogNPd(shiftedG2))
    
    result<-list(euclideanG1=logG1, euclideanG2=logG2)
    
  }else if(type=="PNS"){
    
    # use pns instead of Fre'chet mean
    typeOfSphere<-kurtosisTestFunction(pooledDirection)
    pnsDirection<-pns(pooledDirection,sphere.type = typeOfSphere) #pooled directions
    res_G1<-t(pnsDirection$resmat[,1:nSamplesG1])
    res_G2<-t(pnsDirection$resmat[,(nSamplesG1+1):(nSamplesG1+nSamplesG2)])
    
    result<-list(euclideanG1=res_G1, euclideanG2=res_G2)
    
  }else{
    stop("Please choose the type of analysis i.e., 'PNS' or 'tangent space'!")
  }
  
  return(result)
}

