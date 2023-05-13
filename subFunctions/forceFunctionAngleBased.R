forceFunctionAngleBased <- function(vectorsMatrix) {
  
  vectorsLengths<-apply(X = vectorsMatrix, MARGIN = 1,FUN = myNorm)
  # vectorsUnitDirections<-t(apply(X = vectorsMatrix, MARGIN = 1,FUN = convertVec2unitVec ))
  
  shortestVectors<-which.minn(vectorsLengths,n = 2)
  
  shortestLength<-vectorsLengths[shortestVectors[1]]
  
  v1<-vectorsMatrix[shortestVectors[1],]
  v2<-vectorsMatrix[shortestVectors[2],]
  
  u1<-convertVec2unitVec(v1)
  u2<-convertVec2unitVec(v2)
  
  angle<-acos(pmin(pmax(sum(u1*u2),-1.0), 1.0))
  
  force<-abs(pi/angle-1)
  
  v_mean<-colMeans(rbind(v1,v2))
  
  if(all(v_mean==0)){
    forceDirection<-rep(0,length(v_mean))
  }else{
    forceDirection<-convertVec2unitVec(v_mean)
  }
  
  out<-list(v1=v1,v2=v2,force=force,forceDirection=forceDirection,shortestLength=shortestLength)
  
  return(out)
  
}


forceFunctionAngleBased_4MultiObject <- function(vectorsMatrix,
                                                 labels,
                                                 type='one') {
  
  if(is.null(vectorsMatrix)){
    stop('vectorsMatrix is NULL!')
  }
  
  if(type=='one'){
    if(is.null(dim(vectorsMatrix))){
      out<-list(v1=vectorsMatrix,v2=vectorsMatrix,force=10^4,forceDirection=vectorsMatrix)
    }else{
      
      vectorsLengths<-apply(X = vectorsMatrix, MARGIN = 1,FUN = myNorm)
      # vectorsUnitDirections<-t(apply(X = vectorsMatrix, MARGIN = 1,FUN = convertVec2unitVec ))
      
      shortestVectors<-which.minn(vectorsLengths,n = 2)
      
      shortestLength<-vectorsLengths[shortestVectors[1]]
      
      v1<-vectorsMatrix[shortestVectors[1],]
      v2<-vectorsMatrix[shortestVectors[2],]
      
      u1<-convertVec2unitVec(v1)
      u2<-convertVec2unitVec(v2)
      
      if(labels[shortestVectors[1]]==labels[shortestVectors[2]]){
        
        force<-10^4
        
        v_mean<-colMeans(rbind(v1,v2))
        
        if(all(v_mean==0)){
          forceDirection<-rep(0,length(v_mean))
        }else{
          forceDirection<-convertVec2unitVec(v_mean)
        }
        
        out<-list(v1=v1,v2=v2,force=force,forceDirection=forceDirection,shortestLength=shortestLength)
        
      }else{
        angle<-acos(pmin(pmax(sum(u1*u2),-1.0), 1.0))
        
        force<-abs(pi/angle-1)
        
        v_mean<-colMeans(rbind(v1,v2))
        
        if(all(v_mean==0)){
          forceDirection<-rep(0,length(v_mean))
        }else{
          forceDirection<-convertVec2unitVec(v_mean)
        }
        
        out<-list(v1=v1,v2=v2,force=force,forceDirection=forceDirection,shortestLength=shortestLength)
      }
    }
  }else if(type=='two'){ # in type 2 we treat all the objects as one object
    
    if(is.null(dim(vectorsMatrix))){
      out<-list(v1=vectorsMatrix,v2=vectorsMatrix,force=10^4,forceDirection=vectorsMatrix)
    }else{
      
      vectorsLengths<-apply(X = vectorsMatrix, MARGIN = 1,FUN = myNorm)
      # vectorsUnitDirections<-t(apply(X = vectorsMatrix, MARGIN = 1,FUN = convertVec2unitVec ))
      
      shortestVectors<-which.minn(vectorsLengths,n = 2)
      
      v1<-vectorsMatrix[shortestVectors[1],]
      v2<-vectorsMatrix[shortestVectors[2],]
      
      u1<-convertVec2unitVec(v1)
      u2<-convertVec2unitVec(v2)
      
      angle<-acos(pmin(pmax(sum(u1*u2),-1.0), 1.0))
      
      force<-abs(pi/angle-1)
      
      v_mean<-colMeans(rbind(v1,v2))
      
      if(all(v_mean==0)){
        forceDirection<-rep(0,length(v_mean))
      }else{
        forceDirection<-convertVec2unitVec(v_mean)
      }
      
      out<-list(v1=v1,v2=v2,force=force,forceDirection=forceDirection,shortestLength=shortestLength)
    }
  }else{
    stop('Please specify the type as one or two!')
  }
  
  return(out)
  
}


urchin_AngleOfShortestSpokes <- function(vectorsMatrix) {
  
  vectorsLengths<-apply(X = vectorsMatrix, MARGIN = 1,FUN = myNorm)
  # vectorsUnitDirections<-t(apply(X = vectorsMatrix, MARGIN = 1,FUN = convertVec2unitVec ))
  
  shortestVectors<-which.minn(vectorsLengths,n = 2)
  
  v1<-vectorsMatrix[shortestVectors[1],]
  v2<-vectorsMatrix[shortestVectors[2],]
  
  u1<-convertVec2unitVec(v1)
  u2<-convertVec2unitVec(v2)
  
  angle<-acos(pmin(pmax(sum(u1*u2),-1.0), 1.0))
  
  return(angle)
  
}



urchin_AngleOfShortestSpokes_4MultiObject <- function(vectorsMatrix,
                                                      labels,
                                                      type='one') {
  
  if(is.null(vectorsMatrix)){
    stop('vectorsMatrix is NULL!')
  }
  
  if(type=='one'){
    if(is.null(dim(vectorsMatrix))){
      angle<-0
    }else{
      
      vectorsLengths<-apply(X = vectorsMatrix, MARGIN = 1,FUN = myNorm)
      # vectorsUnitDirections<-t(apply(X = vectorsMatrix, MARGIN = 1,FUN = convertVec2unitVec ))
      
      shortestVectors<-which.minn(vectorsLengths,n = 2)
      
      v1<-vectorsMatrix[shortestVectors[1],]
      v2<-vectorsMatrix[shortestVectors[2],]
      
      u1<-convertVec2unitVec(v1)
      u2<-convertVec2unitVec(v2)
      
      if(labels[shortestVectors[1]]==labels[shortestVectors[2]]){
        
       angle<-0
        
      }else{
        
        angle<-acos(pmin(pmax(sum(u1*u2),-1.0), 1.0))
        
      }
    }
  }else if(type=='two'){ # in type 2 we treat all the objects as one object
    
    if(is.null(dim(vectorsMatrix))){
      angle<-0
    }else{
      
      vectorsLengths<-apply(X = vectorsMatrix, MARGIN = 1,FUN = myNorm)
      # vectorsUnitDirections<-t(apply(X = vectorsMatrix, MARGIN = 1,FUN = convertVec2unitVec ))
      
      shortestVectors<-which.minn(vectorsLengths,n = 2)
      
      v1<-vectorsMatrix[shortestVectors[1],]
      v2<-vectorsMatrix[shortestVectors[2],]
      
      u1<-convertVec2unitVec(v1)
      u2<-convertVec2unitVec(v2)
      
      angle<-acos(pmin(pmax(sum(u1*u2),-1.0), 1.0))
      
    }
  }else{
    stop('Please specify the type as one or two!')
  }
  
  return(angle)
  
}
