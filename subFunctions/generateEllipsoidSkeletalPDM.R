# generate skeletal PDM
# converted c++ code of Zhiyuan Liu (vtkSlicerSkeletalRepresentationInitializerLogic.cxx)

generateEllipsoidSkeletalPDM <- function(rx, ry, nRows=5, nCols=7, ellipseScale=0.9) {

  mra <- rx * ellipseScale
  mrb <- ry * ellipseScale
  
  # compute the skeletal points
  nCrestPoints <- nRows*2 + (nCols-2)*2
  deltaTheta <- 2*pi/nCrestPoints
  skeletal_points_x<-array(NA,dim=c(nRows, nCols))
  skeletal_points_y<-array(NA,dim=c(nRows, nCols))
  r <- 0
  c <- 0
  for(i in 0:nCrestPoints){
    
    theta <- pi - deltaTheta * floor(nRows/2) - deltaTheta*i
    x <- mra * cos(theta)
    y <- mrb * sin(theta)
    
    # these crest points have no inward points (side or corner of the s-rep)
    skeletal_points_x[r+1, c+1] <- x
    skeletal_points_y[r+1, c+1] <- y
    
    if(i < nCols - 1){
      # top row of crest points
      c <- c + 1;
    }else if(i < nCols - 1 + nRows - 1){
      # right side col of crest points ( if the top-left point is the origin)
      r <- r + 1
    }else if(i < nCols - 1 + nRows - 1 + nCols - 1){
      # bottom row of crest points
      c <- c - 1
    }else{
      # left side col of crest points
      r <- r - 1
    }
    
    if((i < nCols - 1 & i > 0) | (i > nCols + nRows - 2 & i < 2*nCols + nRows - 3))
    {
      # compute skeletal points inward
      mx <- (mra * mra - mrb * mrb) * cos(theta) / mra # this is the middle line
      my <- 0
      dx <- x - mx
      dy <- y - my
      numSteps <- floor(nRows/2) # steps from crest point to the skeletal point
      stepSize <- 1 / numSteps # step size on the half side of srep
      for(j in 0:numSteps)
      {
        tempX <- mx + stepSize * j * dx
        tempY <- my + stepSize * j * dy
        if(i < nCols - 1)
        {
          # step from medial to top at current iteration on the top line
          currR <- numSteps - j
          skeletal_points_x[currR+1, c+1-1] <- tempX
          skeletal_points_y[currR+1, c+1-1] <- tempY
        }
        else
        {
          currR <- j + numSteps
          skeletal_points_x[currR+1, c+1+1] <- tempX
          skeletal_points_y[currR+1, c+1+1] <- tempY
        }
        
      }
      
    }
  }
  skelPDM<-cbind(matrix(skeletal_points_x),matrix(skeletal_points_y),rep(0,length(matrix(skeletal_points_x))))
  
  return(skelPDM)
}

# skelPDM<-generateEllipsoidSkeletalPDM(rx = 12,ry = 5,nRows = 7,nCols = 9,ellipseScale = 0.9)
# plot3d(skelPDM,type="s",radius = 0.1,col = "black",expand = 10,box=TRUE,add = TRUE)


generateEllipsoidSkeletalPDM2 <- function(rx, ry, nRows=5, nCols=7, ellipseScale=0.9) {
  
  mra <- rx * ellipseScale
  mrb <- ry * ellipseScale
  
  # compute the skeletal points
  nCrestPoints <- nRows*2 + (nCols-2)*2
  deltaTheta <- 2*pi/nCrestPoints
  numSteps <- floor(nRows/2) # steps from crest point to the skeletal point
  
  # Note that the middle line is a degenerated ellipse, on which there are less boundary point than outer ellipses
  skeletal_points_x<-array(NA,dim=c(nCrestPoints, numSteps+1))
  skeletal_points_y<-array(NA,dim=c(nCrestPoints, numSteps+1))
  r_resample <- 0
  
  for(i in 0:(nCrestPoints-1)){
    
    theta <- pi - deltaTheta*i
    x <- mra * cos(theta)
    y <- mrb * sin(theta)
    # compute skeletal points inward:
    mx <- (mra * mra - mrb * mrb) * cos(theta) / mra # this is the middle line
    my <- 0
    dx <- x - mx
    dy <- y - my
    
    stepSize <- 1.0 / numSteps # step size on the half side of srep
    
    for(j in 0:(numSteps)){
      
      skeletal_points_x[i+1, j+1] <- mx + stepSize * j * dx
      skeletal_points_y[i+1, j+1] <- my + stepSize * j * dy
      
    }
  }
  
  skelPDMTemp<-cbind(matrix(skeletal_points_x),matrix(skeletal_points_y),rep(0,length(matrix(skeletal_points_x))))
  skelPDM<-unique(skelPDMTemp[,1:3]) #remove duplicated points
  return(skelPDM)
}

# skelPDM<-generateEllipsoidSkeletalPDM2(rx = 12,ry = 5,nRows = 5,nCols = 7,ellipseScale = 0.9)
# plot3d(skelPDM,type="s",radius = 0.1,col = "black",expand = 10,box=TRUE,add = TRUE)

