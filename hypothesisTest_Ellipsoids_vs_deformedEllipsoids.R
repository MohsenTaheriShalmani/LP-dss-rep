library(shapes)
library(rgl)
library(matlib)
library(RiemBase)
library(doBy) #for which.maxn and which.minn
library(alphahull) #for alpha convex hull
library(plotrix)
library(Rvcg) #for triangle mesh 
library(Matrix) #for forceSymmetric function
library(Morpho) 
library(alphashape3d) #for alpha 3d
library(ptinpoly) #to check whether a point is inside a 2D or 3D mesh or no
library(fields) #To plot color bar
library(ggplot2)
library(pracma) #for cross product
library(LearnGeom) #for project point to a line
library(rotations)
library(rms)  # for cubic spline


#####################################################################################################
#####################################################################################################

#clear the environment
remove(list=ls())

#####################################################################################################
#####################################################################################################
#set working directory to file location

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

#####################################################################################################
#####################################################################################################
# load functions
if(TRUE){
  source("subFunctions/MathFunctions.R")
  source("subFunctions/readSrepsData.R")
  source("subFunctions/euclideanization.R")
  source("subFunctions/ray_triangle_intersection.R")
  source("subFunctions/cutAndStrechSpokes.R")
  source("subFunctions/forceFunctionAngleBased.R")
  source("subFunctions/frameGenerator.R")
  source("subFunctions/normalsOfSkeletalSheetByTriangles.R")
  source("subFunctions/normalsOfSkeletalSheetBySpline.R")
  source("subFunctions/rotateFrameForwardAndBackward.R")
  source("subFunctions/meanFrames.R")
}

#####################################################################################################
#####################################################################################################
# choose the type of study

#choose the type of study "sizeAndShapeAnalysis" or "shapeAnalysis"
# typeOfStudy<-"shapeAnalysis"              #removing scale
typeOfStudy<-"sizeAndShapeAnalysis"     #preserving scale

# choose euclideanization method typeOfStudy4directions as "PNS" or "tangent space"
# PNS takes more time!
typeOfStudy4directions<-"tangent space"
# typeOfStudy4directions<-"PNS"


#choose type of mean direction
typeOfMeanDirection<-"Frechet"
# typeOfMeanDirection<-"PNS"


# choose typeOfTest as "Parametric" or "Permutation"
typeOfTest<-"Parametric"    # Fast Parametric test is based on normality assumption
# typeOfTest<-"Permutation" # nPerm default is 10000 permutations


#####################################################################################################
#####################################################################################################
# read s-rep data mesh data 

# read ellipsoid meshes

meshPoints <- read.csv(file = paste("files/Ellipsoid.csv",sep = ""),check.names = FALSE, header=TRUE, sep=",")

# connections of triangular mesh
PolygonsCsv <- read.csv("files/Mesh_Polygon.csv")
polyMatrix<-cbind((PolygonsCsv$point1+1),(PolygonsCsv$point2+1),(PolygonsCsv$point3+1))

# create mesh3d for the ellipsoid
spharmPDM<-matrix(meshPoints[[1]], ncol = 3, byrow = TRUE)
verts <- rbind(t(as.matrix(spharmPDM)),1)
trgls <- as.matrix(t(polyMatrix))
tmesh <- tmesh3d(verts, trgls)

#plot ellipsoid
if(TRUE){
  open3d()
  wire3d(tmesh, col="grey",alpha=0.2)  #wire mesh
  shade3d(tmesh, col="white",alpha=0.2)  #surface mech 
}


# read ds-rep data of the ellipsoid
srepsDataEllipsoid<- read.csv(file=paste("files/ellipsoid_skel.csv",sep = ""), header=TRUE, sep=",")

#No. samples and spokes
upSpoeksNumber<-max(srepsDataEllipsoid$SpokesNumber[which(srepsDataEllipsoid$srepNumber==1 & srepsDataEllipsoid$Spoke=='up')])
downSpoeksNumber<-max(srepsDataEllipsoid$SpokesNumber[which(srepsDataEllipsoid$srepNumber==1 & srepsDataEllipsoid$Spoke=='down')])
crestSpoksNumber<-max(srepsDataEllipsoid$SpokesNumber[which(srepsDataEllipsoid$srepNumber==1 & srepsDataEllipsoid$Spoke=='crest')])
nTotalRadii <- upSpoeksNumber + downSpoeksNumber + crestSpoksNumber
skelPointNo <- nTotalRadii-downSpoeksNumber
skelRange<-c(1:downSpoeksNumber,(2*downSpoeksNumber+1):nTotalRadii)

tempEllipsoid<-readSrepsData(srepsData = srepsDataEllipsoid)
SkeletalPDMEllipsoid<-tempEllipsoid$SkeletalPDM
BoundaryPDMEllipsoid<-tempEllipsoid$BoundaryPDM
boundaryPlusSkeletal_Ellipsoid<-tempEllipsoid$boundaryPlusSkeletal

# plot ds-rep
if(TRUE){
  sampleNo<-1
  plot3d(SkeletalPDMEllipsoid[skelRange,,sampleNo],type="s", radius = 0.4,col = "black",expand = 10,box=FALSE,add = TRUE)
  plot3d(BoundaryPDMEllipsoid[,,sampleNo],type="s", radius = 0.4,col = "black",expand = 10,box=FALSE,add = TRUE)
  
  
  srep1<-rbind(SkeletalPDMEllipsoid[,,sampleNo],BoundaryPDMEllipsoid[,,sampleNo])
  plot3d(SkeletalPDMEllipsoid[skelRange,,sampleNo],type="s", size=0.5,col = "blue",expand = 10,box=FALSE,add = TRUE)
  plot3d(BoundaryPDMEllipsoid[103:122,,sampleNo],type="s", radius = 0.2,col = "green",expand = 10,box=FALSE,add = TRUE)
  for (i in 1:nTotalRadii) {
    plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 0.5,col = "red",expand = 10,box=FALSE,add = TRUE)
  }
  plot3d(BoundaryPDMEllipsoid[,,sampleNo],type="s", radius = 0.2,col = "green",expand = 10,box=FALSE,add = TRUE)
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE, main = NULL, sub = NULL, top = T, aspect = FALSE, expand = 1.1)
}

#####################################################################################################
#####################################################################################################
# load simulated ellipsoids

load("files/simulatedEllipsoids_PDM.Rdata")
load("files/simulatedinflatedEllipsoids_PDM.Rdata")

#plot 10 samples of ellipsoids and inflated ellipsoids
if(TRUE){
  for (i in 1:10) {
    verts <- rbind(t(as.matrix(simulatedEllipsoids_PDM[,,i])),1)
    trgls <- as.matrix(t(polyMatrix))
    tmesh <- tmesh3d(verts, trgls)
    # wire3d(tmesh, col=sample(1:10,size = 1),alpha=1)  #wire mesh
    shade3d(tmesh, col="white",alpha=0.2)  #surface mech
  }
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE,
             # main = "Significant spokes' lengths",
             sub = NULL, top = T, aspect = FALSE, expand = 1.1)
  open3d()
  for (i in 1:10) {
    verts <- rbind(t(as.matrix(simulatedinflatedEllipsoids_PDM[,,i])),1)
    trgls <- as.matrix(t(polyMatrix))
    tmesh <- tmesh3d(verts, trgls)
    # wire3d(tmesh, col=sample(1:10,size = 1),alpha=1)  #wire mesh
    shade3d(tmesh, col="white",alpha=0.2)  #surface mech
  }
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE,
             # main = "Significant spokes' lengths",
             sub = NULL, top = T, aspect = FALSE, expand = 1.1)
}


#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################
############################ LP-ds-rep analysis #####################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################
# load ds-rep data of the ellipsoids and the inflated ellipsoids


load("files/skelPlusBoundary_simulatedEllipsoid.Rdata")
load("files/skelPlusBoundary_simulatedInflatedEllipsoid.Rdata")


SkeletalPDMG1<-skelPlusBoundary_simulatedEllipsoid[1:nTotalRadii,,]
BoundaryPDMG1<-skelPlusBoundary_simulatedEllipsoid[(nTotalRadii+1):dim(skelPlusBoundary_simulatedEllipsoid)[1],,]

SkeletalPDMG2<-skelPlusBoundary_simulatedInflatedEllipsoid[1:nTotalRadii,,]
BoundaryPDMG2<-skelPlusBoundary_simulatedInflatedEllipsoid[(nTotalRadii+1):dim(skelPlusBoundary_simulatedEllipsoid)[1],,]


nSamplesG1<-dim(skelPlusBoundary_simulatedEllipsoid)[3]
nSamplesG2<-dim(skelPlusBoundary_simulatedInflatedEllipsoid)[3]

#####################################################################################################
#####################################################################################################
# Find spokes' directions and lengths (radii) after the correspondance adjustment

radii_G1<-array(NA, dim=c(nTotalRadii,nSamplesG1))
for (k in 1:nSamplesG1) {
  for (i in 1:nTotalRadii) {
    radii_G1[i,k]<-norm(SkeletalPDMG1[i,,k]-
                          BoundaryPDMG1[i,,k],type = "2")
  }
}
radii_G2<-array(NA, dim=c(nTotalRadii,nSamplesG2))
for (k in 1:nSamplesG2) {
  for (i in 1:nTotalRadii) {
    radii_G2[i,k]<-norm(SkeletalPDMG2[i,,k]-
                          BoundaryPDMG2[i,,k],type = "2")
  }
}

spokeDirections_G1<-array(NA, dim=c(nTotalRadii,3,nSamplesG1))
for (k in 1:nSamplesG1) {
  for (i in 1:nTotalRadii) {
    spokeDirections_G1[i,,k]<-convertVec2unitVec(BoundaryPDMG1[i,,k]-SkeletalPDMG1[i,,k])
  }
}
spokeDirections_G2<-array(NA, dim=c(nTotalRadii,3,nSamplesG2))
for (k in 1:nSamplesG2) {
  for (i in 1:nTotalRadii) {
    spokeDirections_G2[i,,k]<-convertVec2unitVec(BoundaryPDMG2[i,,k]-SkeletalPDMG2[i,,k])
  }
}

#####################################################################################################
#####################################################################################################
# Plot a ds-rep and its mesh from group1

#plot an inflated ellipsoid with its ds-rep
sampleNo<-10 #choose sampleNo between 1 to nSamplesG2 to see other ds-reps
if(TRUE){
  open3d()
  srep1<-rbind(SkeletalPDMG2[,,sampleNo],BoundaryPDMG2[,,sampleNo])
  plot3d(SkeletalPDMG2[skelRange,,sampleNo],type="s", size=0.5,col = "blue",expand = 10,box=FALSE,add = TRUE)
  for (i in 1:nTotalRadii) {
    plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 2,col = "blue",expand = 10,box=FALSE,add = TRUE)
  }
  # plot mesh and normal vectors of a sample
  verts <- rbind(t(as.matrix(simulatedinflatedEllipsoids_PDM[,,sampleNo])),1)
  trgls <- as.matrix(t(polyMatrix))
  tmesh <- tmesh3d(verts, trgls)
  shade3d(tmesh, col="white",alpha=0.2)  #surface mech
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE, main = NULL, sub = NULL, top = T, aspect = FALSE, expand = 1.1)
}

# convert ds-reps to LP-ds-reps

# Define labels of the frames for the grid of the skeletal sheet
# NB frame 16 is its own parent
framesCenters   <-c(16,13,10,7 ,4 ,1 ,2 ,3 ,19,22,25,28,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,29,30,26,27,23,24,20,21,17,18,14,15,11,12,8 ,9 ,5 ,6 ,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71)
framesParents   <-c(16,16,13,10,7 ,4 ,1 ,2 ,16,19,22,25,28,31,32,28,34,25,36,22,38,19,40,16,42,13,44,10,46,7 ,48,4 ,50,28,29,25,26,22,23,19,20,16,17,13,14,10,11,7 ,8 ,4 ,5 ,3 ,6 ,9 ,12,15,18,21,24,27,30,33,35,37,39,41,43,45,47,49,51)
framesBackPoints<-c(13,16,13,10,7 ,4 ,1 ,2 ,16,19,22,25,28,31,32,28,34,25,36,22,38,19,40,16,42,13,44,10,46,7 ,48,4 ,50,28,29,25,26,22,23,19,20,16,17,13,14,10,11,7 ,8 ,4 ,5 ,3 ,6 ,9 ,12,15,18,21,24,27,30,33,35,37,39,41,43,45,47,49,51)
framesFronts    <-c(19,10,7 ,4 ,1 ,2 ,3 ,52,22,25,28,31,32,33,62,35,63,37,64,39,65,41,66,43,67,45,68,47,69,49,70,51,71,30,61,27,60,24,59,21,58,18,57,15,56,12,55,9 ,54,6 ,53,rep(Inf,20)) #NB! crest frames don't have front point
# number of frames
numberOfFrames<-length(framesCenters)

# Calculate normal vectors of the skeletal sheet
# A general solution to calculate normals could be fitting a spline surface to the skeletal points.
# But spline fitting needs pre-alignment of skeletals to XY plane. Also it add variation to the normals.
# Here we generate normals based on the quadrilateral structure of the skeletal sheet.
#choose a method to find normals "sheetTriangles" or "splineFitting"
method2findNormals<-"sheetTriangles"
# method2findNormals<-"splineFitting"

if(method2findNormals=="sheetTriangles"){
  
  skeletalSheet_G1<-SkeletalPDMG1[skelRange,,]
  skeletalSheet_G2<-SkeletalPDMG2[skelRange,,]
  medialNormals_G1<-array(NA,dim = dim(skeletalSheet_G1))
  pb <- txtProgressBar(style = 3)
  for (i in 1:nSamplesG1) {
    setTxtProgressBar(pb, i/nSamplesG1)
    medialNormals_G1[,,i]<-normalsOfSkeletalSheetByTriangles(skeletalPDM = skeletalSheet_G1[,,i])
  }
  print("Group 1 is done!")
  medialNormals_G2<-array(NA,dim = dim(skeletalSheet_G2))
  pb <- txtProgressBar(style = 3)
  for (i in 1:nSamplesG2) {
    setTxtProgressBar(pb, i/nSamplesG2)
    medialNormals_G2[,,i]<-normalsOfSkeletalSheetByTriangles(skeletalPDM = skeletalSheet_G2[,,i])
  }
  print("Group 2 is done!")
  
}else if(method2findNormals<-"splineFitting"){
  
  # Calculate normals by spline fitting. Takes a minute.
  # To fit spline, we use PCA + GPA to have skeletal sheets along the XY axis
  # Note that the GPA is not necessary, but we use it to make sure all the 
  # normals correspond to the northern of the ellipsoid
  sampleNo<-1 
  srep1<-rbind(SkeletalPDMG1[,,sampleNo],BoundaryPDMG1[,,sampleNo])
  centeredRefSrep<-prcomp(srep1,center = T)$x
  alignedByRefSrep_G1<-array(NA,dim = dim(boundaryPlusSkeletal_G1))
  for (i in 1:nSamplesG1) {
    alignedByRefSrep_G1[,,i]<-procOPA(centeredRefSrep,boundaryPlusSkeletal_G1[,,i],scale = F)$Bhat
  }
  alignedByRefSrep_G2<-array(NA,dim = dim(boundaryPlusSkeletal_G2))
  for (i in 1:nSamplesG2) {
    alignedByRefSrep_G2[,,i]<-procOPA(centeredRefSrep,boundaryPlusSkeletal_G2[,,i],scale = F)$Bhat
  }
  skeletalSheet_G1<-alignedByRefSrep_G1[skelRange,,]
  skeletalSheet_G2<-alignedByRefSrep_G2[skelRange,,]
  
  # Update skeletal PDMs
  SkeletalPDMG1<-alignedByRefSrep_G1[1:nTotalRadii,,]
  SkeletalPDMG2<-alignedByRefSrep_G2[1:nTotalRadii,,]
  
  # Update spokes directions
  # Since we applied PCA, we need to update spokes' directions in global coordinate system
  spokeDirections_G1<-array(NA, dim=c(nTotalRadii,3,nSamplesG1))
  for (k in 1:nSamplesG1) {
    for (i in 1:nTotalRadii) {
      spokeDirections_G1[i,,k]<-convertVec2unitVec(alignedByRefSrep_G1[i+nTotalRadii,,k]-alignedByRefSrep_G1[i,,k])
    }
  }
  spokeDirections_G2<-array(NA, dim=c(nTotalRadii,3,nSamplesG2))
  for (k in 1:nSamplesG2) {
    for (i in 1:nTotalRadii) {
      spokeDirections_G2[i,,k]<-convertVec2unitVec(alignedByRefSrep_G2[i+nTotalRadii,,k]-alignedByRefSrep_G2[i,,k])
    }
  }
  
  #calculate normals
  medialNormals_G1<-array(NA,dim = dim(skeletalSheet_G1))
  pb <- txtProgressBar(style = 3)
  for (i in 1:nSamplesG1) {
    setTxtProgressBar(pb, i/nSamplesG1)
    temp<-normalsOfSkeletalSheetBySpline(centeredSkel = skeletalSheet_G1[,,i])
    medialNormals_G1[,,i]<-temp$medialNormals
  }
  close(pb)
  print("Group 1 is done!")
  medialNormals_G2<-array(NA,dim = dim(skeletalSheet_G2))
  pb <- txtProgressBar(style = 3)
  for (i in 1:nSamplesG2) {
    setTxtProgressBar(pb, i/nSamplesG2)
    temp<-normalsOfSkeletalSheetBySpline(centeredSkel = skeletalSheet_G2[,,i])
    medialNormals_G2[,,i]<-temp$medialNormals
  }
  close(pb)
  print("Group 2 is done!")
  
}else{
  stop("Please specify method2findNormals!")
}


#####################################################################################################
#####################################################################################################
# Calculate frames vectors 

# frames in global coordinate system
if(TRUE){
  framesFirstVectors_G1<-array(NA,dim = c(numberOfFrames,3,nSamplesG1))
  framesSecondVectors_G1<-array(NA,dim = c(numberOfFrames,3,nSamplesG1))
  framesThirdVectors_G1<-array(NA,dim = c(numberOfFrames,3,nSamplesG1))
  pb <- txtProgressBar(style = 3)
  for (i in 1:nSamplesG1) {
    setTxtProgressBar(pb, i/nSamplesG1)
    temp<-frameGenerator(centeredSkel = skeletalSheet_G1[,,i],medialNormals = medialNormals_G1[,,i],
                         framesCenters = framesCenters,framesBackPoints = framesBackPoints,framesFronts = framesFronts)
    
    framesFirstVectors_G1[,,i]<-temp$framesFirstVec
    framesSecondVectors_G1[,,i]<-temp$framesSecondVec
    framesThirdVectors_G1[,,i]<-temp$framesThirdVec
  }
  close(pb)
  print("Group 1 is done!")
  framesFirstVectors_G2<-array(NA,dim = c(numberOfFrames,3,nSamplesG2))
  framesSecondVectors_G2<-array(NA,dim = c(numberOfFrames,3,nSamplesG2))
  framesThirdVectors_G2<-array(NA,dim = c(numberOfFrames,3,nSamplesG2))
  pb <- txtProgressBar(style = 3)
  for (i in 1:nSamplesG2) {
    setTxtProgressBar(pb, i/nSamplesG2)
    temp<-frameGenerator(centeredSkel = skeletalSheet_G2[,,i],medialNormals = medialNormals_G2[,,i],
                         framesCenters = framesCenters,framesBackPoints = framesBackPoints,framesFronts = framesFronts)
    
    framesFirstVectors_G2[,,i]<-temp$framesFirstVec
    framesSecondVectors_G2[,,i]<-temp$framesSecondVec
    framesThirdVectors_G2[,,i]<-temp$framesThirdVec
  }
  close(pb)
  print("Group 2 is done!")
}

# Plot LP-ds-rep
sampleNo<-1 #choose sampleNo between 1 to nSamplesG1 to see other ds-reps
#plot
if(TRUE){
  open3d()
  for (i in 2:numberOfFrames) {
    vectors3d(skeletalSheet_G1[framesCenters[i],,sampleNo],origin = skeletalSheet_G1[framesParents[i],,sampleNo],
              headlength = 0.1,radius = 1/6, col="blue", lwd=1)
  }
  # for (i in framesCenters) {
  #   vectors3d(skeletalSheet_G1[i,,sampleNo]+framesFirstVectors_G1[i,,sampleNo],origin = skeletalSheet_G1[i,,sampleNo],headlength = 0.1,radius = 1/10, col="orange", lwd=1)
  #   vectors3d(skeletalSheet_G1[i,,sampleNo]+framesSecondVectors_G1[i,,sampleNo],origin = skeletalSheet_G1[i,,sampleNo],headlength = 0.1,radius = 1/10, col="orange", lwd=1)
  #   vectors3d(skeletalSheet_G1[i,,sampleNo]+framesThirdVectors_G1[i,,sampleNo],origin = skeletalSheet_G1[i,,sampleNo],headlength = 0.1,radius = 1/10, col="orange", lwd=1)
  # }
  # for (i in 1:nTotalRadii) {
  #   vectors3d(SkeletalPDMG1[i,,sampleNo]+(spokeDirections_G1[i,,sampleNo]*radii_G1[i,sampleNo]),origin = SkeletalPDMG1[i,,sampleNo],headlength = 0.1,radius = 1/10, col="red", lwd=1)
  # }
  for (i in 1:nTotalRadii) {
    plot3d(rbind(SkeletalPDMG1[i,,sampleNo]+(spokeDirections_G1[i,,sampleNo]*radii_G1[i,sampleNo]),
                 SkeletalPDMG1[i,,sampleNo]),type="l",lwd = 2,col = "grey",expand = 10,box=FALSE,add = TRUE)
  }
  verts <- rbind(t(as.matrix(simulatedEllipsoids_PDM[,,sampleNo])),1)
  trgls <- as.matrix(t(polyMatrix))
  tmesh <- tmesh3d(verts, trgls)
  # wire3d(tmesh, col=sample(1:10,size = 1),alpha=1)  #wire mesh
  shade3d(tmesh, col="white",alpha=0.2)  #surface mech
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE, main = NULL, sub = NULL, top = T, aspect = FALSE, expand = 1.1)
  
  open3d()
  for (i in 2:numberOfFrames) {
    vectors3d(skeletalSheet_G2[framesCenters[i],,sampleNo],origin = skeletalSheet_G2[framesParents[i],,sampleNo],
              headlength = 0.1,radius = 1/6, col="blue", lwd=1)
  }
  # for (i in framesCenters) {
  #   vectors3d(skeletalSheet_G2[i,,sampleNo]+framesFirstVectors_G2[i,,sampleNo],origin = skeletalSheet_G2[i,,sampleNo],headlength = 0.1,radius = 1/10, col="orange", lwd=1)
  #   vectors3d(skeletalSheet_G2[i,,sampleNo]+framesSecondVectors_G2[i,,sampleNo],origin = skeletalSheet_G2[i,,sampleNo],headlength = 0.1,radius = 1/10, col="orange", lwd=1)
  #   vectors3d(skeletalSheet_G2[i,,sampleNo]+framesThirdVectors_G2[i,,sampleNo],origin = skeletalSheet_G2[i,,sampleNo],headlength = 0.1,radius = 1/10, col="orange", lwd=1)
  # }
  # for (i in 1:nTotalRadii) {
  #   vectors3d(SkeletalPDMG2[i,,sampleNo]+(spokeDirections_G2[i,,sampleNo]*radii_G2[i,sampleNo]),origin = SkeletalPDMG2[i,,sampleNo],headlength = 0.1,radius = 1/10, col="red", lwd=1)
  # }
  for (i in 1:nTotalRadii) {
    plot3d(rbind(SkeletalPDMG2[i,,sampleNo]+(spokeDirections_G2[i,,sampleNo]*radii_G2[i,sampleNo]),
                 SkeletalPDMG2[i,,sampleNo]),type="l",lwd = 2,col = "grey",expand = 10,box=FALSE,add = TRUE)
  }
  verts <- rbind(t(as.matrix(simulatedinflatedEllipsoids_PDM[,,sampleNo])),1)
  trgls <- as.matrix(t(polyMatrix))
  tmesh <- tmesh3d(verts, trgls)
  # wire3d(tmesh, col=sample(1:10,size = 1),alpha=1)  #wire mesh
  shade3d(tmesh, col="white",alpha=0.2)  #surface mech
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE, main = NULL, sub = NULL, top = T, aspect = FALSE, expand = 1.1)
  
}

# Combine frames vectors to make SO(3) frames in global coordinate system
if(TRUE){
  frames_G1<-array(NA,dim = c(3,3,numberOfFrames,nSamplesG1))
  for (k in 1:nSamplesG1) {
    for (i in framesCenters) {
      frames_G1[,,i,k]<- rbind(framesFirstVectors_G1[i,,k],
                               framesSecondVectors_G1[i,,k],
                               framesThirdVectors_G1[i,,k])
    }
  }
  frames_G2<-array(NA,dim = c(3,3,numberOfFrames,nSamplesG2))
  for (k in 1:nSamplesG2) {
    for (i in framesCenters) {
      frames_G2[,,i,k]<- rbind(framesFirstVectors_G2[i,,k],
                               framesSecondVectors_G2[i,,k],
                               framesThirdVectors_G2[i,,k])
    }
  }
}

# Calculate children frames coordinates based on their parents frames
if(TRUE){
  framesBasedOnParents_G1<-array(NA,dim = c(3,3,numberOfFrames,nSamplesG1))
  for (k in 1:nSamplesG1) {
    for (i in 1:numberOfFrames) {
      k1<-framesCenters[i]
      k2<-framesParents[i]
      framesBasedOnParents_G1[,,k1,k]<-rotateFrameToMainAxes(myFrame = frames_G1[,,k2,k],
                                                             vectors2rotate = frames_G1[,,k1,k])
    } 
  }
  print("Group 1 is done!")
  framesBasedOnParents_G2<-array(NA,dim = c(3,3,numberOfFrames,nSamplesG2))
  for (k in 1:nSamplesG2) {
    for (i in 1:numberOfFrames) {
      k1<-framesCenters[i]
      k2<-framesParents[i]
      framesBasedOnParents_G2[,,k1,k]<-rotateFrameToMainAxes(myFrame = frames_G2[,,k2,k],
                                                             vectors2rotate = frames_G2[,,k1,k])
    } 
  }
  print("Group 2 is done!")
}

# Calculate spokes directions based on their frames
if(TRUE){
  spokesDirectionsBasedOnFrames_G1<-array(NA,dim = c(nTotalRadii,3,nSamplesG1))
  for (k in 1:nSamplesG1) {
    for (i in 1:nTotalRadii) {
      spokeNo<-i # 1<=spokeNo<=nTotalRadii
      frameOfSpokeNo<-NA
      if(spokeNo<=upSpoeksNumber){
        frameOfSpokeNo<-spokeNo
      }else if(spokeNo<=2*upSpoeksNumber){
        frameOfSpokeNo<-(spokeNo-upSpoeksNumber)
      }else{ #crest
        frameOfSpokeNo<-(spokeNo-upSpoeksNumber)
      }
      
      spokesDirectionsBasedOnFrames_G1[i,,k]<-rotateFrameToMainAxes(myFrame = frames_G1[,,frameOfSpokeNo,k],
                                                                    vectors2rotate = spokeDirections_G1[i,,k])
      
    }
  }
  print("Group 1 is done!")
  spokesDirectionsBasedOnFrames_G2<-array(NA,dim = c(nTotalRadii,3,nSamplesG2))
  for (k in 1:nSamplesG2) {
    for (i in 1:nTotalRadii) {
      spokeNo<-i # 1<=spokeNo<=nTotalRadii
      frameOfSpokeNo<-NA
      if(spokeNo<=upSpoeksNumber){
        frameOfSpokeNo<-spokeNo
      }else if(spokeNo<=2*upSpoeksNumber){
        frameOfSpokeNo<-(spokeNo-upSpoeksNumber)
      }else{ #crest
        frameOfSpokeNo<-(spokeNo-upSpoeksNumber)
      }
      
      spokesDirectionsBasedOnFrames_G2[i,,k]<-rotateFrameToMainAxes(myFrame = frames_G2[,,frameOfSpokeNo,k],
                                                                    vectors2rotate = spokeDirections_G2[i,,k])
      
      
    }
  }
  print("Group 2 is done!")
}

# Calculate connection lengths and directions based on their frames
if(TRUE){
  
  connections_G1<-array(NA,dim = c(numberOfFrames,3,nSamplesG1))
  connectionsLengths_G1<-array(NA,dim=c(numberOfFrames,nSamplesG1))
  connectionsBasedOnParentFrames_G1<-array(NA,dim = dim(connections_G1))
  for (k in 1:nSamplesG1) {
    for (i in 1:numberOfFrames) {
      k1<-framesCenters[i]
      k2<-framesParents[i]
      tempVec<-skeletalSheet_G1[k1,,k]- skeletalSheet_G1[k2,,k]
      
      connectionsLengths_G1[k1,k]<-norm(tempVec,type = "2")
      
      if(norm(tempVec,type = "2")==0){
        connections_G1[k1,,k]<-c(0,0,0)
      }else{
        connections_G1[k1,,k]<-convertVec2unitVec(tempVec)
      }
      
      connectionsBasedOnParentFrames_G1[k1,,k]<-rotateFrameToMainAxes(myFrame = frames_G1[,,k2,k],
                                                                      vectors2rotate = connections_G1[k1,,k])
    }
  }
  print("Group 1 is done!")
  connections_G2<-array(NA,dim = c(numberOfFrames,3,nSamplesG2))
  connectionsLengths_G2<-array(NA,dim=c(numberOfFrames,nSamplesG2))
  connectionsBasedOnParentFrames_G2<-array(NA,dim = dim(connections_G2))
  for (k in 1:nSamplesG2) {
    for (i in 1:numberOfFrames) {
      k1<-framesCenters[i]
      k2<-framesParents[i]
      tempVec<-skeletalSheet_G2[k1,,k]- skeletalSheet_G2[k2,,k]
      
      connectionsLengths_G2[k1,k]<-norm(tempVec,type = "2")
      
      if(norm(tempVec,type = "2")==0){
        connections_G2[k1,,k]<-c(0,0,0)
      }else{
        connections_G2[k1,,k]<-convertVec2unitVec(tempVec)
      }
      
      connectionsBasedOnParentFrames_G2[k1,,k]<-rotateFrameToMainAxes(myFrame = frames_G2[,,k2,k],
                                                                      vectors2rotate = connections_G2[k1,,k])
    }
  }
  print("Group 2 is done!")
}

# calculate LP sizes
#NB! length of the center frame is 0 and must be excluded
LP_sizes_G1<-rep(NA,nSamplesG1)
for (i in 1:nSamplesG1) {
  LP_sizes_G1[i]<-exp(mean(c(log(connectionsLengths_G1[-16,i]),log(radii_G1[,i]))))
}
LP_sizes_G2<-rep(NA,nSamplesG2)
for (i in 1:nSamplesG2) {
  LP_sizes_G2[i]<-exp(mean(c(log(connectionsLengths_G2[-16,i]),log(radii_G2[,i]))))
}

#####################################################################################################
#####################################################################################################
#Removing or preserving the scale by LP-size based on the type of study

if(typeOfStudy=="sizeAndShapeAnalysis"){
  
  # sizes_G1 and sizes_G2 are LP size but we use them here for plot and scaling
  
  sizes_G1<-rep(1,nSamplesG1)
  sizes_G2<-rep(1,nSamplesG2)
  
  radiiScaled_G1<-radii_G1 #we don't have scaling in size-and-shape analysis
  radiiScaled_G2<-radii_G2
  
  connectionsLengthsScaled_G1<-connectionsLengths_G1
  connectionsLengthsScaled_G2<-connectionsLengths_G2
  
}else if(typeOfStudy=="shapeAnalysis"){
  
  #NB! length of the center frame is 0 and must be excluded
  sizes_G1<-rep(NA,nSamplesG1)
  for (i in 1:nSamplesG1) {
    sizes_G1[i]<-exp(mean(c(log(connectionsLengths_G1[-16,i]),log(radii_G1[,i]))))
  }
  sizes_G2<-rep(NA,nSamplesG2)
  for (i in 1:nSamplesG2) {
    sizes_G2[i]<-exp(mean(c(log(connectionsLengths_G2[-16,i]),log(radii_G2[,i]))))
  }
  
  radiiScaled_G1<-sweep(radii_G1, 2, sizes_G1, "/") #2 indicate operation "/" on columns
  radiiScaled_G2<-sweep(radii_G2, 2, sizes_G2, "/") 
  
  connectionsLengthsScaled_G1<-sweep(connectionsLengths_G1, 2, sizes_G1, "/") #2 indicate operation "/" on columns
  connectionsLengthsScaled_G2<-sweep(connectionsLengths_G2, 2, sizes_G2, "/") 
  
}

#####################################################################################################
#####################################################################################################
# Calculate mean LP-ds-rep based on Taheri & schulz 2022 article

#calculate mean frames in local and global coordinate systems
if(TRUE){
  
  framesBasedOnParentsVectorized_G1<-array(NA,dim = c(numberOfFrames,9,nSamplesG1))
  for (i in 1:nSamplesG1) {
    for (k in 1:numberOfFrames) {
      framesBasedOnParentsVectorized_G1[k,,i]<-as.vector(t(framesBasedOnParents_G1[,,k,i]))
    }
  }
  framesBasedOnParentsVectorized_G2<-array(NA,dim = c(numberOfFrames,9,nSamplesG2))
  for (i in 1:nSamplesG2) {
    for (k in 1:numberOfFrames) {
      framesBasedOnParentsVectorized_G2[k,,i]<-as.vector(t(framesBasedOnParents_G2[,,k,i]))
    }
  }
  
  meanFramesBasedOnParents_G1<-array(NA, dim = c(3,3,numberOfFrames))
  for (k in framesCenters) {
    if(k==16){
      meanFramesBasedOnParents_G1[,,k]<-rbind(c(0,0,1),c(1,0,0),c(0,1,0))
    }else{
      tempVec<-mean(as.SO3(t(framesBasedOnParentsVectorized_G1[k,,])),type = 'geometric')
      meanFramesBasedOnParents_G1[,,k]<-matrix(tempVec,nrow = 3,byrow = TRUE)
    }
  }
  meanFramesBasedOnParents_G2<-array(NA, dim = c(3,3,numberOfFrames))
  for (k in framesCenters) {
    if(k==16){
      meanFramesBasedOnParents_G2[,,k]<-rbind(c(0,0,1),c(1,0,0),c(0,1,0))
    }else{
      tempVec<-mean(as.SO3(t(framesBasedOnParentsVectorized_G2[k,,])),type = 'geometric')
      meanFramesBasedOnParents_G2[,,k]<-matrix(tempVec,nrow = 3,byrow = TRUE)
    }
  }
  
  meanFramesGlobalCoordinate_G1<-array(NA,dim = dim(meanFramesBasedOnParents_G1))
  meanFramesGlobalCoordinate_G1[,,framesCenters[1]]<-rbind(c(0,0,1),c(1,0,0),c(0,1,0))
  for (k in 2:numberOfFrames) {
    parent_Index<-framesParents[k]
    child_Index<-framesCenters[k]
    updatedParent<-meanFramesGlobalCoordinate_G1[,,parent_Index]
    meanFramesGlobalCoordinate_G1[,,child_Index]<-
      rotateFrameToMainAxesAndRotateBack(myFrame = updatedParent,
                                         vectorsInMainAxes = meanFramesBasedOnParents_G1[,,child_Index])
  }
  meanFramesGlobalCoordinate_G2<-array(NA,dim = dim(meanFramesBasedOnParents_G2))
  meanFramesGlobalCoordinate_G2[,,framesCenters[1]]<-rbind(c(0,0,1),c(1,0,0),c(0,1,0))
  for (k in 2:numberOfFrames) {
    parent_Index<-framesParents[k]
    child_Index<-framesCenters[k]
    updatedParent<-meanFramesGlobalCoordinate_G2[,,parent_Index]
    meanFramesGlobalCoordinate_G2[,,child_Index]<-
      rotateFrameToMainAxesAndRotateBack(myFrame = updatedParent,
                                         vectorsInMainAxes = meanFramesBasedOnParents_G2[,,child_Index])
  }
  print("Done!")
}

# Calculate mean spokes' directions based on frames
if(TRUE){
  meanSpokesDirectionsBasedOnFrames_G1<-array(NA,dim = c(nTotalRadii,3))
  for (i in 1:nTotalRadii) {
    # For extremely concentrated data we use Mardia mean direction 
    pcaTemp<-prcomp(t(spokesDirectionsBasedOnFrames_G1[i,,]))
    if(pcaTemp$sdev[1]<1e-02 | pcaTemp$sdev[2]<1e-02){
      meanSpokesDirectionsBasedOnFrames_G1[i,]<-convertVec2unitVec(colMeans(t(spokesDirectionsBasedOnFrames_G1[i,,])))
    }else if(typeOfMeanDirection=="Frechet"){
      meanSpokesDirectionsBasedOnFrames_G1[i,]<-frechetMean(spokesDirectionsBasedOnFrames_G1[i,,]) 
    }else if(typeOfMeanDirection=="PNS"){
      sphereType<-kurtosisTestFunction(spokesDirectionsBasedOnFrames_G1[i,,])
      meanSpokesDirectionsBasedOnFrames_G1[i,]<-pns(spokesDirectionsBasedOnFrames_G1[i,,],sphere.type = sphereType)$PNS$mean 
    }else{
      stop("Please specify the typeOfMeanDirection by PNS or Frechet!")
    }
  }
  meanSpokesDirectionsBasedOnFrames_G2<-array(NA,dim = c(nTotalRadii,3))
  for (i in 1:nTotalRadii) {
    # For extremely concentrated data we use Mardia mean direction 
    pcaTemp<-prcomp(t(spokesDirectionsBasedOnFrames_G2[i,,]))
    if(pcaTemp$sdev[1]<1e-02 | pcaTemp$sdev[2]<1e-02){
      meanSpokesDirectionsBasedOnFrames_G2[i,]<-convertVec2unitVec(colMeans(t(spokesDirectionsBasedOnFrames_G2[i,,])))
    }else if(typeOfMeanDirection=="Frechet"){
      meanSpokesDirectionsBasedOnFrames_G2[i,]<-frechetMean(spokesDirectionsBasedOnFrames_G2[i,,])
    }else if(typeOfMeanDirection=="PNS"){
      sphereType<-kurtosisTestFunction(spokesDirectionsBasedOnFrames_G2[i,,])
      meanSpokesDirectionsBasedOnFrames_G2[i,]<-pns(spokesDirectionsBasedOnFrames_G2[i,,],sphere.type = sphereType)$PNS$mean 
    }else{
      stop("Please specify the typeOfMeanDirection by PNS or Frechet!")
    }
  }
  print("Done!")
}

# Calculate mean spokes' directions based on global coordinate (using mean frames in global coordinate)
if(TRUE){
  meanSpokesDirectionsGlobalCoordinate_G1<-array(NA,dim = c(nTotalRadii,3))
  for (i in 1:nTotalRadii) {
    spokeNo<-i # 1<=spokeNo<=nTotalRadii
    frameOfSpokeNo<-NA
    if(spokeNo<=upSpoeksNumber){
      frameOfSpokeNo<-spokeNo
    }else if(spokeNo<=2*upSpoeksNumber){
      frameOfSpokeNo<-(spokeNo-upSpoeksNumber)
    }else{ #crest
      frameOfSpokeNo<-(spokeNo-upSpoeksNumber)
    }
    meanSpokesDirectionsGlobalCoordinate_G1[i,]<-
      rotateFrameToMainAxesAndRotateBack(myFrame = meanFramesGlobalCoordinate_G1[,,frameOfSpokeNo],
                                         vectorsInMainAxes = meanSpokesDirectionsBasedOnFrames_G1[i,])
    
  }
  meanSpokesDirectionsGlobalCoordinate_G2<-array(NA,dim = c(nTotalRadii,3))
  for (i in 1:nTotalRadii) {
    spokeNo<-i # 1<=spokeNo<=nTotalRadii
    frameOfSpokeNo<-NA
    if(spokeNo<=upSpoeksNumber){
      frameOfSpokeNo<-spokeNo
    }else if(spokeNo<=2*upSpoeksNumber){
      frameOfSpokeNo<-(spokeNo-upSpoeksNumber)
    }else{ #crest
      frameOfSpokeNo<-(spokeNo-upSpoeksNumber)
    }
    meanSpokesDirectionsGlobalCoordinate_G2[i,]<-
      rotateFrameToMainAxesAndRotateBack(myFrame = meanFramesGlobalCoordinate_G2[,,frameOfSpokeNo],
                                         vectorsInMainAxes = meanSpokesDirectionsBasedOnFrames_G2[i,])
    
  }
}

# Calculate geometric mean of spokes' lengths
if(TRUE){
  radiiMean_G1<-exp(rowMeans(log(radiiScaled_G1)))
  radiiMean_G2<-exp(rowMeans(log(radiiScaled_G2)))
}

# Calculate geometric mean of connections' lengths
if(TRUE){
  meanConnectionsLengths_G1<-exp(rowMeans(log(connectionsLengthsScaled_G1)))
  meanConnectionsLengths_G2<-exp(rowMeans(log(connectionsLengthsScaled_G2)))
}

# Calculate mean connection directions based on frames
if(TRUE){
  meanConnectionsBasedOnParentFrames_G1<-array(NA,dim = c(numberOfFrames,3))
  for (i in 1:numberOfFrames) {
    if(i==framesCenters[1]){
      meanConnectionsBasedOnParentFrames_G1[i,]<-c(0,0,0)
    }else{
      # For extremely concentrated data we use Mardia mean direction 
      pcaTemp<-prcomp(t(connectionsBasedOnParentFrames_G1[i,,]))
      if(pcaTemp$sdev[1]<1e-02 | pcaTemp$sdev[2]<1e-02){
        meanConnectionsBasedOnParentFrames_G1[i,]<-convertVec2unitVec(colMeans(t(connectionsBasedOnParentFrames_G1[i,,])))
      }else if(typeOfMeanDirection=="Frechet"){
        meanConnectionsBasedOnParentFrames_G1[i,]<-frechetMean(connectionsBasedOnParentFrames_G1[i,,]) 
      }else if(typeOfMeanDirection=="PNS"){
        sphereType<-kurtosisTestFunction(connectionsBasedOnParentFrames_G1[i,,])
        meanConnectionsBasedOnParentFrames_G1[i,]<-pns(connectionsBasedOnParentFrames_G1[i,,],sphere.type = sphereType)$PNS$mean  
      }else{
        stop("Please specify the typeOfMeanDirection by PNS or Frechet")
      }
    }
  }
  meanConnectionsBasedOnParentFrames_G2<-array(NA,dim = c(numberOfFrames,3))
  for (i in 1:numberOfFrames) {
    if(i==framesCenters[1]){
      meanConnectionsBasedOnParentFrames_G2[i,]<-c(0,0,0)
    }else{
      # For extremely concentrated data we use Mardia mean direction 
      pcaTemp<-prcomp(t(connectionsBasedOnParentFrames_G2[i,,]))
      if(pcaTemp$sdev[1]<1e-02 | pcaTemp$sdev[2]<1e-02){
        meanConnectionsBasedOnParentFrames_G2[i,]<-convertVec2unitVec(colMeans(t(connectionsBasedOnParentFrames_G2[i,,])))
      }else if(typeOfMeanDirection=="Frechet"){
        meanConnectionsBasedOnParentFrames_G2[i,]<-frechetMean(connectionsBasedOnParentFrames_G2[i,,]) 
      }else if(typeOfMeanDirection=="PNS"){
        sphereType<-kurtosisTestFunction(connectionsBasedOnParentFrames_G2[i,,])
        meanConnectionsBasedOnParentFrames_G2[i,]<-pns(connectionsBasedOnParentFrames_G2[i,,],sphere.type = sphereType)$PNS$mean  
      }else{
        stop("Please specify the typeOfMeanDirection by PNS or Frechet")
      }
    }
  }
  print("Done!")
}

# Calculate mean connection based on global coordinate (using mean frames in global coordinate)
if(TRUE){
  meanConnectionsGlobalCoordinate_G1<-array(NA,dim = c(numberOfFrames,3))
  for (i in 1:numberOfFrames) {
    k1<-framesCenters[i]
    k2<-framesParents[i]
    meanConnectionsGlobalCoordinate_G1[k1,]<-
      rotateFrameToMainAxesAndRotateBack(myFrame = meanFramesGlobalCoordinate_G1[,,k2],
                                         vectorsInMainAxes = meanConnectionsBasedOnParentFrames_G1[k1,])
  }
  meanConnectionsGlobalCoordinate_G2<-array(NA,dim = c(numberOfFrames,3))
  for (i in 1:numberOfFrames) {
    k1<-framesCenters[i]
    k2<-framesParents[i]
    meanConnectionsGlobalCoordinate_G2[k1,]<-
      rotateFrameToMainAxesAndRotateBack(myFrame = meanFramesGlobalCoordinate_G2[,,k2],
                                         vectorsInMainAxes = meanConnectionsBasedOnParentFrames_G2[k1,])
  }
}

#####################################################################################################
#####################################################################################################
# Convert mean LP-ds-rep to a GP-ds-rep

if(TRUE){
  meanPositions_G1<-array(NA,dim = c(numberOfFrames,3))
  meanPositions_G1[16,]<-c(0,0,0)
  for (i in 2:numberOfFrames) {
    k1<-framesCenters[i]
    k2<-framesParents[i]
    meanPositions_G1[k1,]<-meanPositions_G1[k2,]+
      meanConnectionsLengths_G1[k1]*meanConnectionsGlobalCoordinate_G1[k1,]
    
  }
  meanPositions_G2<-array(NA,dim = c(numberOfFrames,3))
  meanPositions_G2[16,]<-c(0,0,0)
  for (i in 2:numberOfFrames) {
    k1<-framesCenters[i]
    k2<-framesParents[i]
    meanPositions_G2[k1,]<-meanPositions_G2[k2,]+
      meanConnectionsLengths_G2[k1]*meanConnectionsGlobalCoordinate_G2[k1,]
    
  }
  meanSpokesTails_G1<-array(NA,dim = c(nTotalRadii,3))
  meanSpokesTips_G1<-array(NA,dim = c(nTotalRadii,3))
  for (i in 1:nTotalRadii) {
    frameOfSpokeNo<-NA
    if(i<=upSpoeksNumber){
      frameOfSpokeNo<-i
    }else if(i<=2*upSpoeksNumber){
      frameOfSpokeNo<-(i-upSpoeksNumber)
    }else{ #crest
      frameOfSpokeNo<-(i-upSpoeksNumber)
    }
    
    meanSpokesTails_G1[i,]<-meanPositions_G1[frameOfSpokeNo,]
    meanSpokesTips_G1[i,]<-meanPositions_G1[frameOfSpokeNo,]+meanSpokesDirectionsGlobalCoordinate_G1[i,]*radiiMean_G1[i]
  }
  meanSpokesTails_G2<-array(NA,dim = c(nTotalRadii,3))
  meanSpokesTips_G2<-array(NA,dim = c(nTotalRadii,3))
  for (i in 1:nTotalRadii) {
    frameOfSpokeNo<-NA
    if(i<=upSpoeksNumber){
      frameOfSpokeNo<-i
    }else if(i<=2*upSpoeksNumber){
      frameOfSpokeNo<-(i-upSpoeksNumber)
    }else{ #crest
      frameOfSpokeNo<-(i-upSpoeksNumber)
    }
    
    meanSpokesTails_G2[i,]<-meanPositions_G2[frameOfSpokeNo,]
    meanSpokesTips_G2[i,]<-meanPositions_G2[frameOfSpokeNo,]+meanSpokesDirectionsGlobalCoordinate_G2[i,]*radiiMean_G2[i]
  }
  print("Done!")
}

# Plot overlaid LP-ds-rep means 
if(TRUE){
  open3d()
  srep1<-rbind(meanSpokesTails_G1,meanSpokesTips_G1)*mean(sizes_G1) #we scale back to the original size by *mean(sizes_G1)
  for (i in 1:nTotalRadii) {
    plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 1.5,col = "blue",expand = 10,box=FALSE,add = TRUE)
  }
  srep2<-rbind(meanSpokesTails_G2,meanSpokesTips_G2)*mean(sizes_G2)
  for (i in 1:nTotalRadii) {
    plot3d(srep2[c(i,(i+nTotalRadii)),],type="l",lwd = 1.5,col = "red",expand = 10,box=FALSE,add = TRUE)
  }
  # legend3d("topright", legend = paste(c('Mean G1', 'Mean G2')), pch = 16, col = c("blue","red"), cex=1, inset=c(0.02))
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             # xlim = c(-10,15),ylim =c(-20,20),zlim = c(-5,5),
             box = TRUE, axes = TRUE, main = "Overlaid mean shapes", sub = NULL,
             top = TRUE, aspect = FALSE, expand = 1.1)
}

#####################################################################################################
#####################################################################################################
# Hypothesis testing

# hypothesis test on LP size
pValues_LP_sizes<-meanDifferenceTest1D(log(LP_sizes_G1),log(LP_sizes_G2),type = typeOfTest) 
cat("pValue of LP sizes is:",pValues_LP_sizes,"\n")
boxplot(LP_sizes_G1, LP_sizes_G2, names = c("CG","SCZ"),main="LP-size")
cat("sd LP size G1:",sd(LP_sizes_G1),"sd LP size G2:",sd(LP_sizes_G2),"\n")
cat("Mean LP size G1:",mean(LP_sizes_G1),"mean LP size G2:",mean(LP_sizes_G2),"\n")


# hypothesis thresholds to ignore extremely concentrated data
if(TRUE){
  thresholdDirections<-pi/200
  thresholdDirections
  thresholdLengths<-mean(c(LP_sizes_G1,LP_sizes_G2))/200
  thresholdLengths  
}


# hypothesis test on spokes' lengths
pValues_TtestRadii<-rep(NA,nTotalRadii)
pb <- txtProgressBar(min = 0, max = nTotalRadii, style = 3) #progress bar
for (i in 1:nTotalRadii) {
  setTxtProgressBar(pb, i) #create progress bar
  
  if(abs(mean(radiiScaled_G1)-mean(radiiScaled_G2))<thresholdLengths){
    pValues_TtestRadii[i]<-1
  }else{
    pValues_TtestRadii[i]<-meanDifferenceTest1D(log(radiiScaled_G1[i,]),
                                                log(radiiScaled_G2[i,]),
                                                type = typeOfTest)  
  }
  
}
# which(pValues_TtestRadii<=0.05)


# hypothesis test on connections' length
pValues_TtestConnectionsLengths<-rep(NA,numberOfFrames)
pb <- txtProgressBar(min = 0, max = numberOfFrames, style = 3) #progress bar
for (i in 1:numberOfFrames) {
  setTxtProgressBar(pb, i) #create progress bar
  if(i==16){
    pValues_TtestConnectionsLengths[i]<-1
  }else if(abs(mean(connectionsLengthsScaled_G1[i,])-mean(connectionsLengthsScaled_G2[i,]))<thresholdLengths){
    pValues_TtestConnectionsLengths[i]<-1
  }else{
    pValues_TtestConnectionsLengths[i]<-meanDifferenceTest1D(log(connectionsLengthsScaled_G1[i,]),
                                                             log(connectionsLengthsScaled_G2[i,]),
                                                             type = typeOfTest)
  }
}
# which(pValues_TtestConnectionsLengths<=0.05)


# hypothesis test on spokes' directions based on local frames
euclideanizedSpokesDirBasedOnFramesG1<-array(NA,dim = c(nSamplesG1,2,nTotalRadii))
euclideanizedSpokesDirBasedOnFramesG2<-array(NA,dim = c(nSamplesG2,2,nTotalRadii))
pValspokesDirectionsBasedOnFrames<-rep(NA,nTotalRadii)
pb <- txtProgressBar(min = 0, max = nTotalRadii, style = 3) #progress bar
for(i in 1:nTotalRadii){
  setTxtProgressBar(pb, i) #create progress bar
  #NB! euclideanization must contain two groups because it uses the pooled mean 
  euclideanizedTemp<-euclideanization(spokesDirectionsBasedOnFrames_G1[i,,],
                                      spokesDirectionsBasedOnFrames_G2[i,,],
                                      type = typeOfStudy4directions)
  
  if(norm(colMeans(euclideanizedTemp$euclideanG1)-colMeans(euclideanizedTemp$euclideanG2),type = '2')<thresholdDirections){
    pValspokesDirectionsBasedOnFrames[i]<-1
  }else{
    pValspokesDirectionsBasedOnFrames[i]<-meanDifferenceTestMultivariate(euclideanizedTemp$euclideanG1,
                                                                         euclideanizedTemp$euclideanG2,
                                                                         type=typeOfTest)
    
  }
  
  euclideanizedSpokesDirBasedOnFramesG1[,,i]<-euclideanizedTemp$euclideanG1
  euclideanizedSpokesDirBasedOnFramesG2[,,i]<-euclideanizedTemp$euclideanG2
  
}
# which(pValspokesDirectionsBasedOnFrames<=0.05)


# hypothesis test on connections' directions based on local frames
euclideanizedConnectionsBasedOnParentFramesG1<-array(0,dim = c(nSamplesG1,2,numberOfFrames))
euclideanizedConnectionsBasedOnParentFramesG2<-array(0,dim = c(nSamplesG2,2,numberOfFrames))
pValConnectionsBasedOnParentFrames<-rep(NA,numberOfFrames)
pb <- txtProgressBar(min = 0, max = numberOfFrames, style = 3) #progress bar
for(i in 1:numberOfFrames){
  setTxtProgressBar(pb, i) #create progress bar
  if(i==16){
    pValConnectionsBasedOnParentFrames[i]<-1
    next
  }
  
  euclideanizedTemp<-euclideanization(connectionsBasedOnParentFrames_G1[i,,],
                                      connectionsBasedOnParentFrames_G2[i,,],
                                      type = typeOfStudy4directions)
  
  if(norm(colMeans(euclideanizedTemp$euclideanG1)-colMeans(euclideanizedTemp$euclideanG2),type = '2')<thresholdDirections){
    
    pValConnectionsBasedOnParentFrames[i]<-1
    
  }else{
    
    pValConnectionsBasedOnParentFrames[i]<-meanDifferenceTestMultivariate(euclideanizedTemp$euclideanG1,
                                                                          euclideanizedTemp$euclideanG2,
                                                                          type=typeOfTest)
    
  }
  
  euclideanizedConnectionsBasedOnParentFramesG1[,,i]<-euclideanizedTemp$euclideanG1
  euclideanizedConnectionsBasedOnParentFramesG2[,,i]<-euclideanizedTemp$euclideanG2
  
}
# which(pValConnectionsBasedOnParentFrames<=0.05)

framesBasedOnParentsVectorized_G1<-array(NA,dim = c(numberOfFrames,9,nSamplesG1))
for (i in 1:nSamplesG1) {
  for (k in 1:numberOfFrames) {
    framesBasedOnParentsVectorized_G1[k,,i]<-as.vector(t(framesBasedOnParents_G1[,,k,i]))
  }
}
framesBasedOnParentsVectorized_G2<-array(NA,dim = c(numberOfFrames,9,nSamplesG2))
for (i in 1:nSamplesG2) {
  for (k in 1:numberOfFrames) {
    framesBasedOnParentsVectorized_G2[k,,i]<-as.vector(t(framesBasedOnParents_G2[,,k,i]))
  }
}

# hypothesis test on frames' normal directions based on parent frames
euclideanizedFrameBasedOnParentG1<-array(0,dim = c(nSamplesG1,3,numberOfFrames))
euclideanizedFrameBasedOnParentG2<-array(0,dim = c(nSamplesG2,3,numberOfFrames))
pValFramesBasedOnParent<-rep(NA,numberOfFrames)
pb <- txtProgressBar(min = 0, max = numberOfFrames, style = 3) #progress bar
for(i in 1:numberOfFrames){
  setTxtProgressBar(pb, i) #create progress bar
  if(i==16){
    pValFramesBasedOnParent[i]<-1
    next
  }
  
  Q4Temp<-as.Q4(as.SO3(t(framesBasedOnParentsVectorized_G1[i,,])))
  Q4Temp2<-matrix(as.numeric(t(Q4Temp)),ncol = 4,byrow = TRUE)
  for (j in 1:dim(Q4Temp2)[1]) {
    Q4Temp2[j,]<-Q4Temp2[j,]/norm(Q4Temp2[j,],type = '2')
  }
  Q4Temp3<-as.Q4(as.SO3(t(framesBasedOnParentsVectorized_G2[i,,])))
  Q4Temp4<-matrix(as.numeric(t(Q4Temp3)),ncol = 4,byrow = TRUE)
  for (j in 1:dim(Q4Temp4)[1]) {
    Q4Temp4[j,]<-Q4Temp4[j,]/norm(Q4Temp4[j,],type = '2')
  }
  
  euclideanizedTemp<-euclideanization(t(Q4Temp2),t(Q4Temp4),type = typeOfStudy4directions)
  
  
  if(norm(colMeans(euclideanizedTemp$euclideanG1)-colMeans(euclideanizedTemp$euclideanG2),type = '2')<thresholdDirections){
    pValFramesBasedOnParent[i]<-1
  }else{
    pValFramesBasedOnParent[i]<-meanDifferenceTestMultivariate(euclideanizedTemp$euclideanG1,
                                                               euclideanizedTemp$euclideanG2,
                                                               type=typeOfTest) 
  }
  
  euclideanizedFrameBasedOnParentG1[,,i]<-euclideanizedTemp$euclideanG1
  euclideanizedFrameBasedOnParentG2[,,i]<-euclideanizedTemp$euclideanG2
  
}
# which(pValFramesBasedOnParent<=0.05)


# plot significant GOPs
pvalues_LP_ds_rep <- c(pValues_TtestRadii,                 
                       pValues_TtestConnectionsLengths,    
                       pValspokesDirectionsBasedOnFrames,  
                       pValConnectionsBasedOnParentFrames, 
                       pValFramesBasedOnParent,            
                       pValues_LP_sizes)

n_s<-nTotalRadii
n_f<-numberOfFrames

length(pvalues_LP_ds_rep)
alpha<-0.05
significantPvalues<-which(pvalues_LP_ds_rep<=alpha)
significantPvalues

#adjust p-values by Benjamini-Hochberg
FDR<-0.15
pvalues_LP_ds_rep_BH<-p.adjust(pvalues_LP_ds_rep,method = "BH")
pvalues_LP_ds_rep_Bonferroni<-p.adjust(pvalues_LP_ds_rep,method = "bonferroni")
significantPvalues_BH<-which(pvalues_LP_ds_rep_BH<=FDR)
significantPvalues_BH
significantPvalues_Bonferroni<-which(pvalues_LP_ds_rep_Bonferroni<=FDR)
significantPvalues_Bonferroni

cat("\n","Percentage of sig raw p-value is:",length(significantPvalues)/length(pvalues_LP_ds_rep),"\n")
cat("\n","Percentage of BH adjusted p-value is:",length(significantPvalues_BH)/length(pvalues_LP_ds_rep),"\n")



# plot by ggplot
df_LP <- data.frame(Type=c(rep("Raw p-value",length(pvalues_LP_ds_rep)),
                           rep("Bonferroni",length(pvalues_LP_ds_rep)),
                           rep("BH",length(pvalues_LP_ds_rep))),
                    ordereOfPvalues=1:length(pvalues_LP_ds_rep),
                    Values=c(sort(pvalues_LP_ds_rep),sort(pvalues_LP_ds_rep_Bonferroni),sort(pvalues_LP_ds_rep_BH)))
p<-ggplot(df_LP, aes(x=ordereOfPvalues, y=Values, group=Type))
p + geom_line(aes(linetype=Type),size=1)+
  # geom_line(aes(linetype=Type, col=Type),size=1)+
  # geom_point(aes(shape=Type),alpha=0.7)+
  geom_hline(yintercept=alpha,linetype="solid", color = "red")+
  geom_hline(yintercept=FDR,linetype="solid", color = "blue")+
  scale_linetype_manual(values=c("solid","dotdash", "dotted")) +
  theme_bw()+
  theme(plot.title = element_text(size = 17, hjust = 0.5),
        legend.text=element_text(size=17),
        legend.title=element_blank(),
        # panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        # legend.title=element_text(size=17),
        legend.position="bottom",
        axis.title=element_text(size=17),
        legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
  guides(colour = guide_legend(title.hjust = 0.5))+
  xlab("Ranking of p-values") + ylab("p-values")


#1
significantRadii<-significantPvalues[which(significantPvalues<=n_s)]
significantRadii
significantRadii_BH<-significantPvalues_BH[which(significantPvalues_BH<=n_s)]
significantRadii_BH
#2
significantConnectionsLengths<-significantPvalues[which(n_s+1<=significantPvalues
                                                        & significantPvalues<=(n_s+n_f))]-n_s
significantConnectionsLengths
significantConnectionsLengths_BH<-significantPvalues_BH[which((n_s+1)<=significantPvalues_BH
                                                              & significantPvalues_BH<=(n_s+n_f))]-n_s
significantConnectionsLengths_BH
#3
significantspokesDirections<-significantPvalues[which((n_s+n_f+1)<=significantPvalues
                                                      & significantPvalues<=(2*n_s+n_f))]-(n_s+n_f)
significantspokesDirections
significantspokesDirections_BH<-significantPvalues_BH[which((n_s+n_f+1)<=significantPvalues_BH
                                                            & significantPvalues_BH<=(2*n_s+n_f))]-(n_s+n_f)
significantspokesDirections_BH
#4
significantConnectionsDirections<-significantPvalues[which((2*n_s+n_f+1)<=significantPvalues
                                                           & significantPvalues<=(2*n_s+2*n_f))]-(2*n_s+n_f)
significantConnectionsDirections
significantConnectionsDirections_BH<-significantPvalues_BH[which((2*n_s+n_f+1)<=significantPvalues_BH
                                                                 & significantPvalues_BH<=(2*n_s+2*n_f))]-(2*n_s+n_f)
significantConnectionsDirections_BH
#5
significantFrame<-significantPvalues[which((2*n_s+2*n_f+1)<=significantPvalues
                                           & significantPvalues<=(2*n_s+3*n_f))]-(2*n_s+2*n_f)
significantFrame
significantFrame_BH<-significantPvalues_BH[which((2*n_s+2*n_f+1)<=significantPvalues_BH
                                                 & significantPvalues_BH<=(2*n_s+3*n_f))]-(2*n_s+2*n_f)
significantFrame_BH


#plot significant GOPs before and after the BH adjustment
if(TRUE){
  
  #1 plot
  srep1<-rbind(meanSpokesTails_G1,meanSpokesTips_G1)*mean(sizes_G1)
  open3d()
  for (i in 1:nTotalRadii) {
    plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 1.5,col = "blue",expand = 10,box=FALSE,add = TRUE)
  }
  for (i in significantRadii) {
    plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 7,col = "red",expand = 10,box=FALSE,add = TRUE)
  }
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE,
             # main = "Significant spokes' lengths",
             sub = NULL, top = T, aspect = FALSE, expand = 1.1)
  # pp <- par3d(no.readonly=TRUE)
  # dput(pp, file="plotView.R", control = "all")
  # pp <- dget("plotView.R")
  # par3d(userMatrix=pp$userMatrix)
  
  #with correction
  open3d()
  for (i in 1:nTotalRadii) {
    plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 1.5,col = "blue",expand = 10,box=FALSE,add = TRUE)
  }
  for (i in significantRadii_BH) {
    plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 7,col = "red",expand = 10,box=FALSE,add = TRUE)
  }
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE,
             # main = "Significant spokes' lengths after BH adjustment",
             sub = NULL, top = T, aspect = FALSE, expand = 1.1)
  # pp <- dget("plotView.R")
  # par3d(userMatrix=pp$userMatrix)
  
  #2 plot
  open3d()
  skelG1_1<-meanSpokesTails_G1[skelRange,]*mean(sizes_G1)
  for (i in 2:numberOfFrames) {
    k1<-framesCenters[i]
    k2<-framesParents[i]
    vectors3d(skelG1_1[k1,],origin = skelG1_1[k2,],headlength = 0.1,radius = 1/6, col="blue", lwd=1)
  }
  for (i in significantConnectionsLengths) {
    k1<-i
    k2<-framesParents[which(framesCenters==i)]
    vectors3d(skelG1_1[k1,],origin = skelG1_1[k2,],headlength = 0.1,radius = 1/6, col="red", lwd=3)
  }
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE,
             # main = "Significant connections' lengths",
             sub = NULL, top = T, aspect = FALSE, expand = 1.1)
  # pp <- dget("plotView.R")
  # par3d(userMatrix=pp$userMatrix)
  
  #with correction
  open3d()
  for (i in 2:numberOfFrames) {
    k1<-framesCenters[i]
    k2<-framesParents[i]
    vectors3d(skelG1_1[k1,],origin = skelG1_1[k2,],headlength = 0.1,radius = 1/6, col="blue", lwd=1)
  }
  for (i in significantConnectionsLengths_BH) {
    k1<-i
    k2<-framesParents[which(framesCenters==i)]
    vectors3d(skelG1_1[k1,],origin = skelG1_1[k2,],headlength = 0.1,radius = 1/6, col="red", lwd=3)
  }
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE,
             # main = "Significant connections' lengths after BH adjustment",
             sub = NULL, top = T, aspect = FALSE, expand = 1.1)
  # pp <- dget("plotView.R")
  # par3d(userMatrix=pp$userMatrix)
  
  #3 plot
  srep1<-rbind(meanSpokesTails_G1,meanSpokesTips_G1)*mean(sizes_G1)
  open3d()
  for (i in 1:nTotalRadii) {
    plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 1.5,col = "blue",expand = 10,box=FALSE,add = TRUE)
  }
  for (i in significantspokesDirections) {
    plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 7,col = "red",expand = 10,box=FALSE,add = TRUE)
  }
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE,
             # main = "Significant spokes' directions",
             sub = NULL, top = T, aspect = FALSE, expand = 1.1)
  # pp <- dget("plotView.R")
  # par3d(userMatrix=pp$userMatrix)
  
  #with correction
  open3d()
  for (i in 1:nTotalRadii) {
    plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 1.5,col = "blue",expand = 10,box=FALSE,add = TRUE)
  }
  for (i in significantspokesDirections_BH) {
    plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 7,col = "red",expand = 10,box=FALSE,add = TRUE)
  }
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE,
             # main = "Significant spokes' directions after BH adjustment",
             sub = NULL, top = T, aspect = FALSE, expand = 1.1)
  # pp <- dget("plotView.R")
  # par3d(userMatrix=pp$userMatrix)
  
  #4 plot
  open3d()
  skelG1_1<-meanSpokesTails_G1[skelRange,]*mean(sizes_G1)
  for (i in 2:numberOfFrames) {
    k1<-framesCenters[i]
    k2<-framesParents[i]
    vectors3d(skelG1_1[k1,],origin = skelG1_1[k2,],headlength = 0.1,radius = 1/6, col="blue", lwd=1)
  }
  for (i in significantConnectionsDirections) {
    k1<-i
    k2<-framesParents[which(framesCenters==i)]
    vectors3d(skelG1_1[k1,],origin = skelG1_1[k2,],headlength = 0.1,radius = 1/6, col="red", lwd=3)
  }
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE,
             # main = "Significant connections' directions",
             sub = NULL, top = T, aspect = FALSE, expand = 1.1)
  # pp <- dget("plotView.R")
  # par3d(userMatrix=pp$userMatrix)
  
  #with correction
  open3d()
  for (i in 2:numberOfFrames) {
    k1<-framesCenters[i]
    k2<-framesParents[i]
    vectors3d(skelG1_1[k1,],origin = skelG1_1[k2,],headlength = 0.1,radius = 1/6, col="blue", lwd=1)
  }
  for (i in significantConnectionsDirections_BH) {
    k1<-i
    k2<-framesParents[which(framesCenters==i)]
    vectors3d(skelG1_1[k1,],origin = skelG1_1[k2,],headlength = 0.1,radius = 1/6, col="red", lwd=3)
  }
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE,
             # main = "Significant connections' directions after BH adjustment",
             sub = NULL, top = T, aspect = FALSE, expand = 1.1)
  # pp <- dget("plotView.R")
  # par3d(userMatrix=pp$userMatrix)
  
  #5 plot
  open3d()
  skeletalSheet<-meanPositions_G1*mean(sizes_G1)
  for (i in 2:numberOfFrames) {
    vectors3d(skeletalSheet[i,]+meanFramesGlobalCoordinate_G1[1,,i],origin = skeletalSheet[i,],headlength = 0.1,radius = 1/10, col="darkblue", lwd=2)
    vectors3d(skeletalSheet[i,]+meanFramesGlobalCoordinate_G1[2,,i],origin = skeletalSheet[i,],headlength = 0.1,radius = 1/10, col="blue", lwd=2)
    vectors3d(skeletalSheet[i,]+meanFramesGlobalCoordinate_G1[3,,i],origin = skeletalSheet[i,],headlength = 0.1,radius = 1/10, col="lightblue", lwd=2)
  }
  for (i in significantFrame) {
    vectors3d(skeletalSheet[i,]+meanFramesGlobalCoordinate_G1[1,,i],origin = skeletalSheet[i,],headlength = 0.1,radius = 1/10, col="red", lwd=7)
    vectors3d(skeletalSheet[i,]+meanFramesGlobalCoordinate_G1[2,,i],origin = skeletalSheet[i,],headlength = 0.1,radius = 1/10, col="red", lwd=7)
    vectors3d(skeletalSheet[i,]+meanFramesGlobalCoordinate_G1[3,,i],origin = skeletalSheet[i,],headlength = 0.1,radius = 1/10, col="red", lwd=7)
  }
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE,
             # main = "Significant frames",
             sub = NULL, top = T, aspect = FALSE, expand = 1.1)
  # pp <- dget("plotView.R")
  # par3d(userMatrix=pp$userMatrix)
  
  #with correction
  open3d()
  for (i in 2:numberOfFrames) {
    vectors3d(skeletalSheet[i,]+meanFramesGlobalCoordinate_G1[1,,i],origin = skeletalSheet[i,],headlength = 0.1,radius = 1/10, col="darkblue", lwd=2)
    vectors3d(skeletalSheet[i,]+meanFramesGlobalCoordinate_G1[2,,i],origin = skeletalSheet[i,],headlength = 0.1,radius = 1/10, col="blue", lwd=2)
    vectors3d(skeletalSheet[i,]+meanFramesGlobalCoordinate_G1[3,,i],origin = skeletalSheet[i,],headlength = 0.1,radius = 1/10, col="lightblue", lwd=2)
  }
  for (i in significantFrame_BH) {
    vectors3d(skeletalSheet[i,]+meanFramesGlobalCoordinate_G1[1,,i],origin = skeletalSheet[i,],headlength = 0.1,radius = 1/10, col="red", lwd=7)
    vectors3d(skeletalSheet[i,]+meanFramesGlobalCoordinate_G1[2,,i],origin = skeletalSheet[i,],headlength = 0.1,radius = 1/10, col="red", lwd=7)
    vectors3d(skeletalSheet[i,]+meanFramesGlobalCoordinate_G1[3,,i],origin = skeletalSheet[i,],headlength = 0.1,radius = 1/10, col="red", lwd=7)
  }
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE,
             # main = "Significant frames after BH adjustment",
             sub = NULL, top = T, aspect = FALSE, expand = 1.1)
  # pp <- dget("plotView.R")
  # par3d(userMatrix=pp$userMatrix)
  
}


#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################
############################ LP-dss-rep (swept skeletal structure) ##################################
#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################


# load dss-reps of the ellipsoids and inflated ellipsoids 

load("files/LP_dss_rep_inflationSimulation_G1.Rdata")
load("files/LP_dss_rep_inflationSimulation_G2.Rdata")

# convert dss-reps to LP-dss-reps

# directions based on frames
nSamplesG1<-length(LP_ds_rep_inflationSimulation_G1)
nSamplesG2<-length(LP_ds_rep_inflationSimulation_G2)
framesCenters<-LP_ds_rep_inflationSimulation_G1[[1]]$framesCenters
parentsIndices<-LP_ds_rep_inflationSimulation_G1[[1]]$parentsIndices
parentsIndices<-LP_ds_rep_inflationSimulation_G1[[1]]$parentsIndices
pointsIndices<-LP_ds_rep_inflationSimulation_G1[[1]]$pointsIndices

frameIndices<-which(!is.na(framesCenters))

numberOfSkelPoints<-length(pointsIndices)


# medial points
medialPoints3D_G1<-array(NA,dim = c(numberOfSkelPoints,3,nSamplesG1))
for (i in 1:nSamplesG1) {
  medialPoints3D_G1[,,i]<-LP_ds_rep_inflationSimulation_G1[[i]]$medialPoints3D
}
medialPoints3D_G2<-array(NA,dim = c(numberOfSkelPoints,3,nSamplesG2))
for (i in 1:nSamplesG2) {
  medialPoints3D_G2[,,i]<-LP_ds_rep_inflationSimulation_G2[[i]]$medialPoints3D
}

#frames
framesFirstVectors_G1<-array(NA,dim = c(numberOfSkelPoints,3,nSamplesG1))
framesSecondVectors_G1<-array(NA,dim = c(numberOfSkelPoints,3,nSamplesG1))
framesThirdVectors_G1<-array(NA,dim = c(numberOfSkelPoints,3,nSamplesG1))
for (i in 1:nSamplesG1) {
  framesFirstVectors_G1[,,i]<-LP_ds_rep_inflationSimulation_G1[[i]]$framesFirstVectors
  framesSecondVectors_G1[,,i]<-LP_ds_rep_inflationSimulation_G1[[i]]$framesSecondVectors
  framesThirdVectors_G1[,,i]<-LP_ds_rep_inflationSimulation_G1[[i]]$framesThirdVectors
}
framesFirstVectors_G2<-array(NA,dim = c(numberOfSkelPoints,3,nSamplesG2))
framesSecondVectors_G2<-array(NA,dim = c(numberOfSkelPoints,3,nSamplesG2))
framesThirdVectors_G2<-array(NA,dim = c(numberOfSkelPoints,3,nSamplesG2))
for (i in 1:nSamplesG2) {
  framesFirstVectors_G2[,,i]<-LP_ds_rep_inflationSimulation_G2[[i]]$framesFirstVectors
  framesSecondVectors_G2[,,i]<-LP_ds_rep_inflationSimulation_G2[[i]]$framesSecondVectors
  framesThirdVectors_G2[,,i]<-LP_ds_rep_inflationSimulation_G2[[i]]$framesThirdVectors
}


# Combine frames vectors to make SO(3) frames in global coordinate system
frames_G1<-array(NA,dim = c(3,3,numberOfSkelPoints,nSamplesG1))
for (k in 1:nSamplesG1) {
  for (i in 1:numberOfSkelPoints) {
    frames_G1[,,i,k]<- rbind(framesFirstVectors_G1[i,,k],
                             framesSecondVectors_G1[i,,k],
                             framesThirdVectors_G1[i,,k]) 
  }
}
frames_G2<-array(NA,dim = c(3,3,numberOfSkelPoints,nSamplesG2))
for (k in 1:nSamplesG2) {
  for (i in 1:numberOfSkelPoints) {
    frames_G2[,,i,k]<- rbind(framesFirstVectors_G2[i,,k],
                             framesSecondVectors_G2[i,,k],
                             framesThirdVectors_G2[i,,k]) 
  }
}

# children frames coordinates based on their parents frames
framesBasedOnParents_G1<-array(NA,dim = c(3,3,numberOfSkelPoints,nSamplesG1))
for (k in 1:nSamplesG1) {
  for (i in frameIndices) {
    k1<-framesCenters[i]
    k2<-parentsIndices[i]
    framesBasedOnParents_G1[,,k1,k]<-rotateFrameToMainAxes(myFrame = frames_G1[,,k2,k],
                                                           vectors2rotate = frames_G1[,,k1,k])
  } 
}
framesBasedOnParents_G2<-array(NA,dim = c(3,3,numberOfSkelPoints,nSamplesG2))
for (k in 1:nSamplesG2) {
  for (i in frameIndices) {
    k1<-framesCenters[i]
    k2<-parentsIndices[i]
    framesBasedOnParents_G2[,,k1,k]<-rotateFrameToMainAxes(myFrame = frames_G2[,,k2,k],
                                                           vectors2rotate = frames_G2[,,k1,k])
  } 
}

# spokes
tipOfCuttedUpSpokes_G1<-array(NA,dim = c(numberOfSkelPoints,3,nSamplesG1))
for (i in 1:nSamplesG1) {
  tipOfCuttedUpSpokes_G1[,,i]<-LP_ds_rep_inflationSimulation_G1[[i]]$tipOfCuttedUpSpokes
}
tipOfCuttedUpSpokes_G2<-array(NA,dim = c(numberOfSkelPoints,3,nSamplesG2))
for (i in 1:nSamplesG2) {
  tipOfCuttedUpSpokes_G2[,,i]<-LP_ds_rep_inflationSimulation_G2[[i]]$tipOfCuttedUpSpokes
}

tipOfCuttedDownSpokes_G1<-array(NA,dim = c(numberOfSkelPoints,3,nSamplesG1))
for (i in 1:nSamplesG1) {
  tipOfCuttedDownSpokes_G1[,,i]<-LP_ds_rep_inflationSimulation_G1[[i]]$tipOfCuttedDownSpokes
}
tipOfCuttedDownSpokes_G2<-array(NA,dim = c(numberOfSkelPoints,3,nSamplesG2))
for (i in 1:nSamplesG2) {
  tipOfCuttedDownSpokes_G2[,,i]<-LP_ds_rep_inflationSimulation_G2[[i]]$tipOfCuttedDownSpokes
}

upSpokeLength_G1<-array(NA, dim=c(numberOfSkelPoints,nSamplesG1))
for (k in 1:nSamplesG1) {
  for (i in frameIndices) {
    upSpokeLength_G1[i,k]<-norm(tipOfCuttedUpSpokes_G1[i,,k]-medialPoints3D_G1[i,,k],type = "2")
  }
}
upSpokeLength_G2<-array(NA, dim=c(numberOfSkelPoints,nSamplesG2))
for (k in 1:nSamplesG2) {
  for (i in frameIndices) {
    upSpokeLength_G2[i,k]<-norm(tipOfCuttedUpSpokes_G2[i,,k]-medialPoints3D_G2[i,,k],type = "2")
  }
}

downSpokeLength_G1<-array(NA, dim=c(numberOfSkelPoints,nSamplesG1))
for (k in 1:nSamplesG1) {
  for (i in frameIndices) {
    downSpokeLength_G1[i,k]<-norm(tipOfCuttedDownSpokes_G1[i,,k]-medialPoints3D_G1[i,,k],type = "2")
  }
}
downSpokeLength_G2<-array(NA, dim=c(numberOfSkelPoints,nSamplesG2))
for (k in 1:nSamplesG2) {
  for (i in frameIndices) {
    downSpokeLength_G2[i,k]<-norm(tipOfCuttedDownSpokes_G2[i,,k]-medialPoints3D_G2[i,,k],type = "2")
  }
}

#spoke directions in GCS
upSpokeDirections_G1<-array(NA, dim=c(numberOfSkelPoints,3,nSamplesG1))
for (k in 1:nSamplesG1) {
  for (i in frameIndices) {
    upSpokeDirections_G1[i,,k]<-convertVec2unitVec(tipOfCuttedUpSpokes_G1[i,,k]-medialPoints3D_G1[i,,k])
  }
}
upSpokeDirections_G2<-array(NA, dim=c(numberOfSkelPoints,3,nSamplesG2))
for (k in 1:nSamplesG2) {
  for (i in frameIndices) {
    upSpokeDirections_G2[i,,k]<-convertVec2unitVec(tipOfCuttedUpSpokes_G2[i,,k]-medialPoints3D_G2[i,,k])
  }
}

downSpokeDirections_G1<-array(NA, dim=c(numberOfSkelPoints,3,nSamplesG1))
for (k in 1:nSamplesG1) {
  for (i in frameIndices) {
    downSpokeDirections_G1[i,,k]<-convertVec2unitVec(tipOfCuttedDownSpokes_G1[i,,k]-medialPoints3D_G1[i,,k])
  }
}
downSpokeDirections_G2<-array(NA, dim=c(numberOfSkelPoints,3,nSamplesG2))
for (k in 1:nSamplesG2) {
  for (i in frameIndices) {
    downSpokeDirections_G2[i,,k]<-convertVec2unitVec(tipOfCuttedDownSpokes_G2[i,,k]-medialPoints3D_G2[i,,k])
  }
}


#spoke directions based on frames
upSpokesDirectionsBasedOnFrames_G1<-array(NA,dim = dim(upSpokeDirections_G1))
for (k in 1:nSamplesG1) {
  for (i in frameIndices) {
    
    upSpokesDirectionsBasedOnFrames_G1[i,,k]<-rotateFrameToMainAxes(myFrame = frames_G1[,,i,k],
                                                                    vectors2rotate = upSpokeDirections_G1[i,,k])
  }
}
upSpokesDirectionsBasedOnFrames_G2<-array(NA,dim = dim(upSpokeDirections_G2))
for (k in 1:nSamplesG2) {
  for (i in frameIndices) {
    
    upSpokesDirectionsBasedOnFrames_G2[i,,k]<-rotateFrameToMainAxes(myFrame = frames_G2[,,i,k],
                                                                    vectors2rotate = upSpokeDirections_G2[i,,k])
  }
}

downSpokesDirectionsBasedOnFrames_G1<-array(NA,dim = dim(downSpokeDirections_G1))
for (k in 1:nSamplesG1) {
  for (i in frameIndices) {
    
    downSpokesDirectionsBasedOnFrames_G1[i,,k]<-rotateFrameToMainAxes(myFrame = frames_G1[,,i,k],
                                                                      vectors2rotate = downSpokeDirections_G1[i,,k])
  }
}
downSpokesDirectionsBasedOnFrames_G2<-array(NA,dim = dim(downSpokeDirections_G2))
for (k in 1:nSamplesG2) {
  for (i in frameIndices) {
    
    downSpokesDirectionsBasedOnFrames_G2[i,,k]<-rotateFrameToMainAxes(myFrame = frames_G2[,,i,k],
                                                                      vectors2rotate = downSpokeDirections_G2[i,,k])
  }
}

# connections
connections_G1<-array(NA,dim = c(numberOfSkelPoints,3,nSamplesG1))
connectionsLengths_G1<-array(NA,dim=c(numberOfSkelPoints,nSamplesG1))
connectionsBasedOnParentFrames_G1<-array(NA,dim = dim(connections_G1))
for (k in 1:nSamplesG1) {
  for (i in 1:length(pointsIndices)) {
    
    k1<-pointsIndices[i]
    k2<-parentsIndices[i]
    
    if(k2==k1){
      connections_G1[k1,,k]<-c(0,0,0)
      connectionsLengths_G1[k1,k]<-0
      connectionsBasedOnParentFrames_G1[k1,,k]<-c(0,0,0)
      next
    }
    
    tempVec<-medialPoints3D_G1[k1,,k]- medialPoints3D_G1[k2,,k]
    
    connectionsLengths_G1[k1,k]<-norm(tempVec,type = "2")
    
    if(norm(tempVec,type = "2")==0){
      connections_G1[k1,,k]<-c(0,0,0)
    }else{
      connections_G1[k1,,k]<-convertVec2unitVec(tempVec)
    }
    
    connectionsBasedOnParentFrames_G1[k1,,k]<-rotateFrameToMainAxes(myFrame = frames_G1[,,k2,k],
                                                                    vectors2rotate = connections_G1[k1,,k])
  }
}
connections_G2<-array(NA,dim = c(numberOfSkelPoints,3,nSamplesG2))
connectionsLengths_G2<-array(NA,dim=c(numberOfSkelPoints,nSamplesG2))
connectionsBasedOnParentFrames_G2<-array(NA,dim = dim(connections_G2))
for (k in 1:nSamplesG2) {
  for (i in 1:length(pointsIndices)) {
    
    k1<-pointsIndices[i]
    k2<-parentsIndices[i]
    if(k2==k1){
      connections_G2[k1,,k]<-c(0,0,0)
      connectionsLengths_G2[k1,k]<-0
      connectionsBasedOnParentFrames_G2[k1,,k]<-c(0,0,0)
      next
    }
    
    tempVec<-medialPoints3D_G2[k1,,k]- medialPoints3D_G2[k2,,k]
    
    connectionsLengths_G2[k1,k]<-norm(tempVec,type = "2")
    
    if(norm(tempVec,type = "2")==0){
      connections_G2[k1,,k]<-c(0,0,0)
    }else{
      connections_G2[k1,,k]<-convertVec2unitVec(tempVec)
    }
    
    connectionsBasedOnParentFrames_G2[k1,,k]<-rotateFrameToMainAxes(myFrame = frames_G2[,,k2,k],
                                                                    vectors2rotate = connections_G2[k1,,k])
  }
}


# calculate LP sizes
skeletal_CentroidIndex<-which(pointsIndices==parentsIndices)

#NB! length of the center frame is 0 and must be excluded
LP_sizes_G1<-rep(NA,nSamplesG1)
for (i in 1:nSamplesG1) {
  LP_sizes_G1[i]<-exp(mean(c(log(connectionsLengths_G1[-skeletal_CentroidIndex,i]),
                             log(upSpokeLength_G1[frameIndices,i]),
                             log(downSpokeLength_G1[frameIndices,i]))))
}
LP_sizes_G2<-rep(NA,nSamplesG2)
for (i in 1:nSamplesG2) {
  LP_sizes_G2[i]<-exp(mean(c(log(connectionsLengths_G2[-skeletal_CentroidIndex,i]),
                             log(upSpokeLength_G2[frameIndices,i]),
                             log(downSpokeLength_G2[frameIndices,i]))))
}


#Removing or preserving the scale by LP-size based on the type of study

if(typeOfStudy=="sizeAndShapeAnalysis"){
  sizes_G1<-rep(1,nSamplesG1)
  sizes_G2<-rep(1,nSamplesG2)
  
  #we don't have scaling in size-and-shape analysis
  
  upSpokeLengthScaled_G1<-upSpokeLength_G1
  upSpokeLengthScaled_G2<-upSpokeLength_G2
  
  downSpokeLengthScaled_G1<-downSpokeLength_G1
  downSpokeLengthScaled_G2<-downSpokeLength_G2
  
  connectionsLengthsScaled_G1<-connectionsLengths_G1
  connectionsLengthsScaled_G2<-connectionsLengths_G2
  
}else if(typeOfStudy=="shapeAnalysis"){
  
  sizes_G1<-rep(NA,nSamplesG1)
  for (i in 1:nSamplesG1) {
    sizes_G1[i]<-exp(mean(c(log(connectionsLengths_G1[-skeletal_CentroidIndex,i]),
                            log(upSpokeLength_G1[frameIndices,i]),
                            log(downSpokeLength_G1[frameIndices,i]))))
  }
  sizes_G2<-rep(NA,nSamplesG2)
  for (i in 1:nSamplesG2) {
    sizes_G2[i]<-exp(mean(c(log(connectionsLengths_G2[-skeletal_CentroidIndex,i]),
                            log(upSpokeLength_G2[frameIndices,i]),
                            log(downSpokeLength_G2[frameIndices,i]))))
  }
  
  upSpokeLengthScaled_G1<-sweep(upSpokeLength_G1, 2, sizes_G1, "/") #2 indicate operation "/" on columns
  upSpokeLengthScaled_G2<-sweep(upSpokeLength_G2, 2, sizes_G2, "/") 
  
  downSpokeLengthScaled_G1<-sweep(downSpokeLength_G1, 2, sizes_G1, "/") #2 indicate operation "/" on columns
  downSpokeLengthScaled_G2<-sweep(downSpokeLength_G2, 2, sizes_G2, "/")
  
  connectionsLengthsScaled_G1<-sweep(connectionsLengths_G1, 2, sizes_G1, "/") #2 indicate operation "/" on columns
  connectionsLengthsScaled_G2<-sweep(connectionsLengths_G2, 2, sizes_G2, "/") 
  
}


# Calculate mean LP-dss-rep for visualization
framesBasedOnParentsVectorized_G1<-array(NA,dim = c(numberOfSkelPoints,9,nSamplesG1))
for (i in 1:nSamplesG1) {
  for (k in frameIndices) {
    framesBasedOnParentsVectorized_G1[k,,i]<-as.vector(t(framesBasedOnParents_G1[,,k,i]))
  }
}
framesBasedOnParentsVectorized_G2<-array(NA,dim = c(numberOfSkelPoints,9,nSamplesG2))
for (i in 1:nSamplesG2) {
  for (k in frameIndices) {
    framesBasedOnParentsVectorized_G2[k,,i]<-as.vector(t(framesBasedOnParents_G2[,,k,i]))
  }
}

meanFramesBasedOnParents_G1<-array(NA, dim = c(3,3,numberOfSkelPoints))
for (k in frameIndices) {
  if(k==skeletal_CentroidIndex){
    meanFramesBasedOnParents_G1[,,k]<-rbind(c(0,0,1),c(1,0,0),c(0,1,0))
  }else{
    tempVec<-mean(as.SO3(t(framesBasedOnParentsVectorized_G1[k,,])),type = 'geometric')
    meanFramesBasedOnParents_G1[,,k]<-matrix(tempVec,nrow = 3,byrow = TRUE)
  }
}
meanFramesBasedOnParents_G2<-array(NA, dim = c(3,3,numberOfSkelPoints))
for (k in frameIndices) {
  if(k==skeletal_CentroidIndex){
    meanFramesBasedOnParents_G2[,,k]<-rbind(c(0,0,1),c(1,0,0),c(0,1,0))
  }else{
    tempVec<-mean(as.SO3(t(framesBasedOnParentsVectorized_G2[k,,])),type = 'geometric')
    meanFramesBasedOnParents_G2[,,k]<-matrix(tempVec,nrow = 3,byrow = TRUE)
  }
}

meanFramesGlobalCoordinate_G1<-array(NA,dim = dim(meanFramesBasedOnParents_G1))
meanFramesGlobalCoordinate_G1[,,skeletal_CentroidIndex]<-rbind(c(0,0,1),c(1,0,0),c(0,1,0))
parents_Temp<-skeletal_CentroidIndex
while (length(parents_Temp)>0) {
  
  allChildrenIndicesTemp<-c()
  for (i in parents_Temp) {
    
    updatedParent<-meanFramesGlobalCoordinate_G1[,,i]
    
    childrenIndices<-na.omit(framesCenters[which(parentsIndices==i)])
    
    allChildrenIndicesTemp<-c(allChildrenIndicesTemp,childrenIndices)
    
    for (child_Index in childrenIndices) {
      meanFramesGlobalCoordinate_G1[,,child_Index]<-
        rotateFrameToMainAxesAndRotateBack(myFrame = updatedParent,
                                           vectorsInMainAxes = meanFramesBasedOnParents_G1[,,child_Index]) 
    }
    
  }
  allChildrenIndicesTemp<-unique(allChildrenIndicesTemp)
  
  parents_Temp<-na.omit(allChildrenIndicesTemp[!allChildrenIndicesTemp %in% parents_Temp])
  
}
meanFramesGlobalCoordinate_G2<-array(NA,dim = dim(meanFramesBasedOnParents_G2))
meanFramesGlobalCoordinate_G2[,,skeletal_CentroidIndex]<-rbind(c(0,0,1),c(1,0,0),c(0,1,0))
parents_Temp<-skeletal_CentroidIndex
while (length(parents_Temp)>0) {
  
  allChildrenIndicesTemp<-c()
  for (i in parents_Temp) {
    
    updatedParent<-meanFramesGlobalCoordinate_G2[,,i]
    
    childrenIndices<-na.omit(framesCenters[which(parentsIndices==i)])
    
    allChildrenIndicesTemp<-c(allChildrenIndicesTemp,childrenIndices)
    
    for (child_Index in childrenIndices) {
      meanFramesGlobalCoordinate_G2[,,child_Index]<-
        rotateFrameToMainAxesAndRotateBack(myFrame = updatedParent,
                                           vectorsInMainAxes = meanFramesBasedOnParents_G2[,,child_Index]) 
    }
    
  }
  allChildrenIndicesTemp<-unique(allChildrenIndicesTemp)
  
  parents_Temp<-na.omit(allChildrenIndicesTemp[!allChildrenIndicesTemp %in% parents_Temp])
  
}

# Calculate mean spokes' directions based on frames
meanUpSpokesDirectionsBasedOnFrames_G1<-array(NA,dim = c(numberOfSkelPoints,3))
for (i in frameIndices) {
  # For extremely concentrated data we use Mardia mean direction 
  pcaTemp<-prcomp(t(upSpokesDirectionsBasedOnFrames_G1[i,,]))
  if(pcaTemp$sdev[1]<1e-02 | pcaTemp$sdev[2]<1e-02){
    meanUpSpokesDirectionsBasedOnFrames_G1[i,]<-convertVec2unitVec(colMeans(t(upSpokesDirectionsBasedOnFrames_G1[i,,])))
  }else if(typeOfMeanDirection=="Frechet"){
    meanUpSpokesDirectionsBasedOnFrames_G1[i,]<-frechetMean(upSpokesDirectionsBasedOnFrames_G1[i,,]) 
  }else if(typeOfMeanDirection=="PNS"){
    sphereType<-kurtosisTestFunction(upSpokesDirectionsBasedOnFrames_G1[i,,])
    meanUpSpokesDirectionsBasedOnFrames_G1[i,]<-pns(upSpokesDirectionsBasedOnFrames_G1[i,,],sphere.type = sphereType)$PNS$mean 
  }else{
    stop("Please specify the typeOfMeanDirection by PNS or Frechet!")
  }
}
meanUpSpokesDirectionsBasedOnFrames_G2<-array(NA,dim = c(numberOfSkelPoints,3))
for (i in frameIndices) {
  # For extremely concentrated data we use Mardia mean direction 
  pcaTemp<-prcomp(t(upSpokesDirectionsBasedOnFrames_G2[i,,]))
  if(pcaTemp$sdev[1]<1e-02 | pcaTemp$sdev[2]<1e-02){
    meanUpSpokesDirectionsBasedOnFrames_G2[i,]<-convertVec2unitVec(colMeans(t(upSpokesDirectionsBasedOnFrames_G2[i,,])))
  }else if(typeOfMeanDirection=="Frechet"){
    meanUpSpokesDirectionsBasedOnFrames_G2[i,]<-frechetMean(upSpokesDirectionsBasedOnFrames_G2[i,,]) 
  }else if(typeOfMeanDirection=="PNS"){
    sphereType<-kurtosisTestFunction(upSpokesDirectionsBasedOnFrames_G2[i,,])
    meanUpSpokesDirectionsBasedOnFrames_G2[i,]<-pns(upSpokesDirectionsBasedOnFrames_G2[i,,],sphere.type = sphereType)$PNS$mean 
  }else{
    stop("Please specify the typeOfMeanDirection by PNS or Frechet!")
  }
}

# Calculate mean spokes' directions based on global coordinate (using mean frames in global coordinate)
meanUpSpokesDirectionsGlobalCoordinate_G1<-array(NA,dim = c(numberOfSkelPoints,3))
for (i in frameIndices) {
  meanUpSpokesDirectionsGlobalCoordinate_G1[i,]<-
    rotateFrameToMainAxesAndRotateBack(myFrame = meanFramesGlobalCoordinate_G1[,,i],
                                       vectorsInMainAxes = meanUpSpokesDirectionsBasedOnFrames_G1[i,])
  
}
meanUpSpokesDirectionsGlobalCoordinate_G2<-array(NA,dim = c(numberOfSkelPoints,3))
for (i in frameIndices) {
  meanUpSpokesDirectionsGlobalCoordinate_G2[i,]<-
    rotateFrameToMainAxesAndRotateBack(myFrame = meanFramesGlobalCoordinate_G2[,,i],
                                       vectorsInMainAxes = meanUpSpokesDirectionsBasedOnFrames_G2[i,])
  
}

#down spokes are along up spokes in opposite directions
meanDownSpokesDirectionsGlobalCoordinate_G1<-(-meanUpSpokesDirectionsGlobalCoordinate_G1)
meanDownSpokesDirectionsGlobalCoordinate_G2<-(-meanUpSpokesDirectionsGlobalCoordinate_G2)

# Calculate geometric mean of up spokes' lengths
meanUpSpokeLength_G1<-exp(rowMeans(log(upSpokeLengthScaled_G1)))
meanUpSpokeLength_G2<-exp(rowMeans(log(upSpokeLengthScaled_G2)))
meanDownSpokeLength_G1<-exp(rowMeans(log(downSpokeLengthScaled_G1)))
meanDownSpokeLength_G2<-exp(rowMeans(log(downSpokeLengthScaled_G2)))

# Calculate geometric mean of connections' lengths
meanConnectionsLengths_G1<-exp(rowMeans(log(connectionsLengthsScaled_G1)))
meanConnectionsLengths_G2<-exp(rowMeans(log(connectionsLengthsScaled_G2)))


# Calculate mean connection directions based on frames
meanConnectionsBasedOnParentFrames_G1<-array(NA,dim = c(numberOfSkelPoints,3))
for (i in 1:numberOfSkelPoints) {
  if(i==skeletal_CentroidIndex){
    meanConnectionsBasedOnParentFrames_G1[i,]<-c(0,0,0)
  }else{
    # For extremely concentrated data we use Mardia mean direction 
    pcaTemp<-prcomp(t(connectionsBasedOnParentFrames_G1[i,,]))
    if(pcaTemp$sdev[1]<1e-02 | pcaTemp$sdev[2]<1e-02){
      meanConnectionsBasedOnParentFrames_G1[i,]<-convertVec2unitVec(colMeans(t(connectionsBasedOnParentFrames_G1[i,,])))
    }else if(typeOfMeanDirection=="Frechet"){
      meanConnectionsBasedOnParentFrames_G1[i,]<-frechetMean(connectionsBasedOnParentFrames_G1[i,,]) 
    }else if(typeOfMeanDirection=="PNS"){
      sphereType<-kurtosisTestFunction(connectionsBasedOnParentFrames_G1[i,,])
      meanConnectionsBasedOnParentFrames_G1[i,]<-pns(connectionsBasedOnParentFrames_G1[i,,],sphere.type = sphereType)$PNS$mean  
    }else{
      stop("Please specify the typeOfMeanDirection by PNS or Frechet")
    }
  }
}
meanConnectionsBasedOnParentFrames_G2<-array(NA,dim = c(numberOfSkelPoints,3))
for (i in 1:numberOfSkelPoints) {
  if(i==skeletal_CentroidIndex){
    meanConnectionsBasedOnParentFrames_G2[i,]<-c(0,0,0)
  }else{
    # For extremely concentrated data we use Mardia mean direction 
    pcaTemp<-prcomp(t(connectionsBasedOnParentFrames_G2[i,,]))
    if(pcaTemp$sdev[1]<1e-02 | pcaTemp$sdev[2]<1e-02){
      meanConnectionsBasedOnParentFrames_G2[i,]<-convertVec2unitVec(colMeans(t(connectionsBasedOnParentFrames_G2[i,,])))
    }else if(typeOfMeanDirection=="Frechet"){
      meanConnectionsBasedOnParentFrames_G2[i,]<-frechetMean(connectionsBasedOnParentFrames_G2[i,,]) 
    }else if(typeOfMeanDirection=="PNS"){
      sphereType<-kurtosisTestFunction(connectionsBasedOnParentFrames_G2[i,,])
      meanConnectionsBasedOnParentFrames_G2[i,]<-pns(connectionsBasedOnParentFrames_G2[i,,],sphere.type = sphereType)$PNS$mean  
    }else{
      stop("Please specify the typeOfMeanDirection by PNS or Frechet")
    }
  }
}

# Calculate mean connection based on global coordinate (using mean frames in global coordinate)
meanConnectionsGlobalCoordinate_G1<-array(NA,dim = c(numberOfSkelPoints,3))
for (i in 1:numberOfSkelPoints) {
  k1<-pointsIndices[i]
  k2<-parentsIndices[i]
  meanConnectionsGlobalCoordinate_G1[k1,]<-
    rotateFrameToMainAxesAndRotateBack(myFrame = meanFramesGlobalCoordinate_G1[,,k2],
                                       vectorsInMainAxes = meanConnectionsBasedOnParentFrames_G1[k1,])
}
meanConnectionsGlobalCoordinate_G2<-array(NA,dim = c(numberOfSkelPoints,3))
for (i in 1:numberOfSkelPoints) {
  k1<-pointsIndices[i]
  k2<-parentsIndices[i]
  meanConnectionsGlobalCoordinate_G2[k1,]<-
    rotateFrameToMainAxesAndRotateBack(myFrame = meanFramesGlobalCoordinate_G2[,,k2],
                                       vectorsInMainAxes = meanConnectionsBasedOnParentFrames_G2[k1,])
}

# Convert mean LP-dss-rep to a GP-dss-rep for visualization

meanPositions_G1<-array(NA,dim = c(numberOfSkelPoints,3))
meanPositions_G1[skeletal_CentroidIndex,]<-c(0,0,0)
parents_Temp<-skeletal_CentroidIndex
while (length(parents_Temp)>0) {
  
  allChildrenIndicesTemp<-c()
  for (parent_Index in parents_Temp) {
    
    childrenIndices<-which(parentsIndices==parent_Index)
    
    allChildrenIndicesTemp<-c(allChildrenIndicesTemp,childrenIndices)
    
    for (child_Index in childrenIndices) {
      meanPositions_G1[child_Index,]<-meanPositions_G1[parent_Index,]+
        meanConnectionsLengths_G1[child_Index]*meanConnectionsGlobalCoordinate_G1[child_Index,]
    }
    
  }
  
  allChildrenIndicesTemp<-unique(allChildrenIndicesTemp)
  
  parents_Temp<-na.omit(allChildrenIndicesTemp[!allChildrenIndicesTemp %in% parents_Temp])
  
}
meanPositions_G2<-array(NA,dim = c(numberOfSkelPoints,3))
meanPositions_G2[skeletal_CentroidIndex,]<-c(0,0,0)
parents_Temp<-skeletal_CentroidIndex
while (length(parents_Temp)>0) {
  
  allChildrenIndicesTemp<-c()
  for (parent_Index in parents_Temp) {
    
    childrenIndices<-which(parentsIndices==parent_Index)
    
    allChildrenIndicesTemp<-c(allChildrenIndicesTemp,childrenIndices)
    
    for (child_Index in childrenIndices) {
      meanPositions_G2[child_Index,]<-meanPositions_G2[parent_Index,]+
        meanConnectionsLengths_G2[child_Index]*meanConnectionsGlobalCoordinate_G2[child_Index,]
    }
    
  }
  
  allChildrenIndicesTemp<-unique(allChildrenIndicesTemp)
  
  parents_Temp<-na.omit(allChildrenIndicesTemp[!allChildrenIndicesTemp %in% parents_Temp])
  
}


# mean spokes in GCS
meanUpSpokesTips_G1<-array(NA,dim = c(numberOfSkelPoints,3))
for (i in frameIndices) {
  meanUpSpokesTips_G1[i,]<-meanPositions_G1[i,]+
    meanUpSpokesDirectionsGlobalCoordinate_G1[i,]*meanUpSpokeLength_G1[i]
}
meanDownSpokesTips_G1<-array(NA,dim = c(numberOfSkelPoints,3))
for (i in frameIndices) {
  meanDownSpokesTips_G1[i,]<-meanPositions_G1[i,]+
    meanDownSpokesDirectionsGlobalCoordinate_G1[i,]*meanDownSpokeLength_G1[i]
}

meanUpSpokesTips_G2<-array(NA,dim = c(numberOfSkelPoints,3))
for (i in frameIndices) {
  meanUpSpokesTips_G2[i,]<-meanPositions_G2[i,]+
    meanUpSpokesDirectionsGlobalCoordinate_G2[i,]*meanUpSpokeLength_G2[i]
}
meanDownSpokesTips_G2<-array(NA,dim = c(numberOfSkelPoints,3))
for (i in frameIndices) {
  meanDownSpokesTips_G2[i,]<-meanPositions_G2[i,]+
    meanDownSpokesDirectionsGlobalCoordinate_G2[i,]*meanDownSpokeLength_G2[i]
}

# Plot overlaid LP-ds-rep means of PD and CG
#spine
numberOf2DspokePoints<-4
gorupIndicesOfVeins<-list()
k<-1
i<-1
while (k<(numberOfSkelPoints-2)) {
  gorupIndicesOfVeins[[i]]<-k:(k+(numberOf2DspokePoints+2))
  k<-k+(2*numberOf2DspokePoints-1)
  i<-i+1
}

spineIndices<-c()
for (i in 1:length(gorupIndicesOfVeins)) {
  spineIndices<-c(spineIndices,gorupIndicesOfVeins[[i]][1])
}
spineIndices<-c(numberOfSkelPoints-1,spineIndices,numberOfSkelPoints)

# plot mean LP-dss-rep ellipsoids vs mean LP-dss-rep inflated ellipsoids
if(TRUE){
  open3d()
  #connections
  for (j in 1:nrow(meanPositions_G1)) {
    k1<-pointsIndices[j]
    k2<-parentsIndices[j]
    if(k2==k1){
      next
    }
    vectors3d(meanPositions_G1[k1,],origin = meanPositions_G1[k2,],
              headlength = 0.2,radius = 1/10, col="blue", lwd=1)
  }
  #spokes
  for (i in frameIndices) {
    plot3d(rbind(meanPositions_G1[i,],meanUpSpokesTips_G1[i,]),type="l",col = "blue",expand = 10,box=FALSE,add = TRUE)
  }
  for (i in frameIndices) {
    plot3d(rbind(meanPositions_G1[i,],meanDownSpokesTips_G1[i,]),type="l",col = "red",expand = 10,box=FALSE,add = TRUE)
  }
  #spine
  plot3d(meanPositions_G1[spineIndices,],type="l",lwd=4,col = "blue",expand = 10,box=FALSE,add = TRUE)
  # # frames
  # for (i in frameIndices) {
  #   vectors3d(meanPositions_G1[i,]+meanFramesGlobalCoordinate_G1[1,,i],origin = meanPositions_G1[i,],headlength = 0.1,radius = 1/10, col="orange", lwd=1)
  #   vectors3d(meanPositions_G1[i,]+meanFramesGlobalCoordinate_G1[2,,i],origin = meanPositions_G1[i,],headlength = 0.1,radius = 1/10, col="orange", lwd=1)
  #   vectors3d(meanPositions_G1[i,]+meanFramesGlobalCoordinate_G1[3,,i],origin = meanPositions_G1[i,],headlength = 0.1,radius = 1/10, col="orange", lwd=1)
  # }
  
  # plot mean LP-ds-rep G2
  #connections
  for (j in 1:nrow(meanPositions_G2)) {
    k1<-pointsIndices[j]
    k2<-parentsIndices[j]
    if(k2==k1){
      next
    }
    vectors3d(meanPositions_G2[k1,],origin = meanPositions_G2[k2,],
              headlength = 0.2,radius = 1/10, col="red", lwd=1)
  }
  #spokes
  for (i in frameIndices) {
    plot3d(rbind(meanPositions_G2[i,],meanUpSpokesTips_G2[i,]),type="l",col = "blue",expand = 10,box=FALSE,add = TRUE)
  }
  for (i in frameIndices) {
    plot3d(rbind(meanPositions_G2[i,],meanDownSpokesTips_G2[i,]),type="l",col = "red",expand = 10,box=FALSE,add = TRUE)
  }
  #spine
  plot3d(meanPositions_G2[spineIndices,],type="l",lwd=4,col = "red",expand = 10,box=FALSE,add = TRUE)
  # # frames
  # for (i in frameIndices) {
  #   
  #   vectors3d(meanPositions_G2[i,]+meanFramesGlobalCoordinate_G2[1,,i],origin = meanPositions_G2[i,],headlength = 0.1,radius = 1/10, col="orange", lwd=1)
  #   vectors3d(meanPositions_G2[i,]+meanFramesGlobalCoordinate_G2[2,,i],origin = meanPositions_G2[i,],headlength = 0.1,radius = 1/10, col="orange", lwd=1)
  #   vectors3d(meanPositions_G2[i,]+meanFramesGlobalCoordinate_G2[3,,i],origin = meanPositions_G2[i,],headlength = 0.1,radius = 1/10, col="orange", lwd=1)
  # }
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE, main = "Mean shapes",
             sub = NULL, top = T, aspect = FALSE, expand = 1.1)
  
}

#####################################################################################################
#####################################################################################################
# LP-dss-rep Hypothesis testing

# hypothesis test on LP size
pValues_LP_sizes<-meanDifferenceTest1D(log(LP_sizes_G1),log(LP_sizes_G2),type = typeOfTest) 
cat("pValue of LP sizes is:",pValues_LP_sizes,"\n")
boxplot(LP_sizes_G1, LP_sizes_G2, names = c("CG","PD"),main="LP-size")
cat("Mean LP size G1:",mean(LP_sizes_G1),"mean LP size G2:",mean(LP_sizes_G2),"\n")
cat("sd LP size G1:",sd(LP_sizes_G1),"sd LP size G2:",sd(LP_sizes_G2),"\n")

# hypothesis thresholds to ignore extremely concentrated data
if(TRUE){
  thresholdDirections<-pi/50
  thresholdDirections
  thresholdLengths<-mean(c(LP_sizes_G1,LP_sizes_G2))/200
  thresholdLengths  
}

# hypothesis test on spokes' lengths
pValues_TtestUpSpokesLengths<-rep(NA,numberOfSkelPoints)
pb <- txtProgressBar(min = 0, max = numberOfSkelPoints, style = 3) #progress bar
for (i in frameIndices) {
  setTxtProgressBar(pb, i) #create progress bar
  
  if(abs(mean(upSpokeLengthScaled_G1[i,])-mean(upSpokeLengthScaled_G2[i,]))<thresholdLengths){
    pValues_TtestUpSpokesLengths[i]<-1
  }else{
    pValues_TtestUpSpokesLengths[i]<-meanDifferenceTest1D(log(upSpokeLengthScaled_G1[i,]),
                                                          log(upSpokeLengthScaled_G2[i,]),
                                                          type = typeOfTest)   
  }
  
}
sigPvalIndicesTemp<-which(!is.na(pValues_TtestUpSpokesLengths) & pValues_TtestUpSpokesLengths<=0.05)
# for (i in sigPvalIndicesTemp) {
#   boxplot(upSpokeLengthScaled_G1[i,], upSpokeLengthScaled_G2[i,], names = c("G1","G2"),main="up Spokes Lengths")
# }

pValues_TtestDownSpokesLengths<-rep(NA,numberOfSkelPoints)
pb <- txtProgressBar(min = 0, max = numberOfSkelPoints, style = 3) #progress bar
for (i in frameIndices) {
  setTxtProgressBar(pb, i) #create progress bar
  
  if(abs(mean(downSpokeLengthScaled_G1[i,])-mean(downSpokeLengthScaled_G2[i,]))<thresholdLengths){
    pValues_TtestUpSpokesLengths[i]<-1
  }else{
    pValues_TtestDownSpokesLengths[i]<-meanDifferenceTest1D(log(downSpokeLengthScaled_G1[i,]),
                                                            log(downSpokeLengthScaled_G2[i,]),
                                                            type = typeOfTest)
  }
}
sigPvalIndicesTemp<-which(!is.na(pValues_TtestDownSpokesLengths) & pValues_TtestDownSpokesLengths<=0.05)
# for (i in sigPvalIndicesTemp) {
#   boxplot(downSpokeLengthScaled_G1[i,], downSpokeLengthScaled_G2[i,], names = c("G1","G2"),main="down Spokes Lengths")
# }

# hypothesis test on connections' length
pValues_TtestConnectionsLengths<-rep(NA,numberOfSkelPoints)
pb <- txtProgressBar(min = 0, max = numberOfSkelPoints, style = 3) #progress bar
for (i in 1:numberOfSkelPoints) {
  setTxtProgressBar(pb, i) #create progress bar
  if(i==skeletal_CentroidIndex){
    pValues_TtestConnectionsLengths[i]<-1
  }else if(abs(mean(connectionsLengthsScaled_G1[i,])-mean(connectionsLengthsScaled_G2[i,]))<thresholdLengths){
    
    pValues_TtestConnectionsLengths[i]<-1
    
  }else{
    pValues_TtestConnectionsLengths[i]<-meanDifferenceTest1D(log(connectionsLengthsScaled_G1[i,]),
                                                             log(connectionsLengthsScaled_G2[i,]),
                                                             type = typeOfTest)
  }
}
sigPvalIndicesTemp<-which(pValues_TtestConnectionsLengths<=0.05)
# for (i in sigPvalIndicesTemp) {
#   boxplot(connectionsLengthsScaled_G1[i,], connectionsLengthsScaled_G2[i,], names = c("G1","G2"),main="connections Lengths")
# }

# we skip this part as spoke directions are the same as frames first elements!!!!

# # hypothesis test on spokes' directions based on local frames
# euclideanizedUpSpokesDirBasedOnFramesG1<-array(NA,dim = c(nSamplesG1,2,numberOfSkelPoints))
# euclideanizedUpSpokesDirBasedOnFramesG2<-array(NA,dim = c(nSamplesG2,2,numberOfSkelPoints))
# pValspokesDirectionsBasedOnFrames<-rep(NA,numberOfSkelPoints)
# pb <- txtProgressBar(min = 0, max = numberOfSkelPoints, style = 3) #progress bar
# for(i in frameIndices){
#   setTxtProgressBar(pb, i) #create progress bar
#   #NB! euclideanization must contain two groups because it uses the pooled mean 
#   euclideanizedTemp<-euclideanization(upSpokesDirectionsBasedOnFrames_G1[i,,],
#                                       upSpokesDirectionsBasedOnFrames_G2[i,,],
#                                       type = typeOfStudy4directions)
#   
#   pValUpSpokesDirectionsBasedOnFrames[i]<-meanDifferenceTestMultivariate(euclideanizedTemp$euclideanG1,
#                                                                        euclideanizedTemp$euclideanG2,
#                                                                        type=typeOfTest)
#   
#   euclideanizedUpSpokesDirBasedOnFramesG1[,,i]<-euclideanizedTemp$euclideanG1
#   euclideanizedUpSpokesDirBasedOnFramesG2[,,i]<-euclideanizedTemp$euclideanG2
#   
# }
# # which(pValspokesDirectionsBasedOnFrames<=0.05)


# hypothesis test on connections' directions based on local frames
euclideanizedConnectionsBasedOnParentFramesG1<-array(0,dim = c(nSamplesG1,2,numberOfSkelPoints))
euclideanizedConnectionsBasedOnParentFramesG2<-array(0,dim = c(nSamplesG2,2,numberOfSkelPoints))
pValConnectionsBasedOnParentFrames<-rep(NA,numberOfSkelPoints)
pb <- txtProgressBar(min = 0, max = numberOfSkelPoints, style = 3) #progress bar
for(i in 1:numberOfSkelPoints){
  setTxtProgressBar(pb, i) #create progress bar
  if(i==skeletal_CentroidIndex){
    pValConnectionsBasedOnParentFrames[i]<-1
    next
  }
  
  
  euclideanizedTemp<-euclideanization(connectionsBasedOnParentFrames_G1[i,,],
                                      connectionsBasedOnParentFrames_G2[i,,],
                                      type = typeOfStudy4directions)
  
  if(norm(colMeans(euclideanizedTemp$euclideanG1)-colMeans(euclideanizedTemp$euclideanG2),type = '2')<thresholdDirections){
    pValConnectionsBasedOnParentFrames[i]<-1
  }else{
    pValConnectionsBasedOnParentFrames[i]<-meanDifferenceTestMultivariate(euclideanizedTemp$euclideanG1,
                                                                          euclideanizedTemp$euclideanG2,
                                                                          type=typeOfTest) 
  }
  
  euclideanizedConnectionsBasedOnParentFramesG1[,,i]<-euclideanizedTemp$euclideanG1
  euclideanizedConnectionsBasedOnParentFramesG2[,,i]<-euclideanizedTemp$euclideanG2
  
}
# plot(euclideanizedConnectionsBasedOnParentFramesG1[,,25],col='blue',xlim = c(-pi,pi),ylim =  c(-0.5,0.5))
# par(new=TRUE)
# plot(euclideanizedConnectionsBasedOnParentFramesG2[,,25],col='red',xlim = c(-0.5,0.5),ylim =  c(-0.5,0.5))
# which(pValConnectionsBasedOnParentFrames<=0.05)

framesBasedOnParentsVectorized_G1<-array(NA,dim = c(numberOfSkelPoints,9,nSamplesG1))
for (i in 1:nSamplesG1) {
  for (k in 1:numberOfSkelPoints) {
    framesBasedOnParentsVectorized_G1[k,,i]<-as.vector(t(framesBasedOnParents_G1[,,k,i]))
  }
}
framesBasedOnParentsVectorized_G2<-array(NA,dim = c(numberOfSkelPoints,9,nSamplesG2))
for (i in 1:nSamplesG2) {
  for (k in 1:numberOfSkelPoints) {
    framesBasedOnParentsVectorized_G2[k,,i]<-as.vector(t(framesBasedOnParents_G2[,,k,i]))
  }
}

# hypothesis test on frames' normal directions based on parent frames
euclideanizedFrameBasedOnParentG1<-array(0,dim = c(nSamplesG1,3,numberOfSkelPoints))
euclideanizedFrameBasedOnParentG2<-array(0,dim = c(nSamplesG2,3,numberOfSkelPoints))
pValFramesBasedOnParent<-rep(NA,numberOfSkelPoints)
pb <- txtProgressBar(min = 0, max = numberOfSkelPoints, style = 3) #progress bar
for(i in frameIndices){
  setTxtProgressBar(pb, i) #create progress bar
  if(i==skeletal_CentroidIndex){
    pValFramesBasedOnParent[i]<-1
    next
  }
  
  Q4Temp<-as.Q4(as.SO3(t(framesBasedOnParentsVectorized_G1[i,,])))
  Q4Temp2<-matrix(as.numeric(t(Q4Temp)),ncol = 4,byrow = TRUE)
  for (j in 1:dim(Q4Temp2)[1]) {
    Q4Temp2[j,]<-Q4Temp2[j,]/norm(Q4Temp2[j,],type = '2')
  }
  Q4Temp3<-as.Q4(as.SO3(t(framesBasedOnParentsVectorized_G2[i,,])))
  Q4Temp4<-matrix(as.numeric(t(Q4Temp3)),ncol = 4,byrow = TRUE)
  for (j in 1:dim(Q4Temp4)[1]) {
    Q4Temp4[j,]<-Q4Temp4[j,]/norm(Q4Temp4[j,],type = '2')
  }
  
  euclideanizedTemp<-euclideanization(t(Q4Temp2),t(Q4Temp4),type = typeOfStudy4directions)
  
  
  if(norm(colMeans(euclideanizedTemp$euclideanG1)-colMeans(euclideanizedTemp$euclideanG2),type = '2')<thresholdDirections){
    pValFramesBasedOnParent[i]<-1
  }else{
    pValFramesBasedOnParent[i]<-meanDifferenceTestMultivariate(euclideanizedTemp$euclideanG1,
                                                               euclideanizedTemp$euclideanG2,
                                                               type=typeOfTest) 
  }
  
  euclideanizedFrameBasedOnParentG1[,,i]<-euclideanizedTemp$euclideanG1
  euclideanizedFrameBasedOnParentG2[,,i]<-euclideanizedTemp$euclideanG2
  
}
# for (i in frameIndices) {
#   plot3d(euclideanizedFrameBasedOnParentG1[,,i],type="p",col = "blue",expand = 10,box=FALSE,add = TRUE)
#   plot3d(euclideanizedFrameBasedOnParentG2[,,i],type="p",col = "red",expand = 10,box=FALSE,add = TRUE)
#   
# }
# which(!is.na(pValFramesBasedOnParent) & pValFramesBasedOnParent<=0.05)

# plot significant GOPs
pvalues_LP_ds_rep <- c(pValues_TtestUpSpokesLengths,
                       pValues_TtestDownSpokesLengths,
                       pValues_TtestConnectionsLengths,    
                       # pValUpSpokesDirectionsBasedOnFrames,  
                       pValConnectionsBasedOnParentFrames, 
                       pValFramesBasedOnParent,            
                       pValues_LP_sizes)

n_s<-numberOfSkelPoints
length(pvalues_LP_ds_rep)

alpha<-0.05
significantPvalues<-which(!is.na(pvalues_LP_ds_rep) & pvalues_LP_ds_rep<=alpha)
significantPvalues

#adjust p-values by Benjamini-Hochberg
FDR<-0.15
pvalues_LP_ds_rep_BH<-p.adjust(pvalues_LP_ds_rep,method = "BH")
pvalues_LP_ds_rep_Bonferroni<-p.adjust(pvalues_LP_ds_rep,method = "bonferroni")
significantPvalues_BH<-which(pvalues_LP_ds_rep_BH<=FDR)
significantPvalues_BH
significantPvalues_Bonferroni<-which(pvalues_LP_ds_rep_Bonferroni<=FDR)
significantPvalues_Bonferroni

cat("\n","Percentage of sig raw p-value is:",length(significantPvalues)/length(pvalues_LP_ds_rep),"\n")
cat("\n","Percentage of BH adjusted p-value is:",length(significantPvalues_BH)/length(pvalues_LP_ds_rep),"\n")



# plot by ggplot
df_LP <- data.frame(Type=c(rep("Raw p-value",sum(!is.na(pvalues_LP_ds_rep))),
                           rep("Bonferroni",sum(!is.na(pvalues_LP_ds_rep))),
                           rep("BH",sum(!is.na(pvalues_LP_ds_rep)))),
                    ordereOfPvalues=1:sum(!is.na(pvalues_LP_ds_rep)),
                    Values=c(sort(pvalues_LP_ds_rep[!is.na(pvalues_LP_ds_rep)]),
                             sort(pvalues_LP_ds_rep_Bonferroni[!is.na(pvalues_LP_ds_rep_Bonferroni)]),
                             sort(pvalues_LP_ds_rep_BH[!is.na(pvalues_LP_ds_rep_BH)])))
p<-ggplot(df_LP, aes(x=ordereOfPvalues, y=Values, group=Type))
p + geom_line(aes(linetype=Type),size=1)+
  # geom_line(aes(linetype=Type, col=Type),size=1)+
  # geom_point(aes(shape=Type),alpha=0.7)+
  geom_hline(yintercept=alpha,linetype="solid", color = "red")+
  geom_hline(yintercept=FDR,linetype="solid", color = "blue")+
  scale_linetype_manual(values=c("solid","dotdash", "dotted")) +
  theme_bw()+
  theme(plot.title = element_text(size = 17, hjust = 0.5),
        legend.text=element_text(size=17),
        legend.title=element_blank(),
        # panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        # legend.title=element_text(size=17),
        legend.position="bottom",
        axis.title=element_text(size=17),
        legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
  guides(colour = guide_legend(title.hjust = 0.5))+
  xlab("Ranking of p-values") + ylab("p-values")


#1
significantUpSpokesLengths<-significantPvalues[which(significantPvalues<=n_s)]
significantUpSpokesLengths
significantUpSpokesLengths_BH<-significantPvalues_BH[which(significantPvalues_BH<=n_s)]
significantUpSpokesLengths_BH
#2
significantDownSpokesLengths<-significantPvalues[which(n_s+1<=significantPvalues
                                                       & significantPvalues<=(2*n_s))]-n_s
significantDownSpokesLengths
significantDownSpokesLengths_BH<-significantPvalues_BH[which(n_s+1<=significantPvalues_BH
                                                             & significantPvalues_BH<=(2*n_s))]-n_s
significantDownSpokesLengths_BH
#2
significantConnectionsLengths<-significantPvalues[which((2*n_s+1)<=significantPvalues
                                                        & significantPvalues<=(3*n_s))]-(2*n_s)
significantConnectionsLengths
significantConnectionsLengths_BH<-significantPvalues_BH[which((2*n_s+1)<=significantPvalues_BH
                                                              & significantPvalues_BH<=(3*n_s))]-(2*n_s)
significantConnectionsLengths_BH
#3
significantConnectionsDirections<-significantPvalues[which((3*n_s+1)<=significantPvalues
                                                           & significantPvalues<=(4*n_s))]-(3*n_s)
significantConnectionsDirections
significantConnectionsDirections_BH<-significantPvalues_BH[which((3*n_s+1)<=significantPvalues_BH
                                                                 & significantPvalues_BH<=(4*n_s))]-(3*n_s)
significantConnectionsDirections_BH
#5
significantFrame<-significantPvalues[which((4*n_s+1)<=significantPvalues
                                           & significantPvalues<=(5*n_s))]-(4*n_s)
significantFrame
significantFrame_BH<-significantPvalues_BH[which((4*n_s+1)<=significantPvalues_BH
                                                 & significantPvalues_BH<=(5*n_s))]-(4*n_s)
significantFrame_BH


if(TRUE){
  #plot significant GOPs before and after the BH adjustment
  #1 spokes lengths
  open3d()
  # plot3d(meanPositions_G1,type="s",radius = 0.2,col = "blue",expand = 10,box=FALSE,add = TRUE)
  for (i in frameIndices) {
    plot3d(rbind(meanPositions_G1[i,],meanUpSpokesTips_G1[i,]),lwd = 1.5,type="l",col = "blue",expand = 10,box=FALSE,add = TRUE)
    plot3d(rbind(meanPositions_G1[i,],meanDownSpokesTips_G1[i,]),lwd = 1.5,type="l",col = "blue",expand = 10,box=FALSE,add = TRUE)
  }
  for (i in significantUpSpokesLengths) {
    plot3d(rbind(meanPositions_G1[i,],meanUpSpokesTips_G1[i,]),type="l",lwd = 7,col = "red",expand = 10,box=FALSE,add = TRUE)
  }
  for (i in significantDownSpokesLengths) {
    plot3d(rbind(meanPositions_G1[i,],meanDownSpokesTips_G1[i,]),type="l",lwd = 7,col = "red",expand = 10,box=FALSE,add = TRUE)
  }
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE,
             # main = "Significant spokes' lengths",
             sub = NULL, top = T, aspect = FALSE, expand = 1.1)
  # pp <- par3d(no.readonly=TRUE)
  # dput(pp, file="plotView.R", control = "all")
  # pp <- dget("plotView.R")
  # par3d(userMatrix=pp$userMatrix)
  
  #with correction
  open3d()
  # plot3d(meanPositions_G1,type="s",radius = 0.2,col = "blue",expand = 10,box=FALSE,add = TRUE)
  for (i in frameIndices) {
    plot3d(rbind(meanPositions_G1[i,],meanUpSpokesTips_G1[i,]),lwd = 1.5,type="l",col = "blue",expand = 10,box=FALSE,add = TRUE)
    plot3d(rbind(meanPositions_G1[i,],meanDownSpokesTips_G1[i,]),lwd = 1.5,type="l",col = "blue",expand = 10,box=FALSE,add = TRUE)
  }
  for (i in significantUpSpokesLengths_BH) {
    plot3d(rbind(meanPositions_G1[i,],meanUpSpokesTips_G1[i,]),type="l",lwd = 7,col = "red",expand = 10,box=FALSE,add = TRUE)
  }
  for (i in significantDownSpokesLengths_BH) {
    plot3d(rbind(meanPositions_G1[i,],meanDownSpokesTips_G1[i,]),type="l",lwd = 7,col = "red",expand = 10,box=FALSE,add = TRUE)
  }
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE,
             # main = "Significant spokes' lengths after BH adjustment",
             sub = NULL, top = T, aspect = FALSE, expand = 1.1)
  # pp <- dget("plotView.R")
  # par3d(userMatrix=pp$userMatrix)
  
  
  # spokes directions 
  open3d()
  # plot3d(meanPositions_G1,type="s",radius = 0.2,col = "blue",expand = 10,box=FALSE,add = TRUE)
  for (i in frameIndices) {
    plot3d(rbind(meanPositions_G1[i,],meanUpSpokesTips_G1[i,]),lwd = 1.5,type="l",col = "blue",expand = 10,box=FALSE,add = TRUE)
    plot3d(rbind(meanPositions_G1[i,],meanDownSpokesTips_G1[i,]),lwd = 1.5,type="l",col = "blue",expand = 10,box=FALSE,add = TRUE)
  }
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE,
             # main = "Significant spokes' directions",
             sub = NULL, top = T, aspect = FALSE, expand = 1.1)
  # pp <- dget("plotView.R")
  # par3d(userMatrix=pp$userMatrix)
  
  # with correction 
  open3d()
  # plot3d(meanPositions_G1,type="s",radius = 0.2,col = "blue",expand = 10,box=FALSE,add = TRUE)
  for (i in frameIndices) {
    plot3d(rbind(meanPositions_G1[i,],meanUpSpokesTips_G1[i,]),lwd = 1.5,type="l",col = "blue",expand = 10,box=FALSE,add = TRUE)
    plot3d(rbind(meanPositions_G1[i,],meanDownSpokesTips_G1[i,]),lwd = 1.5,type="l",col = "blue",expand = 10,box=FALSE,add = TRUE)
  }
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE,
             # main = "Significant spokes' directions after BH adjustment",
             sub = NULL, top = T, aspect = FALSE, expand = 1.1)
  # pp <- dget("plotView.R")
  # par3d(userMatrix=pp$userMatrix)
  
  #2 connections lengths
  open3d()
  for (j in 1:nrow(meanPositions_G1)) {
    k1<-pointsIndices[j]
    k2<-parentsIndices[j]
    if(k2==k1){
      next
    }
    vectors3d(meanPositions_G1[k1,],origin = meanPositions_G1[k2,],
              headlength = 0.2,radius = 1/6, col="blue", lwd=1)
  }
  for (j in significantConnectionsLengths) {
    k1<-pointsIndices[j]
    k2<-parentsIndices[j]
    if(k2==k1){
      next
    }
    vectors3d(meanPositions_G1[k1,],origin = meanPositions_G1[k2,],
              headlength = 0.2,radius = 1/6, col="red", lwd=5)
  }
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE,
             # main = "Significant connections' lengths",
             sub = NULL, top = T, aspect = FALSE, expand = 1.1)
  # pp <- dget("plotView.R")
  # par3d(userMatrix=pp$userMatrix)
  
  #with correction
  open3d()
  for (j in 1:nrow(meanPositions_G1)) {
    k1<-pointsIndices[j]
    k2<-parentsIndices[j]
    if(k2==k1){
      next
    }
    vectors3d(meanPositions_G1[k1,],origin = meanPositions_G1[k2,],
              headlength = 0.2,radius = 1/6, col="blue", lwd=1)
  }
  for (j in significantConnectionsLengths_BH) {
    k1<-pointsIndices[j]
    k2<-parentsIndices[j]
    if(k2==k1){
      next
    }
    vectors3d(meanPositions_G1[k1,],origin = meanPositions_G1[k2,],
              headlength = 0.2,radius = 1/6, col="red", lwd=5)
  }
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE,
             # main = "Significant connections' lengths after BH adjustment",
             sub = NULL, top = T, aspect = FALSE, expand = 1.1)
  # pp <- dget("plotView.R")
  # par3d(userMatrix=pp$userMatrix)
  
  #3 connections directions
  open3d()
  for (j in 1:nrow(meanPositions_G1)) {
    k1<-pointsIndices[j]
    k2<-parentsIndices[j]
    if(k2==k1){
      next
    }
    vectors3d(meanPositions_G1[k1,],origin = meanPositions_G1[k2,],
              headlength = 0.2,radius = 1/6, col="blue", lwd=1)
  }
  for (j in significantConnectionsDirections) {
    k1<-pointsIndices[j]
    k2<-parentsIndices[j]
    if(k2==k1){
      next
    }
    vectors3d(meanPositions_G1[k1,],origin = meanPositions_G1[k2,],
              headlength = 0.2,radius = 1/6, col="red", lwd=5)
  }
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE,
             # main = "Significant connections' directions",
             sub = NULL, top = T, aspect = FALSE, expand = 1.1)
  # pp <- dget("plotView.R")
  # par3d(userMatrix=pp$userMatrix)
  
  #with correction
  open3d()
  for (j in 1:nrow(meanPositions_G1)) {
    k1<-pointsIndices[j]
    k2<-parentsIndices[j]
    if(k2==k1){
      next
    }
    vectors3d(meanPositions_G1[k1,],origin = meanPositions_G1[k2,],
              headlength = 0.2,radius = 1/6, col="blue", lwd=1)
  }
  for (j in significantConnectionsDirections_BH) {
    k1<-pointsIndices[j]
    k2<-parentsIndices[j]
    if(k2==k1){
      next
    }
    vectors3d(meanPositions_G1[k1,],origin = meanPositions_G1[k2,],
              headlength = 0.2,radius = 1/6, col="red", lwd=5)
  }
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE,
             # main = "Significant connections' directions after BH adjustment",
             sub = NULL, top = T, aspect = FALSE, expand = 1.1)
  # pp <- dget("plotView.R")
  # par3d(userMatrix=pp$userMatrix)
  
  #5 frames
  open3d()
  vectors3d(meanPositions_G1[frameIndices,]+t(meanFramesGlobalCoordinate_G1[1,,frameIndices]),origin = meanPositions_G1[frameIndices,],headlength = 0.2,radius = 1/10, col="blue", lwd=1)
  vectors3d(meanPositions_G1[frameIndices,]+t(meanFramesGlobalCoordinate_G1[2,,frameIndices]),origin = meanPositions_G1[frameIndices,],headlength = 0.2,radius = 1/10, col="blue", lwd=1)
  vectors3d(meanPositions_G1[frameIndices,]+t(meanFramesGlobalCoordinate_G1[3,,frameIndices]),origin = meanPositions_G1[frameIndices,],headlength = 0.2,radius = 1/10, col="blue", lwd=1)
  for (i in significantFrame) {
    vectors3d(meanPositions_G1[i,]+meanFramesGlobalCoordinate_G1[1,,i],origin = meanPositions_G1[i,],headlength = 0.2,radius = 1/10, col="red", lwd=5)
    vectors3d(meanPositions_G1[i,]+meanFramesGlobalCoordinate_G1[2,,i],origin = meanPositions_G1[i,],headlength = 0.2,radius = 1/10, col="red", lwd=5)
    vectors3d(meanPositions_G1[i,]+meanFramesGlobalCoordinate_G1[3,,i],origin = meanPositions_G1[i,],headlength = 0.2,radius = 1/10, col="red", lwd=5)
  }
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE,
             # main = "Significant frames",
             sub = NULL, top = T, aspect = FALSE, expand = 1.1)
  # pp <- dget("plotView.R")
  # par3d(userMatrix=pp$userMatrix)
  
  #with correction
  open3d()
  vectors3d(meanPositions_G1[frameIndices,]+t(meanFramesGlobalCoordinate_G1[1,,frameIndices]),origin = meanPositions_G1[frameIndices,],headlength = 0.2,radius = 1/10, col="blue", lwd=1)
  vectors3d(meanPositions_G1[frameIndices,]+t(meanFramesGlobalCoordinate_G1[2,,frameIndices]),origin = meanPositions_G1[frameIndices,],headlength = 0.2,radius = 1/10, col="blue", lwd=1)
  vectors3d(meanPositions_G1[frameIndices,]+t(meanFramesGlobalCoordinate_G1[3,,frameIndices]),origin = meanPositions_G1[frameIndices,],headlength = 0.2,radius = 1/10, col="blue", lwd=1)
  for (i in significantFrame_BH) {
    vectors3d(meanPositions_G1[i,]+meanFramesGlobalCoordinate_G1[1,,i],origin = meanPositions_G1[i,],headlength = 0.2,radius = 1/10, col="red", lwd=5)
    vectors3d(meanPositions_G1[i,]+meanFramesGlobalCoordinate_G1[2,,i],origin = meanPositions_G1[i,],headlength = 0.2,radius = 1/10, col="red", lwd=5)
    vectors3d(meanPositions_G1[i,]+meanFramesGlobalCoordinate_G1[3,,i],origin = meanPositions_G1[i,],headlength = 0.2,radius = 1/10, col="red", lwd=5)
  }
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE,
             # main = "Significant frames after BH adjustment",
             sub = NULL, top = T, aspect = FALSE, expand = 1.1)
  # pp <- dget("plotView.R")
  # par3d(userMatrix=pp$userMatrix)
  
}



