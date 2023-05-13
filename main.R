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
# functions 

source("subFunctions/MathFunctions.R")
source("subFunctions/euclideanization.R")
source("subFunctions/ray_triangle_intersection.R")
source("subFunctions/cutAndStrechSpokes.R")
source("subFunctions/forceFunctionAngleBased.R")
source("subFunctions/frameGenerator.R")
source("subFunctions/normalsOfSkeletalSheetBySpline.R")
source("subFunctions/normalsOfSkeletalSheetByConnections.R")
source("subFunctions/rotateFrameForwardAndBackward.R")
source("subFunctions/meanFrames.R")
source("fit_slabular_LP_ds_rep_CrossSectionsNormal2Spine.R")

#####################################################################################################
#####################################################################################################
# load mesh

meshPoints <- read.csv(file = paste("files/Mesh.csv",sep = ""),check.names = FALSE, header=TRUE, sep=",")

# connections of triangular mesh
PolygonsCsv <- read.csv("files/Mesh_Polygon.csv")
polyMatrix<-cbind((PolygonsCsv$point1+1),(PolygonsCsv$point2+1),(PolygonsCsv$point3+1))

spharmPDM_Mesh<-matrix(meshPoints[[1]], ncol = 3, byrow = TRUE)


if(TRUE){
  open3d()
  verts <- rbind(t(as.matrix(spharmPDM_Mesh)),1)
  trgls <- as.matrix(t(polyMatrix))
  tmesh <- tmesh3d(verts, trgls)
  wire3d(tmesh, col="grey",alpha=0.2)  #wire mesh
  shade3d(tmesh, col="white",alpha=0.2)  #surface mesh
}



