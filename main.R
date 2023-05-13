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
source("fit_LP_dss_rep.R")

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


fited_LP_ds_rep<-fit_LP_dss_rep(tmesh=tmesh,
                                plotting=TRUE,
                                numberOfCoverPoints3D=1000000,
                                k_Ring3D=10, # for dividing 2D and 3D object
                                lambda3D=0.5, # for dividing 2D and 3D object
                                k_Ring2D=5, # for dividing 2D and 3D object
                                lambda2D=2.5, # for dividing 2D and 3D object
                                sphereResolution=1, # choose 1,2,or 3 for the resolution of the urchin
                                circleResolution=24, # resolution of the 2D urchin
                                circleResolutionHigh=150, # resolution of the 2D urchin
                                urchinRadius=0.5, # resolution of the 2D urchin
                                thresholdDistance2D=0.2,
                                polyDegree3D=4,
                                polyDegree2D=2,
                                alpha1=2, # for alpha convex hull
                                numberOfPoints4alphaHull=5000,
                                numberOf2DspokePoints=4,
                                numberOfSpanialPoints=20,
                                numberOfSpokes4Interpolation=5,
                                rotationGap=10, #to fit best medial curve
                                cross_sections_visualization=FALSE)


