This repository is associated with the bellow publication.

## Article
Mohsen Taheri, Stephen M. Pizer, Jörn Schulz et al. Fitting the Discrete Swept Skeletal Representation to Slabular Objects, 18 May 2023, PREPRINT available at Research Square


### Link

[https://doi.org/10.21203/rs.3.rs-2927062/v1]

### Cite:
```
@article{taheri2023fitting,
  title={Fitting the Discrete Swept Skeletal Representation to Slabular Objects},
  author={Taheri, Mohsen and Pizer, Stephen M and Schulz, J{\"o}rn},
  year={2023}
}
```
************************************************************************
*************************** main.R description *************************
************************************************************************
There are two R.scripts as "main.R" and "hypothesisTest_Ellipsoids_vs_deformedEllipsoids.R"

Fitting Locally Parameterized Swept Skeletal Structure (LP-dss-rep) to slabular objects

fit_LP_dss_rep function is in the main.R script
The input values of fit_LP_dss_rep function is as follows

Arguments                     values               Details
--------------------------------------------------------------------------------------------------------------------------------------------
tmesh                          mesh3d object        3D triangle mesh object
plotting                       TRUE/FALSE           Plotting all the steps of the model fitting procedure
numberOfCoverPoints3D          integer			   Number of uniformly distributed points 
k_Ring3D                       integer              For dividing the 3D object
lambda3D                       real number          Penalization parameter for dividing the 3D object into the top and bottom parts
k_Ring2D                       integer              For dividing the 2D object
lambda2D                       real number          Penalization parameter for dividing the 3D object into the top and bottom parts
sphereResolution               1,2,or 3             3D urchins' resolution
circleResolution               integer              2D urchins' resolution
polyDegree3D                   integer              3D polynomial degree
polyDegree2D                   integer              2D polynomial degree
alpha1                         real number          Alpha value for convex hull
numberOf2DspokePoints          integer              Number of points on each vein
numberOfSpanialPoints          integer              Number of spinal points
cross_sections_visualization   TRUE/FALSE           For visualizing the cross-sections


************************************************************************
********** hypothesisTest_Ellipsoids_vs_deformedEllipsoids.R ***********
************************************************************************
LP-ds-rep hypothesis testing versus LP-dss-rep hypothesis testing to detect local dissimilarities 
between a group of ellipsoids versus a group of inflated ellipsoids 

Please run the R script line by line


![sweptSkeletalStructure](https://github.com/MohsenTaheriShalmani/LP-dss-rep/assets/19237855/4882edbf-8245-4994-9de4-d38b7daecb76)



