There are two R.scripts as "main.R" and "hypothesisTest_Ellipsoids_vs_deformedEllipsoids.R"


************************************************************************
*************************** main.R discription *************************
************************************************************************
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
LP-ds-rep hypothesis testing versus LP-dss-rep hypothesis testing to detect local dissimilariteis 
between a group of ellipsoids versus a group of inflated ellipsoids 

Please run the R script line by line
	
	
	
	
	
