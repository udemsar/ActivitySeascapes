# ActivitySeascapes
----------------------------------------------
R Code implementing the activity seascape calculation for trajectories 

Author: Urska Demsar, University of St Andrews, UK, urska.demsar@st-andrews.ac.uk

21 May 2018

----------------------------------------------
This code is free and open and provided as supplementary information to the following article:

Papastamatiou YP, Watanabe YY, Demšar U, Leos-Barajas V, Langrock R, Bradley D, Weng K, Lowe C, Friedlander A and Caselle J, 2018, Activity seascapes highlight central place refuging strategies in marine predators that never stop swimming. 
https://movementecologyjournal.biomedcentral.com/articles/10.1186/s40462-018-0127-3

The code is provided as is and can be changed, but we would appreciate if you would please, cite the paper if you use the code as well as send us a message - we will be delighted if anyone is interested in our work.

----------------------------------------------
Code uses and implements methods for stacked space-time densities and Brownian bridges. 

For stacked space-time density algorithm see:

Demsar U, Buchin K, van Loon EE and Shamoun-Baranes J, 2015, 
Stacked space-time densities: a geovisualisation approach to explore 
dynamics of space use over time. GeoInformatica, 19(1):85-115.

------------------------------------------------

How to run the code:

The main file is StackedSpaceTimeDensities_WithActivityProbability.R, which specifies values for input files, parameters and outputs. The file runs on example data which are provided as part of this repository. These data include test trajectories and a daily time series of activity probability.

For use on other data and with other input parameters, the top section of the file StackedSpaceTimeDensities_WithActivityProbability.R file ('Parameters to be set/changed by the user') should be appropriately edited.

