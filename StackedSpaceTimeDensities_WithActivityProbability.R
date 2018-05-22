# ------------------------------------------------
# Activity seascapes based on stacked space-time densities around a set of trajectories
# ========================
# Code provided as is and can be used or modified freely. 
#
# Code belongs to the following paper:
#
# Papastamatiou YP, Watanabe YY, Demsar U, Leos-Barajas V, 
# Langrock R, Bradley D, Weng K, Lowe C, Friedlander A and 
# Caselle J, 2018, Activity seascapes highlight central place 
# refuging strategies in marine predators that never stop swimming. 
# Under Review, Movement Ecology. 
#
# Space-time densities algorithm is from this paper:
#
# Demsar U, Buchin K, van Loon EE and Shamoun-Baranes J, 2015, 
# Stacked space-time densities: a geovisualisation approach to explore 
# dynamics of space use over time. GeoInformatica, 19(1):85-115.
# DOI 10.1007/s10707-014-0207-5
# ------------------------------------------------
# Author: Urska Demsar
# University of St Andrews
# St Andrews, Scotland, UK
# http://udemsar.com
# ------------------------------------------------
# This is the main file that implements the activity seascape calculation
# ------------------------------------------------
# Supporting R files:
# DecayFunctions.R
# DensityAroundOnePoint.R
# DensityAroundTrajectory.R
# SquaredDistance2points.R
# WriteCSVtableFromFourRectangular3DArrays.R
# ------------------------------------------------
# Input 1: list of trajectories
# This must be a csv file including the following four variables:
# ID, xcoord, ycoord, zcoord, where
# - ID - ID of each trajectory
# - X - x coordinate of each trajectory point in some projected
# coordinate system (in metres)
# - Y - x coordinate of each trajectory point in some projected
# coordinate system (in metres)
# - T - time of the day
#
# Input 2: csv file with probability of activity at each second in a day
# ------------------------------------------------
# The following also need to be manually specified: 
# - volume extent (min/max X, min/max Y, min/max T)
# - voxel size
# - kernel size 
# - kernel type
# ------------------------------------------------
# Output 1: a csv file with density, activity probability and activity seascape
# given as volumes, i.e. with five columns:
# x, y, t, density value, activity seascape value (density x actProb)
#
# Output 2: the activity probability volume in a separate file
#
# Ouput 3: space-time density for each individual shark
#
# Note: you will need volumetric visualisation software to visualise these outputs,
# we recommend either Voxler (free trial available, 
# http://www.goldensoftware.com/products/voxler) or
# or Paraview (free and open source, https://www.paraview.org/)
# ------------------------------------------------


# ------------------------------------------------
# Parameters to be set/changed by the user 

# Set working directory
setwd('.')

# Set voxel size in metres 
voxelSize <- 100

# Set kernel size in metres
kernelSize <- 600

# Set spatio-temporal extent of the volume
# This is done manually and not by reading the min/max coordinates from the
# trajectories file, so the program can be used for densities calculated
# for e.g. different days, where each day the trajectories have different
# spatial extents, but we want to add them up at the end into one general 
# space-time cube in which all density volumes must have the same extent.

# X & Y axis - extent in metres for the text data set:
minXcoord <- 814200 # westernmost point
maxXcoord <- 832800 # easternmost point
minYcoord <- 647900 # southernmost point
maxYcoord <- 653200 # northernmost point

# Z axis - temporal extent, in seconds across 24hrs
minZcoord <- 0 # start point in time
maxZcoord <- 14400 # end point in time

# Select kernel type: 1-linear, 2-bisquare, 3-Gaussian, 4-Brownian bridges
#method <- 1
#method <- 2
#method <- 3
method <- 4

# If Brownian bridges, the user sets sigmas 1 and 2
if (method == 4) {
  sigma11 <- 1000 # for sharks
  sigma12 <- 300 # for sharks
  #sigma11 <- 3 # to be set by the user - this is just for testing
  #sigma12 <- 1 # to be set by the user - this is just for testing
} else {
  sigma11 <- 0
  sigma12 <- 0
}

# Name of input and output files

# Input
trajectories.file <- 'Shark_Trajectories_Test.csv' # Input 1 

# File with activity probability
activity.file <- 'ActivityProbability_Test.csv' # Input 2

# Output
output.file <- 'Shark_Density_Test.csv'

# ------------------------------------------------
# General setup & read data

# Read fuctions in additional files
# Functions to calculate density for each trajectory or around one point 
source("DensityAroundTrajectory.R")
source("DensityAroundOnePoint.R") 
# Functions for saving density as as csv coordinate file
source("WriteCSVtableFromRectangular3DArrays.R") 
# Functions for geometric and kernel calculations
source("SquaredDistance2points.R")
source("DecayFunctions.R")

# Package for 3D plotting & colour schemes
library(plot3D)
library(RColorBrewer)

# Read trajectory data
dfAll <- read.csv(trajectories.file,stringsAsFactors=FALSE)
head(dfAll)

# Read activity probability data
activProb <- read.csv(activity.file,stringsAsFactors=FALSE )
head(activProb)

# ------------------------------------------------
# Build three 3D arrays of x, y, z coordinates and initialise the total density volume

startX <- floor(minXcoord/voxelSize)*voxelSize
startY <- floor(minYcoord/voxelSize)*voxelSize
startZ <- floor(minZcoord/voxelSize)*voxelSize
endX <- ceiling(maxXcoord/voxelSize)*voxelSize
endY <- ceiling(maxYcoord/voxelSize)*voxelSize
endZ <- ceiling(maxZcoord/voxelSize)*voxelSize

xvoxels <- (endX-startX)/voxelSize
yvoxels <- (endY-startY)/voxelSize
zvoxels <- (endZ-startZ)/voxelSize

# Build the volume
x <- seq(startX,endX,by=voxelSize)
y <- seq(startY,endY,by=voxelSize)
z <- seq(startZ,endZ,by=voxelSize)
M <- mesh(x,y,z)
xcoord <- M$x
ycoord <- M$y
zcoord <- M$z

# Initialise the total density as zeros everywhere in a 3D array of the same 
# size as the xcoord/ycoord/vcoord volumes

totalDensity <- array(data=0,dim=dim(xcoord))

# ------------------------------------------------
# Build an activity probability volume to multiply the density with to get the 
# activity seascape

# Create activity probability volume

activityVolume <- array(data=0,dim=dim(xcoord))
actL <- dim(activProb)

for (i in 1:actL[1]) {
  # find locations in zcoord where values are the same as in activProb
  ind <- which(zcoord==activProb[i,1])
  # Change these locations in activityVolume into relevant probabilities
  activityVolume[ind] <- activProb[i,2]
} # for (i in 1:length(activProb))

# Save activity probability volume for visual display later
WriteCSVtableFromFourRectangular3DArrays(xcoord,ycoord,zcoord,activityVolume,'ActivityProbabilityVolume.csv')

# ------------------------------------------------
# Calculate daily space-time density for all sharks

listOfSharks <- unique(dfAll$ID_shark)
noOfSharks <- length(listOfSharks)

# Loop across all sharks
for (shk in 1:noOfSharks) {

  # Initialise individual shark density as zeros everywhere in a 3D array of the same 
  # size as the xcoord/ycoord/vcoord volumes
  sharkDensity <- array(data=0,dim=dim(xcoord))
  
  # find where this shark shk is in the table
  sharkIndices <- which(dfAll$ID_shark == listOfSharks[shk])

  # extract only trajectories of this one shark
  trajData <- dfAll[sharkIndices,]
  head(trajData)

  # ------------------------------------------------
  # Count the number of different trajectories for this shark 

  # listOfTrajIDs = list of unique ID values in first column of trajData 
  # listOfTrajLengths = how many points there are in each trajectory
  # trajNo = number of trajectories
  
  listOfTrajIDs <- unique(trajData$ID_day)
  trajNo <- length(listOfTrajIDs)
  listOfTrajLengths <- vector(mode="numeric",length=trajNo)

  for (i in 1:trajNo) {
    # find indices where this trajectory is in dfAll table
    trajIndices <- which(trajData$ID_day == listOfTrajIDs[i])
  
    # assign this to the vector of lengths
    listOfTrajLengths[i] <- length(trajIndices)
  }

  # Loop through trajData using the list of unique traj IDs
  # and at each step extract the trajectory and calculate the density.

  currentRow <- 1

  for (i in 1:trajNo) {
    # Find ID of the current trajectory
    currentTrajID <- listOfTrajIDs[i]
    currentTrajLength <- listOfTrajLengths[i]

    # currentLine will be a submatrix in trajData starting at currentRow with
    # currentTrajLength rows and of columns 2, 3, 4, which have coordinate
    # data
    currentLine <- trajData[currentRow:(currentRow+currentTrajLength-1),3:5]
    currentRow <- currentRow+currentTrajLength
  
    # Calculate density for current line and plot the line and the density
    print('*********')
    print(paste(sprintf('Shark no.: %d',shk), sprintf('out of %d', noOfSharks)))
    print(paste(sprintf('Its trajectory: %d',i), sprintf('out of %d',trajNo)))
    print(paste(sprintf('Density for trajectory no.: %d',currentTrajID),sprintf('with this many points: %d', currentTrajLength) ))

    # Note about data: as these are acoustic detections, some sharks may 
    # have only one point in a daily trajectory, meaning that there won't be a segment for that day:
    # If that is the case, current Line will only have 1 point, so we need to calculate
    # density around one point only.
      
    if (currentTrajLength > 1) {
      # if daily trajectory has more than one point, calculate density for each segment
      density <- DensityAroundTrajectory(currentLine,kernelSize,kernelBoundary,xcoord,ycoord,zcoord,method,voxelSize,sigma11,sigma12)   

    } else {
      # if daily trajectory has only one point, calculate density around the point      
      density <- DensityAroundOnePoint(currentLine,kernelSize,kernelBoundary,xcoord,ycoord,zcoord,method,voxelSize,sigma11,sigma12)
    
      } # if (currentTrajLength > 1)
    
    # Add density for this row to total shark density 
    sharkDensity <- sharkDensity+density
    
  } # for i in trajNo

  # Normalise shark density with no of trajectories (days)
  sharkDensity  <- sharkDensity/trajNo

  # Multiply shark density with activity probability
  sharkActivityDensity <- sharkDensity * activityVolume
  
  
  # Export density of this individual shark
  sharkDensity.file <- paste(c(listOfSharks[shk],"_Density.csv"),collapse="")
  WriteCSVtableFromFiveRectangular3DArrays(xcoord,ycoord,zcoord,sharkDensity,sharkActivityDensity,sharkDensity.file)
  
  # Add this one individual's density to total Density

  totalDensity <- totalDensity + sharkDensity
  
} # for (shk in 1:noOfSharks)

# ===== 21 March 2016: 
# At the end we need to normalise the total density also w/no. of sharks

totalDensity  <- totalDensity/noOfSharks

# Multiply shark density with activity probability
totalActivityDensity <- totalDensity * activityVolume

# ------------------------------------------------
# Quickly test plot the results to check what happened

# Generate colour scheme
# palPurple <- brewer.pal(9,'Purples')
# colPal <- colorRampPalette(palPurple[2:9])(length(unique(totalDensity)))

# scatter3D(xcoord,ycoord,zcoord,colvar=totalDensity,pch=".",cex=2,col=colPal)

# ------------------------------------------------
# Export density into csv file that can be imported into Voxler.

print('Exporting total density into a csv file.')
WriteCSVtableFromFiveRectangular3DArrays(xcoord,ycoord,zcoord,totalDensity,totalActivityDensity,output.file)

# ------------------------------------------------


