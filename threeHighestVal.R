#This function takes a n-by-2 matrix and looks for the highest three values in the right-hand
#column, then outputs those values with the corresponding value in the right-hand column in
# a vector of the form (left1 right1 left2 right2 left3 right3).
threeHighestVal <- function(mSMatrix) {
  
  resultVector <- numeric(6)
  
  for (i in 0:2) {
    maxIntensityLocation <- which.max(mSMatrix[,2])
    maxMass <- mSMatrix[maxIntensityLocation,1]
    maxIntensity <- mSMatrix[maxIntensityLocation,2]
    
    resultVector[1 + 2*i] <- maxMass
    resultVector[2 + 2*i] <- maxIntensity
    
    mSMatrix[maxIntensityLocation,2] <- 0}
  
  return(resultVector)
}