#This function accepts the theoretical three highest peaks and the 
#experimental three highest peaks. It returns TRUE if the main peaks
#are less than 0.05 apart and returns FALSE if they are more than
#0.05 apart.
compareThreeHighest <- function(matrix1,matrix2) {
  eps = 0.05
  if(abs(matrix1[1] - matrix2[1]) < eps){
    return(TRUE)
  }
  else{
    return(FALSE)
  }
}