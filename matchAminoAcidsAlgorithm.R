#This function matches two "threeHighestVal" vectors and outputs the sequence that matches and either "Sequence Matches"
# or "Sequence [+- 1 or 2 Amino acids] Matches" or "No Match Found" and a sequence with all zeroes.
matchAminoAcidsAlgorithm <- function(v, dataMatrix) {
  
  aminoAcids = c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
  
  #This is the case where the original sequence and test data match
  if(compareThreeHighest(aminoAcidToThreeHighest(v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15],v[16],v[17],v[18],v[19],v[20]), dataMatrix)) {
      return(c(aminoAcidToThreeHighest(v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15],v[16],v[17],v[18],v[19],v[20]), "Sequence Matches"))
  }

  #This tests the case where the experimental data is +1 amino acid
  for (i in 1:20) {
    v[i] = v[i] + 1
      if(compareThreeHighest(aminoAcidToThreeHighest(v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15],v[16],v[17],v[18],v[19],v[20]), dataMatrix)) {
        letter = aminoAcids[i]
        str = paste("+1", letter, "Matches.")
        return(aminoAcidToThreeHighest(v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15],v[16],v[17],v[18],v[19],v[20]), str)
      }
    v[i] = v[i] - 1
  }
  
  #This tests the case where the experimental data is -1 amino acid (but this can only happen if there is one or more amino acids)
  for (i in 1:20) {
    if(v[i] >= 1) {
      
      v[i] = v[i] - 1
      if(compareThreeHighest(aminoAcidToThreeHighest(v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15],v[16],v[17],v[18],v[19],v[20]), dataMatrix)) {
        
        letter = aminoAcids[i]
        str = paste("-1", letter, "Matches.")
        return(aminoAcidToThreeHighest(v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15],v[16],v[17],v[18],v[19],v[20]), str)
      }
      v[i] = v[i] + 1
    }
  }
  
  #This tests the case where the experimental data is +2 amino acids
  for (i in 1:20) {
    v[i] = v[i] + 2
    if(compareThreeHighest(aminoAcidToThreeHighest(v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15],v[16],v[17],v[18],v[19],v[20]), dataMatrix)) {
      letter = aminoAcids[i]
      str = paste("+2", letter, "Matches.")
      return(aminoAcidToThreeHighest(v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15],v[16],v[17],v[18],v[19],v[20]), str)
    }
    v[i] = v[i] - 2
  }
  
  #This tests the case where the experimental data is -2 amino acids (but this can only happen if there is two or more amino acids)
  for (i in 1:20) {
    if(v[i] >= 2) {
      v[i] = v[i] - 2
      if(compareThreeHighest(aminoAcidToThreeHighest(v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15],v[16],v[17],v[18],v[19],v[20]), dataMatrix)) {
        
        letter = aminoAcids[i]
        str = paste("-2", letter, "Matches.")
        return(aminoAcidToThreeHighest(v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15],v[16],v[17],v[18],v[19],v[20]), str)
      }
      v[i] = v[i] - 2
    }
  }
  #This is the default return if no sequence matches
  return(c(0,0,0,0,0,0,"No Match Found"))
}