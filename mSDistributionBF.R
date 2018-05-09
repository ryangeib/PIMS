#This function inputs the number of C, H, O, N, and S and outputs a theoretical mass spectrum of a compound 
#with this composition as a matrix with the masses in the left-hand column and intesnities in the right-hand column
mSDistributionBF <- function(v,w,x,y,z,intDP) {
  
  #store the molecular weight and isotopic abundances of the isotopes of each atom
  C12 <- c(0.98893,12.0000)	
  C13 <- c(0.01107,13.0034)
  
  H1 <- c(0.9985,1.0078)
  H2 <- c(0.00015,2.0141)
  
  N14 <- c(0.99634,14.0030)
  N15 <- c(0.00366,15.00001)
  
  O16 <- c(0.99759,15.9949)
  O17 <- c(0.00037,16.9991)
  O18 <- c(0.00204,17.9991)
  
  S32 <- c(0.9500,31.9720)
  S33 <- c(0.0076,32.9715)
  S34 <- c(0.0422,33.9671)
  
  dp = intDP
  
  #Record maximum and minimum mass possible of the molecule
  minMass <- round(C12[2]*v + H1[2]*w + N14[2]*x + O16[2]*y + S32[2]*z, digits = dp)
  maxMass <- round(C13[2]*v + H2[2]*w + N15[2]*x + O18[2]*y + S34[2]*z, digits = dp)
  
  #Create nBins to make a vector as large as possible with entries truncated to the hundreths place
  nBins = (maxMass - minMass)*100 + 1
  
  #Create mass string that is the length that the user desires
  mass <- numeric(nBins)
  mass[1] <- minMass
  for(i in 1:(nBins)){
    mass[i+1] <- mass[1] + i/100.0
  }
  
  #############################################################################################################################################
  
  #Fix initial probability
  prob <- numeric(nBins)
  prob[1] <- (C12[1]^v)*(H1[1]^w)*(N14[1]^x)*(O16[1]^y)*(S32[1]^z)
  
  ############## M+1 Peak calculation #############
  
    #Probability of ONE C13
    if(v > 0){
      molMass <- round(C12[2]*(v-1) + C13[2]*1 + H1[2]*w + N14[2]*x + O16[2]*y + S32[2]*z, digits = dp)
      massDif <- molMass - minMass
      massEnt <- massDif*100
      prob[massEnt] <- (C12[1]^(v-1))*(C13[1]^1)*(H1[1]^w)*(N14[1]^x)*(O16[1]^y)*(S32[1]^z) + prob[massEnt]}

    #Probability of ONE H2
    if(w > 0){
      molMass <- round(C12[2]*v + H1[2]*(w-1) + H2[2]*1 + N14[2]*x + O16[2]*y + S32[2]*z, digits = dp)
      massDif <- molMass - minMass
      massEnt <- massDif*100
      prob[massEnt] <- (C12[1]^v)*(H1[1]^(w-1))*(H2[1]^1)*(N14[1]^x)*(O16[1]^y)*(S32[1]^z) + prob[massEnt]}

    #Probability of ONE N15
    if(x > 0){
      molMass <- round(C12[2]*v + H1[2]*w + N14[2]*(x-1) + N15[2]*1 + O16[2]*y + S32[2]*z, digits = dp)
      massDif <- molMass - minMass
      massEnt <- massDif*100
      prob[massEnt] <- (C12[1]^v)*(H1[1]^w)*(N14[1]^(x-1))*(N15[1]^1)*(O16[1]^y)*(S32[1]^z) + prob[massEnt]}
  
    #Probability of ONE O17
    if(y > 0){
      molMass <- round(C12[2]*v + H1[2]*w + N14[2]*x + O16[2]*(y-1) + O17[2]*1 + S32[2]*z, digits = dp)
      massDif <- molMass - minMass
      massEnt <- massDif*100
      prob[massEnt] <- (C12[1]^v)*(H1[1]^w)*(N14[1]^x)*(O16[1]^(y-1))*(O17[1]^1)*(S32[1]^z) + prob[massEnt]}

    #Probability of ONE S33
    if(z > 0){
      molMass <- round(C12[2]*v + H1[2]*w + N14[2]*x + O16[2]*y + S32[2]*(z-1) + S33[2]^1, digits = dp)
      massDif <- molMass - minMass
      massEnt <- massDif*100
      prob[massEnt] <- (C12[1]^v)*(H1[1]^w)*(N14[1]^x)*(O16[1]^(y))*(S32[1]^(z-1))*(S33[1]^1) + prob[massEnt]}
  
  ############## M+2 Peak calculation #############

    #Probability of TWO C13
    if(v > 1){
      molMass <- round(C12[2]*(v-2) + C13[2]*2 + H1[2]*w + N14[2]*x + O16[2]*y + S32[2]*z, digits = dp)
      massDif <- molMass - minMass
      massEnt <- massDif*100
      prob[massEnt] <- (C12[1]^(v-2))*(C13[1]^2)*(H1[1]^w)*(N14[1]^x)*(O16[1]^y)*(S32[1]^z) + prob[massEnt]}

    #Probability of TWO H2
    if(w > 1){
      molMass <- round(C12[2]*v + H1[2]*(w-2) + H2[2]*2 + N14[2]*x + O16[2]*y + S32[2]*z, digits = dp)
      massDif <- molMass - minMass
      massEnt <- massDif*100
      prob[massEnt] <- (C12[1]^v)*(H1[1]^(w-2))*(H2[1]^2)*(N14[1]^x)*(O16[1]^y)*(S32[1]^z) + prob[massEnt]}

    #Probability of TWO N15
    if(x > 1){
      molMass <- round(C12[2]*v + H1[2]*w + N14[2]*(x-2) + N15[2]*2 + O16[2]*y + S32[2]*z, digits = dp)
      massDif <- molMass - minMass
      massEnt <- massDif*100
      prob[massEnt] <- (C12[1]^v)*(H1[1]^w)*(N14[1]^(x-2))*(N15[1]^2)*(O16[1]^y)*(S32[1]^z) + prob[massEnt]}

    #Probability of TWO O17 or ONE O18
    if(y > 1){
      molMass <- round(C12[2]*v + H1[2]*w + N14[2]*x + O16[2]*(y-2) + O17[2]*2 + S32[2]*z, digits = dp)
      massDif <- molMass - minMass
      massEnt <- massDif*100
      prob[massEnt] <- (C12[1]^v)*(H1[1]^w)*(N14[1]^x)*(O16[1]^(y-2))*(O17[1]^2)*(S32[1]^z) + prob[massEnt]}
    if(y > 0){
      molMass <- round(C12[2]*v + H1[2]*w + N14[2]*x + O16[2]*(y-1) + O18[2]*1 + S32[2]*z, digits = dp)
      massDif <- molMass - minMass
      massEnt <- massDif*100
      prob[massEnt] <- (C12[1]^v)*(H1[1]^w)*(N14[1]^x)*(O16[1]^(y-1))*(O18[1]^1)*(S32[1]^z) + prob[massEnt]}

    #Probability of TWO S33 or ONE S34
    if(z > 1){
      molMass <- round(C12[2]*v + H1[2]*w + N14[2]*x + O16[2]*y + S32[2]*(z-2) + S33[2]^2, digits = dp)
      massDif <- molMass - minMass
      massEnt <- massDif*100
      prob[massEnt] <- (C12[1]^v)*(H1[1]^w)*(N14[1]^x)*(O16[1]^y)*(S32[1]^(z-2))*(S33[1]^2) + prob[massEnt]}
    if(z > 0){
      molMass <- round(C12[2]*v + H1[2]*w + N14[2]*x + O16[2]*y + S32[2]*(z-1) + S34[2]^1, digits = dp)
      massDif <- molMass - minMass
      massEnt <- massDif*100
      prob[massEnt] <- (C12[1]^v)*(H1[1]^w)*(N14[1]^x)*(O16[1]^y)*(S32[1]^(z-1))*(S34[1]^1) + prob[massEnt]}

    #Probability of ONE C13 and ONE H2
    if(v > 0 & w > 0){
      molMass <- round(C12[2]*(v-1) + C13[2]*1 + H1[2]*(w-1) + H2[2]*1 + N14[2]*x + O16[2]*y + S32[2]*z, digits = dp)
      massDif <- molMass - minMass
      massEnt <- massDif*100
      prob[massEnt] <- (C12[1]^(v-1))*(C13[1]^1)*(H1[1]^(w-1))*(H2[1]^1)*(N14[1]^x)*(O16[1]^y)*(S32[1]^z) + prob[massEnt]}

    #Probability of ONE C13 and ONE N15
    if(v > 0 & x > 0){
      molMass <- round(C12[2]*(v-1) + C13[2]*1 + H1[2]*w + N14[2]*(x-1) + N15[2]*1 + O16[2]*y + S32[2]*z, digits = dp)
      massDif <- molMass - minMass
      massEnt <- massDif*100
      prob[massEnt] <- (C12[1]^(v-1))*(C13[1]^1)*(H1[1]^w)*(N14[1]^(x-1))*(N15[1]^1)*(O16[1]^y)*(S32[1]^z) + prob[massEnt]}

    #Probability of ONE C13 and ONE O17
    if(v > 0 & y > 0){
      molMass <- round(C12[2]*(v-1) + C13[2]*1 + H1[2]*w + N14[2]*x + O16[2]*(y-1) + O17[2]*1 + S32[2]*z, digits = dp)
      massDif <- molMass - minMass
      massEnt <- massDif*100
      prob[massEnt] <- (C12[1]^(v-1))*(C13[1]^1)*(H1[1]^w)*(N14[1]^x)*(O16[1]^(y-1))*(O17[1]^1)*(S32[1]^z) + prob[massEnt]}

    #Probability of ONE C13 and ONE O17
    if(v > 0 & z > 0){
      molMass <- round(C12[2]*(v-1) + C13[2]*1 + H1[2]*w + N14[2]*x + O16[2]*y + S32[2]*(z-1) + S33[2]*1, digits = dp)
      massDif <- molMass - minMass
      massEnt <- massDif*100
      prob[massEnt] <- (C12[1]^(v-1))*(C13[1]^1)*(H1[1]^w)*(N14[1]^x)*(O16[1]^y)*(S32[1]^(z-1))*(S33[1]^1) + prob[massEnt]}
  
    #Probability of ONE H2 and ONE N15
    if(w > 0 & x > 0){
      molMass <- round(C12[2]*v + H1[2]*(w-1) + H2[2]*1 + N14[2]*(x-1) + N15[2]*1 + O16[2]*y + S32[2]*z, digits = dp)
      massDif <- molMass - minMass
      massEnt <- massDif*100
      prob[massEnt] <- (C12[1]^v)*(H1[1]^(w-1))*(H2[1]^1)*(N14[1]^(x-1))*(N15[1]^1)*(O16[1]^y)*(S32[1]^z) + prob[massEnt]}

    #Probability of ONE H2 and ONE O17
    if(w > 0 & y > 0){
      molMass <- round(C12[2]*v + H1[2]*(w-1) + H2[2]*1 + N14[2]*x + O16[2]*(y-1) + O17[2]*1 + S32[2]*z, digits = dp)
      massDif <- molMass - minMass
      massEnt <- massDif*100
      prob[massEnt] <- (C12[1]^v)*(H1[1]^(w-1))*(H2[1]^1)*(N14[1]^x)*(O16[1]^(y-1))*(O17[1]^1)*(S32[1]^z) + prob[massEnt]}

    #Probability of ONE H2 and ONE S33
    if(w > 0 & y > 0){
      molMass <- round(C12[2]*v + H1[2]*(w-1) + H2[2]*1 + N14[2]*x + O16[2]*y + S32[2]*(z-1) + S33[2]*1, digits = dp)
      massDif <- molMass - minMass
      massEnt <- massDif*100
      prob[massEnt] <- (C12[1]^v)*(H1[1]^(w-1))*(H2[1]^1)*(N14[1]^x)*(O16[1]^y)*(S32[1]^(z-1))*(S33[1]*1) + prob[massEnt]}

    #Probability of ONE N15 and ONE O17
    if(x > 0 & y > 0){
      molMass <- round(C12[2]*v + H1[2]*w + N14[2]*(x-1) + N15[2]*1 + O16[2]*(y-1) + O17[2]*1 + S32[2]*z, digits = dp)
      massDif <- molMass - minMass
      massEnt <- massDif*100
      prob[massEnt] <- (C12[1]^v)*(H1[1]^w)*(N14[1]^(x-1))*(N15[1]^1)*(O16[1]^(y-1))*(O17[1]^1)*(S32[1]^z) + prob[massEnt]}

    #Probability of ONE N15 and ONE O17
    if(x > 0 & z > 0){
      molMass <- round(C12[2]*v + H1[2]*w + N14[2]*(x-1) + N15[2]*1 + O16[2]*y + S32[2]*(z-1) + S33[2]*1, digits = dp)
      massDif <- molMass - minMass
      massEnt <- massDif*100
      prob[massEnt] <- (C12[1]^v)*(H1[1]^w)*(N14[1]^(x-1))*(N15[1]^1)*(O16[1]^y)*(S32[1]^(z-1))*(S33[1]^1) + prob[massEnt]}
  
  #prob[nBins + 1] = 0
  maxEntry <- which.max(prob)
  
  #Translate probability into intensities
  intensities = numeric(nBins)
  for(i in 1:(nBins + 1)){
    intensities[i] <- round((prob[i]/prob[maxEntry])*100, digits = intDP)}
  
  #Result makes matrix of masses with relative intensities
  result <- matrix(c(mass,intensities), nrow <- length(mass))
  
  return(result)
}