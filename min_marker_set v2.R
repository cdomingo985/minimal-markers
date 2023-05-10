"""
MINIMUM MARKER SETS 

min_marker_set.R

Program purpose: 
    
    From a SNP large dataset genotypied collection of varieties a list of minimum sets 
    of markers that could unequivocally identify each of the varities 

Authors:    
    Reig Valiente, Juan Luis
    Domingo Velasco, Concha
    Garcia-Romeral, Julia

    
Version 1.0

Date: 2/10/2022

"""


###################################################################################################
##############################               PRE-STEPS               ##############################
###################################################################################################

# Load required packages --------------------------------------------------------------------------
library(plyr)

# Load data ---------------------------------------------------------------------------------------
VarietieSet <- read.csv("matriz_1317_trasp_filtr.txt", stringsAsFactors = F, sep = "\t")

# Parameters to tune the analysis ------------------------------------------------------------------
# Limit the set number to analyse. If TRUE after making the combinations will discard the excess of sets
limitofSets2check <- TRUE 
# Number of sets that will limit after making the combinations
NumberofSetsLimit <- 100000 
# If TRUE the sets will be random choose, if FALSE will choose the first sets from the list generated
SelectRamdom <- TRUE 
# If TRUE will remove the failed genotypied positions (they have to be annotated as NA)
IncludeDeletionsinAnalysis <- TRUE 


###################################################################################################
##########################               MARKER FILTERING               ###########################
###################################################################################################

# Remove or not the deletions ---------------------------------------------------------------------
if(IncludeDeletionsinAnalysis==TRUE){
  VarietieSet[is.na(VarietieSet)] <- "del"
}else{
  VarietieSet <- VarietieSet[,apply(VarietieSet, 2, function(x) !(TRUE %in% is.na(x)))]    
}

# Filtering markers that cannot distinguish between varieties -------------------------------------
Markers2Keep <- apply(VarietieSet,2,function(x){length(unique(x)) > 1})
print(paste("Total markers in set are",ncol(VarietieSet), ",",sum(Markers2Keep), "will keep",
            ncol(VarietieSet) - sum(Markers2Keep),"will be removed."))

MarkerFilteredMarkerVarSet <- VarietieSet[,Markers2Keep]
print(paste("Total markers in set are", ncol(MarkerFilteredMarkerVarSet)))


###################################################################################################
##########################               VARIETY FILTERING               ##########################
###################################################################################################


# Filtering varieties that cannot be distinguish -----------------------------------------------------
# If there are varieties that cannot be distinguished we keep only the first
VarMarkerFilteredSet <- MarkerFilteredMarkerVarSet
NonDiferencialbleTaxons <- list()
Var1 <-1
while(Var1 < nrow(VarMarkerFilteredSet)){
  Taxon2Remove <- c()
  for(Var2 in (Var1+1):nrow(VarMarkerFilteredSet)){
    if(!(FALSE %in% (VarMarkerFilteredSet[Var1,] == VarMarkerFilteredSet[Var2,]))){
      print(paste("I cant differenciate", row.names(VarMarkerFilteredSet)[Var1]," from ",
                  row.names(VarMarkerFilteredSet)[Var2],sep=""))
      Taxon2Remove <- append(Taxon2Remove,Var2)
    }
  }
  if(length(Taxon2Remove)>0){
    NonDiferencialbleTaxons <- append(NonDiferencialbleTaxons,
                                      list(c(row.names(VarMarkerFilteredSet)[Var1],
                                             row.names(VarMarkerFilteredSet)[Taxon2Remove])))
    VarMarkerFilteredSet <- VarMarkerFilteredSet[-Taxon2Remove,]
    Var1 <- Var1+1
  }else{
    Var1 <- Var1+1
  }
}

# Pairs of varieties distinguished by markers -----------------------------------------------------
Varcombinations <- combn(c(1:nrow(VarMarkerFilteredSet)),2,simplify = F)
DiffMarkersInPairsComparisions <- apply(VarMarkerFilteredSet,2, function(x){
  unlist(lapply(Varcombinations, function(y){
    x[y[1]] != x[y[2]]
  }))
} )
rownames(DiffMarkersInPairsComparisions) <- unlist(lapply(Varcombinations, 
                                                          function(x){
                                                            paste(rownames(VarMarkerFilteredSet)[x[1]],
                                                                  "vs",rownames(VarMarkerFilteredSet)[x[2]])}))
colnames(DiffMarkersInPairsComparisions) <- colnames(VarMarkerFilteredSet)


###################################################################################################
############               GENERATION OF THE MINIMAL SETS OF MARKERS               ################
###################################################################################################

# Redundant markers -------------------------------------------------------------------------------
# Which markers that distinguished the same pair of varieties and choose the first marker we find 
RedundantMarkers <- list()
FilteredDiffMarkersInPairsComparisions <- DiffMarkersInPairsComparisions
i <- 1
while(i < ncol(FilteredDiffMarkersInPairsComparisions)){
  Markers2Remove <- c()
  for(x in (i+1):ncol(FilteredDiffMarkersInPairsComparisions)){
    if(!(FALSE %in% (FilteredDiffMarkersInPairsComparisions[,i]==FilteredDiffMarkersInPairsComparisions[,x]))){
      Markers2Remove <- append(Markers2Remove,x)
    }
  }
  if(length(Markers2Remove) > 0 ){
    RedundantMarkers <- append(RedundantMarkers,
                               list(c(colnames(FilteredDiffMarkersInPairsComparisions)[i],
                                      colnames(FilteredDiffMarkersInPairsComparisions)[Markers2Remove])))
    FilteredDiffMarkersInPairsComparisions <- FilteredDiffMarkersInPairsComparisions[,-Markers2Remove]
    i <- i+1
  }else{
    i <- i+1
  }
}
print(paste("MarkerFilteredMarkerVarSet contains ",ncol(MarkerFilteredMarkerVarSet), 
            "markers, after removing redundant markers",ncol(FilteredDiffMarkersInPairsComparisions),
            "will remain"))
print(paste(ncol(MarkerFilteredMarkerVarSet) - ncol(FilteredDiffMarkersInPairsComparisions),
            "have been removed"))

# Pairs of varieties distinguished by markers -----------------------------------------------------
# Extract again the markers that distinguish the pairs of varieties, but in the variety set were the
# redundant markers are not present
OrderedFDMPC <- FilteredDiffMarkersInPairsComparisions[order(apply(FilteredDiffMarkersInPairsComparisions, 1, sum)),]

firstTRUE <- matrix(which(OrderedFDMPC[1,] == TRUE),ncol = 1)
SecondTRUE <- matrix(which(OrderedFDMPC[2,] == TRUE),ncol = 1)
FirstTRUEin2 <- matrix(firstTRUE[which(firstTRUE[,1] %in% SecondTRUE[,1]),],ncol = 1)
FirstNotTRUEin2 <- matrix(firstTRUE[which(!(firstTRUE[,1] %in% SecondTRUE[,1])),],ncol = 1)

SetsCombinationNotTRUE_TRUE <- as.matrix(expand.grid(FirstNotTRUEin2,SecondTRUE))
FirstTRUEin2 <- cbind(FirstTRUEin2,rep(NA,nrow(FirstTRUEin2)))
CombineFirstTRUEin2SetsCombinationsNotTRUE_TRUE <- rbind(FirstTRUEin2, SetsCombinationNotTRUE_TRUE)

Sets <- CombineFirstTRUEin2SetsCombinationsNotTRUE_TRUE

UniqueTrue = (apply(OrderedFDMPC, 1, function(x){sum(x,na.rm = TRUE)})) < 2
Init = sum(UniqueTrue, na.rm = TRUE) + 1

print(paste("Starting to create Minimum Set, number of pairs is",nrow(OrderedFDMPC)))
print("Current Pair being analysed is:")
for(i in Init:nrow(OrderedFDMPC)){
  print(i)
  MarkersTRUEinNextPair <- matrix(which(OrderedFDMPC[i,] == TRUE),ncol = 1)
  SetsTRUEinNextPair <- Sets[which(apply(Sets, 1,function(x){TRUE %in% (x %in% MarkersTRUEinNextPair)})),]
  SetsNotTRUEinNextPair <- Sets[which(apply(Sets, 1,function(x){!(TRUE %in% (x %in% MarkersTRUEinNextPair))})),]
  if(class(SetsNotTRUEinNextPair)== "integer"){
    SetsNotTRUEinNextPair <- matrix(SetsNotTRUEinNextPair,nrow = 1,ncol = length(SetsNotTRUEinNextPair))
  }
  if(nrow(SetsNotTRUEinNextPair)> 0){
    PositionNoTRUEinNexPaircomb_TRUE <- as.matrix(expand.grid(seq_len(nrow(SetsNotTRUEinNextPair)),
                                                              MarkersTRUEinNextPair))
    SetCombinations <- cbind(SetsNotTRUEinNextPair[PositionNoTRUEinNexPaircomb_TRUE[,1],],
                             PositionNoTRUEinNexPaircomb_TRUE[,2])
    SetsTRUEinNextPair <- cbind(SetsTRUEinNextPair,rep(NA,nrow(SetsTRUEinNextPair)))
    Sets <- rbind(SetsTRUEinNextPair, SetCombinations)  
  }
  if(limitofSets2check == TRUE && nrow(Sets)>NumberofSetsLimit){
    if(SelectRamdom == TRUE){
      Sets <- Sets[sample(1:NumberofSetsLimit,replace = FALSE),]
    }else{
      Sets <- Sets[c(1:NumberofSetsLimit),]
    }
  }
}

SetsLength <- apply(Sets,1,function(x)sum(!is.na(x)))
minlengthSets <- min(SetsLength)
SetsMinLength <- Sets[SetsLength == minlengthSets,]
markersORderedFDMPC <- colnames(OrderedFDMPC)

if(sum(SetsLength == minlengthSets)==1){
  SetsMinimos <- markersORderedFDMPC[SetsMinLength[!is.na(SetsMinLength)]]
}else{
  SetsMinimos <- lapply(seq_len(nrow(SetsMinLength)), function(i) markersORderedFDMPC[SetsMinLength[i,!is.na(SetsMinLength[i,])]])
}

# Pairs of varieties that can only be distinguished by one marker ------------------------------------
# Adding those unique markers to the minimum sets (if they are not there yet) that are the unique one ables to distinguished between some pairs of varieties ------------
UniqueMatrix = OrderedFDMPC[1: (Init-1),]
MarkersUniqueList = list()
for ( i in 1:(Init-1)){
  MarkersUniqueList <- append(MarkersUniqueList,list(c(row.names(UniqueMatrix)[i],names(which(UniqueMatrix[i,] == TRUE)))))
}

for ( i in 1:(Init-1)){
  Marker = MarkersUniqueList[[i]][2] 
  
  for (j in 1:length(SetsMinimos)){
    In = Marker %in% SetsMinimos[[j]]
    
    if (In == FALSE){
      SetsMinimos[[j]] = c(SetsMinimos[[j]], Marker)
    }
  }
}

# Saving the final results -------------------------------------------------------------------------------
#List of the final minimum sets of markers
capture.output(SetsMinimos, file = "ListaSetsMinimos.txt")
#List of the redundant markers
capture.output(RedundantMarkers, file = "RedundantMarkers.txt")
#List of the varieties that cannot be distinguised
capture.output(NonDiferencialbleTaxons,file="TaxonesNodiferenciables.txt")
#Varieties that can only be differenciate between one marker
capture.output(MarkersUniqueList, file = "MarcadoresUnicos.txt")
