# ==============================================================================
# Assortative-mating function
# ==============================================================================
#
# Forms mating pairs according to genetic similarity at population-informative
# loci 201–280.
#
# Arguments
# ---------
# deme:
#   AlphaSimR population object containing the individuals available for mating.
#
# threshold:
#   Maximum permitted genetic distance between two parents. A pair is accepted
#   when its summed genotype mismatch across loci 201–280 is strictly smaller
#   than this threshold.
#
# nProgeny:
#   Number of offspring produced by each accepted parental pair.
#
# SP:
#   AlphaSimR simulation-parameter object passed explicitly to makeCross().
#
# Procedure
# ---------
# The mating procedure is repeated for up to three rounds. In each round:
#
#   1. Individuals are randomly paired.
#   2. Genetic distance is calculated between the members of each pair.
#   3. Pairs below the compatibility threshold reproduce.
#   4. Parents that reproduced are removed before the next round.
#
# The function returns a merged AlphaSimR population containing all offspring.
# If no compatible pairs are found during any round, it returns NULL.
# ==============================================================================

mating3 <- function(deme, threshold, nProgeny, SP) {
  
  # Store offspring populations produced during individual mating rounds.
  offspring_list <- list()
  
  # Allow unmatched individuals to be paired again for up to three rounds.
  for (round in 1:3) {
    
    # At least four individuals are required to continue forming pairs.
    if (deme@nInd < 4) break
    
    # Extract genotypes at the population-informative loci used to determine
    # mating compatibility.
    genoM <- pullSegSiteGeno(deme)[, 201:280]
    
    # Randomize the order of individuals before forming mating pairs.
    indices <- sample(1:deme@nInd)
    
    # Determine the maximum number of complete pairs.
    nPairs <- floor(length(indices) / 2)
    
    if (nPairs == 0) break
    
    # Arrange the randomly ordered individual indices into two-parent pairs.
    pairs <- matrix(indices[1:(2 * nPairs)], ncol = 2, byrow = TRUE)
    
    # Calculate genetic distance within each pair as the summed absolute
    # difference in genotype scores across loci 201–280.
    distances <- rowSums(abs(genoM[pairs[, 1], ] - genoM[pairs[, 2], ]))
    
    # Retain pairs whose genetic distance is below the mating threshold.
    compatible_idx <- which(distances < threshold)
    
    # Continue to the next round when no compatible pairs are found.
    if (length(compatible_idx) == 0) next
    
    compatible_pairs <- pairs[compatible_idx, ,drop = FALSE]
    
    # Produce offspring from all compatible parental pairs.
    off <- makeCross( deme,crossPlan = compatible_pairs, simParam = SP,
      nProgeny = nProgeny)
    
    # Store offspring from the current mating round.
    offspring_list[[round]] <- off
    
    # Identify parents that successfully reproduced.
    usedP <- unique(c(off@mother, off@father))
    
    # Remove used parents so that each parent reproduces at most once.
    # Remaining individuals may be paired again in the next round.
    deme <- deme[!(deme@id %in% usedP)]
  }
  
  # Merge offspring produced across all successful mating rounds.
  if (length(offspring_list) > 0) {
    deme <- mergePops(offspring_list)
  } else {
    # Return NULL when no compatible pair reproduced.
    return(NULL)
  }
}