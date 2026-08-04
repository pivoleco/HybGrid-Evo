# ==============================================================================
# CREATE INITIAL HYBGRIDR SIMULATION ENVIRONMENT
# ==============================================================================
#
# Purpose
# -------
# Create two parental populations for the spatial hybridization simulation:
#
#   Population A:
#     Adapted through loci 101–200 and fixed for allele 1 at loci 201–280.
#
#   Population B:
#     Adapted through loci 1–100 and fixed for allele 0 at loci 201–280.
#
# The populations undergo separate burn-in periods under selection. The complete
# initialization procedure is repeated until their mean genotypic values differ
# by no more than the specified tolerance.
#
# The resulting objects are saved to an .RData environment file that can be
# loaded by the main simulation script.
#
# Locus structure
# ---------------
#   1–100:   fitness loci associated with population B
#   101–200: fitness loci associated with population A
#   201–230: epistatic loci
#   201–280: population-informative/classification loci
#   281–300: additional neutral loci
#
# Required package
# ----------------
#   AlphaSimR
#
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Package check
# ------------------------------------------------------------------------------

if (!requireNamespace("AlphaSimR", quietly = TRUE)) {
  stop(
    "Package 'AlphaSimR' is required.",
    "\nInstall it before running this script.",
    call. = FALSE
  )
}


# ------------------------------------------------------------------------------
# 2. Simulation parameters
# ------------------------------------------------------------------------------

base_seed <- 1223L
set.seed(base_seed)

# Main spatial simulation parameters
grid_size       <- 20L
n_gen           <- 50L
n_iter          <- 10L
migration_rate  <- 0.30
n_ind           <- 100L
preference_prob <- 0.10
nProgeny        <- 10L
threshold       <- 20

# Genome configuration
n_chromosomes <- 1L
seg_sites     <- 300L
map_interval  <- 0.5

# Burn-in configuration
burnin_generations <- 100L
burnin_n_crosses   <- 50L
burnin_n_progeny   <- 5L

# Distribution of additive effects
meanA <- 0
sdA   <- 0.06

# Distribution of dominance degrees
meanD <- 0.50
varD  <- 0.05

# Distribution used to generate epistatic effects during initialization
meanEpi <- 0.01
sdEpi   <- 0.006

# Target mean genotypic value before burn-in
target_genotypic_mean <- 0

# Required similarity between mean post-burn-in genotypic values
fitness_tolerance <- 0.01

# Prevent an unlimited repeat loop if suitable populations are not obtained
max_initialization_attempts <- 1000L

# Output environment
environment_file <- file.path("data", "envi.RData")


# ------------------------------------------------------------------------------
# 3. Validate basic parameters
# ------------------------------------------------------------------------------

if (seg_sites < 300L) {
  stop("`seg_sites` must be at least 300 because the model uses loci 1–300.",
    call. = FALSE)
}

if (n_ind < 2L) { stop("`n_ind` must be at least 2.", call. = FALSE)
}

if (migration_rate < 0 || migration_rate > 1) {
  stop("`migration_rate` must be between 0 and 1.", call. = FALSE)
}

if (preference_prob < 0 || preference_prob > 1) {
  stop("`preference_prob` must be between 0 and 1.", call. = FALSE)
}


# ------------------------------------------------------------------------------
# 4. Helper: create the founder population and genetic map
# ------------------------------------------------------------------------------

create_founder_population <- function(n_ind,n_chromosomes,seg_sites,map_interval) {
  
  # Generate founder haplotypes using MaCS.
  founder_raw <- AlphaSimR::runMacs(nInd = n_ind, nChr = n_chromosomes,
    segSites = seg_sites)
  
  # Create an initial SimParam object so that founder haplotypes can be read.
  SP_initial <- AlphaSimR::SimParam$new(founder_raw)
  founder_population <- AlphaSimR::newPop(founder_raw,simParam = SP_initial)
  founder_haplotypes <- AlphaSimR::pullSegSiteHaplo(founder_population,
    simParam = SP_initial)
  
  # Define an evenly spaced genetic map.
  #
  # With 300 loci and a spacing of 0.5 cM, positions range from
  # 0 to 149.5 cM.
  genetic_positions <- seq( from = 0, by = map_interval,length.out = seg_sites)
  genetic_map <- list(genetic_positions)
  haplotype_map <- list(founder_haplotypes)
  mapped_founders <- AlphaSimR::newMapPop(genMap = genetic_map,
    haplotypes = haplotype_map, inbred = FALSE, ploidy = 2L)
  
  # Create the final simulation-parameter object based on the custom map.
  SP <- AlphaSimR::SimParam$new(mapped_founders)
  breed <- AlphaSimR::newPop(mapped_founders, simParam = SP)
  
  list( SP = SP,breed = breed,founder_population = mapped_founders)
}


# ------------------------------------------------------------------------------
# 5. Helper: generate the genetic architecture
# ------------------------------------------------------------------------------

create_genetic_architecture <- function(
    mean_additive,
    sd_additive,
    mean_dominance,
    var_dominance,
    mean_epistasis,
    sd_epistasis
) {
  
  # Construct 15 epistatic pairs.
  #
  # Each of loci 1–15 in the epistatic block is paired once with one of
  # loci 16–30. These indices are relative to the epistatic block
  # corresponding to genomic loci 201–230.
  pairsE <- cbind(
    sample(1:15, size = 15, replace = FALSE),
    sample(16:30, size = 15, replace = FALSE)
  )
  
  # Additive effects for the 200 primary fitness loci.
  add_effA <- rnorm(n = 100,mean = mean_additive, sd = sd_additive)
  add_effB <- rnorm(n = 100,mean = mean_additive, sd = sd_additive)
  
  # Dominance effects are proportional to the absolute additive effects.
  #
  # The random multiplier represents the dominance degree.
  dominance_degreeA <- rnorm(n = 100,mean = mean_dominance,sd = sqrt(var_dominance))
  dominance_degreeB <- rnorm(n = 100,mean = mean_dominance,sd = sqrt(var_dominance))
  
  dom_effA <- abs(add_effA) * dominance_degreeA
  dom_effB <- abs(add_effB) * dominance_degreeB
  
  add_eff <- dom_eff <- numeric(200)
  add_eff[1:100]   <- add_effB
  add_eff[101:200] <- add_effA
  
  dom_eff[1:100]   <- dom_effB
  dom_eff[101:200] <- dom_effA
  
  # Epistatic effects used when balancing the initial parental populations.
  epi_eff <- rnorm( n = 15, mean = mean_epistasis, sd = sd_epistasis)
  
  list(pairsE = pairsE, add_eff = add_eff, dom_eff = dom_eff, epi_eff = epi_eff)
}


# ------------------------------------------------------------------------------
# 6. Helper: calculate complete genetic components
# ------------------------------------------------------------------------------

calculate_genetic_components <- function(population, add_eff, dom_eff,epi_eff,
    pairsE, SP) {
  
  genotypes <- AlphaSimR::pullSegSiteGeno(population, simParam = SP)
  
  # Main additive and dominance loci.
  geno_main <- genotypes[, 1:200, drop = FALSE]
  
  # Epistatic genotypes are centered from 0/1/2 to -1/0/1.
  geno_epi <- genotypes[, 201:230, drop = FALSE] - 1L
  
  additive_value <- drop((geno_main - 1L) %*% add_eff)
  dominance_value <- drop((geno_main == 1L) %*% dom_eff)
  
  pair_products <-
    geno_epi[, pairsE[, 1], drop = FALSE] *
    geno_epi[, pairsE[, 2], drop = FALSE]
  
  epistatic_value <- drop(pair_products %*% epi_eff)
  
  list(
    additive = additive_value,
    dominance = dominance_value,
    epistatic = epistatic_value,
    total = additive_value + dominance_value + epistatic_value
  )
}


# ------------------------------------------------------------------------------
# 7. Helper: burn in one parental population
# ------------------------------------------------------------------------------

burn_in_population <- function(
    population,
    selected_loci,
    additive_effects,
    dominance_effects,
    population_intercept,
    n_generations,
    n_crosses,
    n_progeny,
    target_population_size,
    SP
) {
  
  # Record genome-wide heterozygosity after each burn-in generation.
  heterozygosity <- numeric(n_generations)
  
  for (generation in seq_len(n_generations)) {
    
    # Generate a larger offspring population before selection.
    population <- AlphaSimR::randCross( population,nCrosses = n_crosses,
      nProgeny = n_progeny,simParam = SP)
    
    genotypes_all <- AlphaSimR::pullSegSiteGeno(population, simParam = SP)
    
    parental_genotypes <- genotypes_all[, selected_loci,drop = FALSE]
    additive_value <- drop((parental_genotypes - 1L) %*% additive_effects)
    dominance_value <- drop((parental_genotypes == 1L) %*% dominance_effects)
    fitness <-  population_intercept +  additive_value +  dominance_value
    n_selected <- min(target_population_size, population@nInd)
    selected_indices <- order(fitness,decreasing = TRUE )[seq_len(n_selected)]
    population <- population[selected_indices]
    
    # This is genome-wide heterozygosity across all 300 loci.
    selected_genotypes <- AlphaSimR::pullSegSiteGeno( population, simParam = SP)
    
    heterozygosity[generation] <- mean( selected_genotypes == 1L, na.rm = TRUE)
  }
  
  list(
    population = population,
    heterozygosity = heterozygosity
  )
}


# ------------------------------------------------------------------------------
# 8. Helper: compare final mean genotypic values
# ------------------------------------------------------------------------------

symmetric_relative_difference <- function(value_A, value_B) {
  denominator <- abs(value_A) + abs(value_B)
  if (denominator <= .Machine$double.eps) { return(0)
  }
  2 * abs(value_A - value_B) / denominator
}


# ------------------------------------------------------------------------------
# 9. Create and balance parental populations
# ------------------------------------------------------------------------------

initialization_successful <- FALSE

for (initialization_attempt in seq_len(max_initialization_attempts)) {
  
  message(
    "Initialization attempt ",
    initialization_attempt,
    "/",
    max_initialization_attempts
  )
  
  # ---------------------------------------------------------------------------
  # 9.1 Create a new founder genome and simulation parameters
  # ---------------------------------------------------------------------------
  
  founder_setup <- create_founder_population(
    n_ind = n_ind,
    n_chromosomes = n_chromosomes,
    seg_sites = seg_sites,
    map_interval = map_interval
  )
  
  SP <- founder_setup$SP
  breed <- founder_setup$breed
  
  # ---------------------------------------------------------------------------
  # 9.2 Generate additive, dominance, and epistatic effects
  # ---------------------------------------------------------------------------
  
  architecture <- create_genetic_architecture(
    mean_additive = meanA,
    sd_additive = sdA,
    mean_dominance = meanD,
    var_dominance = varD,
    mean_epistasis = meanEpi,
    sd_epistasis = sdEpi
  )
  
  pairsE  <- architecture$pairsE
  add_eff <- architecture$add_eff
  dom_eff <- architecture$dom_eff
  epi_eff <- architecture$epi_eff
  
  # ---------------------------------------------------------------------------
  # 9.3 Create parental populations
  # ---------------------------------------------------------------------------
  
  breedA <- breed
  breedB <- breed
  
  classification_loci <- 201:280
  
  # Population A is fixed for allele 1 at the population-informative loci.
  breedA <- AlphaSimR::editGenome(
    pop = breedA,
    ind = seq_len(breedA@nInd),
    chr = rep(1L, length(classification_loci)),
    segSites = classification_loci,
    allele = 1,
    simParam = SP
  )
  
  # Population B is fixed for allele 0 at the population-informative loci.
  breedB <- AlphaSimR::editGenome(
    pop = breedB,
    ind = seq_len(breedB@nInd),
    chr = rep(1L, length(classification_loci)),
    segSites = classification_loci,
    allele = 0,
    simParam = SP
  )
  
  # ---------------------------------------------------------------------------
  # 9.4 Calculate an intercept that centers the initial genotypic value
  # ---------------------------------------------------------------------------
  
  initial_components <- calculate_genetic_components(
    population = breedA,
    add_eff = add_eff,
    dom_eff = dom_eff,
    epi_eff = epi_eff,
    pairsE = pairsE,
    SP = SP
  )
  
  pop_mean <- target_genotypic_mean - mean(initial_components$total)
  
  # ---------------------------------------------------------------------------
  # 9.5 Burn in population B
  # ---------------------------------------------------------------------------
  
  burnin_B <- burn_in_population(
    population = breedB,
    selected_loci = 1:100,
    additive_effects = add_eff[1:100],
    dominance_effects = dom_eff[1:100],
    population_intercept = pop_mean,
    n_generations = burnin_generations,
    n_crosses = burnin_n_crosses,
    n_progeny = burnin_n_progeny,
    target_population_size = n_ind,
    SP = SP
  )
  
  breedB <- burnin_B$population
  hbB <- burnin_B$heterozygosity
  
  # ---------------------------------------------------------------------------
  # 9.6 Burn in population A
  # ---------------------------------------------------------------------------
  
  burnin_A <- burn_in_population(
    population = breedA,
    selected_loci = 101:200,
    additive_effects = add_eff[101:200],
    dominance_effects = dom_eff[101:200],
    population_intercept = pop_mean,
    n_generations = burnin_generations,
    n_crosses = burnin_n_crosses,
    n_progeny = burnin_n_progeny,
    target_population_size = n_ind,
    SP = SP
  )
  
  breedA <- burnin_A$population
  hbA <- burnin_A$heterozygosity
  
  # ---------------------------------------------------------------------------
  # 9.7 Recalculate complete post-burn-in genotypic values
  # ---------------------------------------------------------------------------
  
  components_A <- calculate_genetic_components(
    population = breedA,
    add_eff = add_eff,
    dom_eff = dom_eff,
    epi_eff = epi_eff,
    pairsE = pairsE,
    SP = SP
  )
  
  components_B <- calculate_genetic_components(
    population = breedB,
    add_eff = add_eff,
    dom_eff = dom_eff,
    epi_eff = epi_eff,
    pairsE = pairsE,
    SP = SP
  )
  
  gen_valA <- pop_mean + components_A$total
  gen_valB <- pop_mean + components_B$total
  
  fit_totA <- gen_valA
  fit_totB <- gen_valB
  
  mean_fitness_A <- mean(fit_totA)
  mean_fitness_B <- mean(fit_totB)
  
  relative_fitness_difference <- symmetric_relative_difference(
    mean_fitness_A,
    mean_fitness_B
  )
  
  message(
    sprintf(
      paste0(
        "Mean A = %.5f | Mean B = %.5f | ",
        "relative difference = %.4f"
      ),
      mean_fitness_A,
      mean_fitness_B,
      relative_fitness_difference
    )
  )
  
  if (relative_fitness_difference <= fitness_tolerance) {
    initialization_successful <- TRUE
    break
  }
}


# ------------------------------------------------------------------------------
# 10. Confirm successful initialization
# ------------------------------------------------------------------------------

if (!initialization_successful) {
  stop(
    "No pair of sufficiently similar parental populations was obtained after ",
    max_initialization_attempts,
    " initialization attempts.",
    call. = FALSE
  )
}

message(
  "Suitable parental populations found after ",
  initialization_attempt,
  " attempt(s)."
)

# ------------------------------------------------------------------------------
# 11. Split effects into population-specific vectors
# ------------------------------------------------------------------------------

# Loci 1–100 contribute to adaptation of population B.
add_effB <- add_eff[1:100]
dom_effB <- dom_eff[1:100]

# Loci 101–200 contribute to adaptation of population A.
add_effA <- add_eff[101:200]
dom_effA <- dom_eff[101:200]

# ------------------------------------------------------------------------------
# 12. Record reproducibility information
# ------------------------------------------------------------------------------

created_at <- Sys.time()
session_info <- utils::sessionInfo()

# ------------------------------------------------------------------------------
# 13. Save the simulation environment
# ------------------------------------------------------------------------------

dir.create( dirname(environment_file), recursive = TRUE,showWarnings = FALSE)

save(  # AlphaSimR objects
  SP, breedA, breedB,
  
  # Main simulation settings
  grid_size, n_gen, n_iter, migration_rate, n_ind, preference_prob, nProgeny,
  threshold,
  
  # Genetic architecture
  seg_sites,  pairsE,  add_eff, dom_eff, add_effA, dom_effA, add_effB, dom_effB,
  pop_mean,
  
  # Burn-in settings and diagnostics
  burnin_generations,  burnin_n_crosses,burnin_n_progeny, hbA, hbB,
  mean_fitness_A, mean_fitness_B, relative_fitness_difference,
  initialization_attempt,
  
  # Reproducibility
  base_seed, created_at, session_info,
  
  file = environment_file
)

message("Environment saved to: ", normalizePath(environment_file, mustWork = TRUE))