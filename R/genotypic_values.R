# ==============================================================================
# Extract individual genotypic values from HybGridR haplotype output
# ==============================================================================

## SETTINGS --------------------------------------------------------------------

project_dir <- ".../HybGrid-Evo"   # change 

initial_isolation <- 10L        #(choose 10,100,300,500)
epistasis_type <- "bdm"         #(choose bdm or pos)
epistasis_magnitude <- "001"    #(choose 01,005,001)
preference_prob     <- 0.10       # choose 0.10, 0.50, 0.90
environment_type <- "gradient"  #(choose gradient or homogeneous)
n_iterations        <- 10L

input_dir <- file.path(project_dir, "output",
    sprintf("isolation_%d_%s_%s_%s_%s", initial_isolation, epistasis_type,
   epistasis_magnitude,sprintf("%02d", round(preference_prob * 10)), environment_type))

environment_file <- file.path(project_dir, "data", "environments",
                              sprintf("envi%d.RData", initial_isolation))


## LOAD MODEL PARAMETERS --------------------------------------------------------

model_env <- new.env()
load(environment_file, envir = model_env)

pairsE   <- get("pairsE", envir = model_env)
add_effA <- get("add_effA", envir = model_env)
dom_effA <- get("dom_effA", envir = model_env)
add_effB <- get("add_effB", envir = model_env)
dom_effB <- get("dom_effB", envir = model_env)

grid_size <- as.integer(get("grid_size", envir = model_env))
n_gen     <- as.integer(get("n_gen", envir = model_env))
n_ind     <- as.integer(get("n_ind", envir = model_env))
max_ind   <- as.integer(floor(n_ind / 2))


## EPISTATIC EFFECTS ------------------------------------------------------------

EPI_EFFECTS <- list(
  "01" = c(
    -0.020733651,-0.018625038,-0.006291431,-0.004034495,-0.008043853,
    -0.005442798,-0.003021580,-0.015429728,-0.003002657,-0.011406430,
    -0.010789055,-0.015517526,-0.012286533,-0.008870254,-0.017678590),
  "005" = c(
    -0.0008622194,-0.0009184874,-0.0007796196,-0.0008385601,-0.0013560065,
    -0.0002860874,-0.0012618780,-0.0003284061,-0.0007399826,-0.0009307868,
    -0.0005277535,-0.0012485223,-0.0001434993,-0.0008188440,-0.0018042037),
  "001" = c(
    -0.004008850,-0.003938429,-0.004865292,-0.005568119,-0.005254795,
    -0.004558156,-0.004817784,-0.004825235,-0.005426624,-0.005845740,
    -0.004434927,-0.005463951,-0.005877341,-0.005674150,-0.005675997)
)

base_eff <- EPI_EFFECTS[[epistasis_magnitude]]
EF <- if (epistasis_type == "bdm") abs(base_eff) else base_eff


## LOCI -------------------------------------------------------------------------

loci_B     <- 1:100
loci_A     <- 101:200
loci_epi   <- 201:230
loci_class <- 201:280

needed_cols <- sort(unique(c(loci_B, loci_A, loci_epi, loci_class)))

idx_B     <- match(loci_B, needed_cols)
idx_A     <- match(loci_A, needed_cols)
idx_epi   <- match(loci_epi, needed_cols)
idx_class <- match(loci_class, needed_cols)

p1 <- pairsE[, 1L]
p2 <- pairsE[, 2L]


## OUTPUT ARRAYS -----------------------------------------------------------------

pop_levels <- c("A", "B", "H") # A- introduced, B- native, H- Hybrids

dims <- c(n_iterations, n_gen, grid_size, grid_size, max_ind, 3L)
dn <- list(NULL, NULL, NULL, NULL, NULL, pop_levels)

gen_val <- array(NA_real_, dim = dims, dimnames = dn)
addval  <- array(NA_real_, dim = dims, dimnames = dn)
domval  <- array(NA_real_, dim = dims, dimnames = dn)
epival  <- array(NA_real_, dim = dims, dimnames = dn)


## CALCULATE VALUES --------------------------------------------------------------

for (i in seq_len(n_iterations)) {
  for (j in seq_len(n_gen)) {
    message(sprintf("Iteration %d/%d | Generation %d/%d",i, n_iterations, j, n_gen))
    haps <- readRDS(file.path(input_dir, sprintf("iter_%03d", i),
                              sprintf("gen_%04d.rds", j)))
    
    for (x in seq_len(grid_size)) {
      for (y in seq_len(grid_size)) {
        
        hap <- haps[[x, y]]
        if (is.null(hap) || nrow(hap) < 2L) next
        if (ncol(hap) < max(needed_cols))
          stop(sprintf("Only %d loci found in iter %d, gen %d, deme [%d,%d]; classification requires %d.",
                       ncol(hap), i, j, x, y, max(needed_cols)))
        
        odd <- seq.int(1L, nrow(hap) - 1L, by = 2L)
        even <- odd + 1L
        gens <- hap[odd, needed_cols, drop = FALSE] +  hap[even, needed_cols, drop = FALSE]
        
        ## Population classification
        cls <- rowSums(gens[, idx_class, drop = FALSE]) /
          (2 * length(idx_class))
        
        ids_list <- list(
          A = cls >= 0.9,
          B = cls <= 0.1,
          H = cls > 0.1 & cls < 0.9)
        
        ## Genetic components
        genoA <- gens[, idx_A, drop = FALSE]
        genoB <- gens[, idx_B, drop = FALSE]
        genoE <- gens[, idx_epi, drop = FALSE] - 1L
        
        add_A <- drop((genoA - 1L) %*% add_effA)
        dom_A <- drop((genoA == 1L) %*% dom_effA)
        
        add_B <- drop((genoB - 1L) %*% add_effB)
        dom_B <- drop((genoB == 1L) %*% dom_effB)
        
        epi_val <- drop(
          (genoE[, p1, drop = FALSE] * genoE[, p2, drop = FALSE]) %*% EF
        )
        
        ## Environmental weighting -- identical to simulation
        if (environment_type == "homogeneous") {
          
          weighted_add <- add_B + add_A
          weighted_dom <- dom_B + dom_A
          
        } else if (environment_type == "gradient") {
          
          weight_A <- (x - 1L) / (grid_size - 1L)
          weight_B <- 1 - weight_A
          weighted_add <- weight_B * add_B + weight_A * add_A
          weighted_dom <- weight_B * dom_B + weight_A * dom_A
          
        } else {
          stop("environment_type must be 'homogeneous' or 'gradient'")
        }
        
        fit_pop <- weighted_add + weighted_dom + epi_val
        
        ## Store A, B and hybrids
        for (p in seq_along(ids_list)) {
          idx <- ids_list[[p]]
          n <- sum(idx)
          
          if (n == 0L) next
          
          pos <- seq_len(n)
          
          gen_val[i,j,x,y,pos,p] <- fit_pop[idx]
          addval[i,j,x,y,pos,p]  <- weighted_add[idx]
          domval[i,j,x,y,pos,p]  <- weighted_dom[idx]
          epival[i,j,x,y,pos,p]  <- epi_val[idx]
        }
      }
    }
  }
}


## CONVERT TO DATA FRAME ---------------------------------------------------------

ind <- which(!is.na(gen_val), arr.ind = TRUE)

out <- data.frame(
  iteration       = ind[,1],
  generation      = ind[,2],
  deme_x          = ind[,3],
  deme_y          = ind[,4],
  individual      = ind[,5],
  class           = pop_levels[ind[,6]],
  genotypic_value = gen_val[ind],
  additive_value  = addval[ind],
  dominance_value = domval[ind],
  epistatic_value = epival[ind]
)

saveRDS(out, file.path(input_dir, "genotypic_values.rds"), compress = TRUE)
