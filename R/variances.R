library(dplyr)

# ==============================================================================
# Genetic and genic variances
# ==============================================================================

## SETTINGS --------------------------------------------------------------------

project_dir <- ".../HybGrid-Evo"

initial_isolation   <- 10L        # choose 10, 100, 300, 500
epistasis_type      <- "bdm"      # choose "bdm" or "pos"
epistasis_magnitude <- "001"      # choose "01", "005", "001"
preference_prob     <- 0.10       # choose 0.10, 0.50, 0.90
environment_type    <- "gradient" # choose "gradient" or "homogeneous"

n_iterations <- 10L

environment_file <- file.path(project_dir, "data", "environments",
                              sprintf("envi%d.RData", initial_isolation))
scenario_name <- sprintf("isolation_%d_%s_%s_%s_%s", initial_isolation, epistasis_type,
                         epistasis_magnitude,sprintf("%02d", round(preference_prob * 10)), environment_type)
source_folder <- file.path(project_dir, "output",scenario_name)
save_folder   <- source_folder

## LOAD MODEL PARAMETERS --------------------------------------------------------

model_env <- new.env()
load(environment_file, envir = model_env)

add_effA <- get("add_effA", envir = model_env)
dom_effA <- get("dom_effA", envir = model_env)
add_effB <- get("add_effB", envir = model_env)
dom_effB <- get("dom_effB", envir = model_env)

grid_size <- as.integer(get("grid_size", envir = model_env))
n_gen     <- as.integer(get("n_gen", envir = model_env))


## LOCI -------------------------------------------------------------------------

pop_levels <- c("A", "B", "H")

loci_main  <- 1:200
loci_class <- 201:280
needed_cols <- sort(unique(c(loci_main, loci_class)))
idx_main    <- match(loci_main, needed_cols)
idx_class   <- match(loci_class, needed_cols)


## HELPERS ----------------------------------------------------------------------

safe_wmean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  weighted.mean(x[ok], w[ok])
}

# Genic variance at one x position
calc_genic <- function(geno, x, environment_type) {
  if (nrow(geno) == 0L)
    return(c(genicVA = NA_real_, genicVD = NA_real_))
  if (environment_type == "homogeneous") {
    # New simulation: full B + full A contribution
    add_eff <- c(add_effB, add_effA)
    dom_eff <- c(dom_effB, dom_effA)
    
  } else if (environment_type == "gradient") {
    weight_A <- (x - 1L) / (grid_size - 1L)
    weight_B <- 1 - weight_A
    add_eff <- c( weight_B * add_effB,weight_A * add_effA)
    dom_eff <- c( weight_B * dom_effB,weight_A * dom_effA)
    
  } else {
    stop("environment_type must be 'homogeneous' or 'gradient'")
  }
  p <- colMeans(geno) / 2
  q <- 1 - p
  genicVA <- sum(2 * p * q * (add_eff + dom_eff * (q - p))^2, na.rm = TRUE)
  genicVD <- sum((2 * p * q * dom_eff)^2,na.rm = TRUE)
  c(genicVA = genicVA, genicVD = genicVD)
}

# Genic variance across x positions
calc_genic_population <- function(geno, xpos, environment_type) {
  if (nrow(geno) == 0L || length(xpos) == 0L)
    return(c(genicVA = NA_real_, genicVD = NA_real_))
  xs <- sort(unique(xpos))
  out <- lapply(xs, function(xx) {
    ids <- xpos == xx
    st <- calc_genic( geno = geno[ids, , drop = FALSE], x = xx,
      environment_type = environment_type)
    data.frame( x = xx, n = sum(ids), genicVA = st["genicVA"],genicVD = st["genicVD"])
  })
  
  out <- bind_rows(out)
  c(genicVA = safe_wmean(out$genicVA, out$n),
    genicVD = safe_wmean(out$genicVD, out$n))
}

## OUTPUT ARRAYS ----------------------------------------------------------------

genVA <- genVD <- genVG <- array(NA_real_,
  dim = c(n_iterations, n_gen, 3),dimnames = list(iteration = seq_len(n_iterations),
    generation = seq_len(n_gen),population = pop_levels))

genicVA <- genicVD <- array(NA_real_,
  dim = c(n_iterations, n_gen, 3), dimnames = list(iteration = seq_len(n_iterations),
    generation = seq_len(n_gen),population = pop_levels))

# ==============================================================================
# 1. REALIZED GENETIC VARIANCES
# ==============================================================================

gv_file <- file.path(save_folder, "genotypic_values.rds")

if (!file.exists(gv_file))
  stop("genotypic_values.rds not found: ", gv_file)
dats <- readRDS(gv_file)

realized_var <- dats %>%
  # Selection varies only along x, so calculate within x positions first
  group_by(iteration, generation, class, deme_x) %>%
  summarise(n = n(),
    genVA = if (n() > 1L)
      var(additive_value, na.rm = TRUE)
    else NA_real_,
    genVD = if (n() > 1L)
      var(dominance_value, na.rm = TRUE)
    else NA_real_,
    genVG = if (n() > 1L)
      var(genotypic_value, na.rm = TRUE)
    else NA_real_,
    .groups = "drop"
  ) %>%
  
  # Population-size weighted mean across x positions
  group_by(iteration, generation, class) %>%
  
  summarise(
    genVA = safe_wmean(genVA, n),
    genVD = safe_wmean(genVD, n),
    genVG = safe_wmean(genVG, n),
    .groups = "drop"
  )

for (r in seq_len(nrow(realized_var))) {
  
  i   <- realized_var$iteration[r]
  j   <- realized_var$generation[r]
  pop <- realized_var$class[r]
  
  if (!pop %in% pop_levels) next
  
  genVA[i, j, pop] <- realized_var$genVA[r]
  genVD[i, j, pop] <- realized_var$genVD[r]
  genVG[i, j, pop] <- realized_var$genVG[r]
}

# ==============================================================================
# 2. GENIC VARIANCES
# ==============================================================================

for (i in seq_len(n_iterations)) {
  cat("Iteration", i, "/", n_iterations, "\n")
  for (j in seq_len(n_gen)) {
    
    hap_file <- file.path(
      source_folder,sprintf("iter_%03d", i),sprintf("gen_%04d.rds", j))
    
    if (!file.exists(hap_file)) next
    
    haps <- readRDS(hap_file)
    
    gens_list <- list()
    x_list <- list()
    k <- 0L
    
    for (x in seq_len(grid_size)) {
      for (y in seq_len(grid_size)) {
        hap <- haps[[x, y]]
        if (is.null(hap) || nrow(hap) < 2L) next
        if (ncol(hap) < max(needed_cols))
          stop(sprintf(
            "Only %d loci found in iter %d, gen %d, deme [%d,%d]; %d required.",
            ncol(hap), i, j, x, y, max(needed_cols)
          ))
        
        odd <- seq.int(1L, nrow(hap) - 1L, by = 2L)
        even <- odd + 1L
        gens <- hap[odd, needed_cols, drop = FALSE] +
          hap[even, needed_cols, drop = FALSE]
        
        k <- k + 1L
        gens_list[[k]] <- gens
        x_list[[k]] <- rep(x, nrow(gens))
      }
    }
    
    if (k == 0L) next
    
    gens_all <- do.call(rbind, gens_list)
    xpos <- unlist(x_list, use.names = FALSE)
    
    ## Population classification
    cls <- rowMeans(gens_all[, idx_class, drop = FALSE]) / 2
    
    ## Main QTL genotypes
    geno <- gens_all[, idx_main, drop = FALSE]
    ids <- list(
      A = cls >= 0.9,
      B = cls <= 0.1,
      H = cls > 0.1 & cls < 0.9)
    
    for (pop in pop_levels) {
      
      id <- ids[[pop]]
      if (!any(id)) next
      st <- calc_genic_population(geno = geno[id, , drop = FALSE],xpos = xpos[id],
        environment_type = environment_type)
      
      genicVA[i, j, pop] <- st["genicVA"]
      genicVD[i, j, pop] <- st["genicVD"]
    }
  }
}


# ==============================================================================
# 3. SAVE
# ==============================================================================

V <- list(genV = list(VA = genVA,VD = genVD,VG = genVG),genicV = list(
    VA = genicVA, VD = genicVD))
output_file <- file.path(save_folder, "genetic_variances.rds")
saveRDS(V, output_file, compress = TRUE)