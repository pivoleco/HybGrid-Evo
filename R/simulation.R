# ==============================================================================
# HybGridR: spatial hybridization simulation
# ==============================================================================
#
# Purpose
# -------
# Simulate migration, reproduction, selection, and hybridization between two
# parental populations distributed across a two-dimensional grid.
#
# Required input files
# --------------------
# The .RData file specified in `config$environment_file` must contain:
#
#   breedA, breedB, grid_size, n_gen, migration_rate, nProgeny, n_ind,
#   pairsE, add_effA, dom_effA, add_effB, dom_effB, threshold,
#   and either SP or SimParam.
#
# The R script specified in `config$mating_script` must define `reproduction()`
# and any helper functions used by it. It is sourced into the same environment
# as the objects loaded from the .RData file.
#
# Output
# ------
# One directory per iteration:
#
#   output/haplotypes/iter_001/gen_0001.rds
#   output/haplotypes/iter_001/gen_0002.rds
#   ...
#
# Each generation file contains a grid-sized list matrix. Every occupied cell
# contains the haplotype matrix returned by AlphaSimR::pullSegSiteHaplo().
#
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. User configuration
# ------------------------------------------------------------------------------

initial_isolation <- 10L        #(choose 10,100,300,500)
epistasis_type <- "bdm"         #(choose bdm or pos)
epistasis_magnitude <- "001"    #(choose 01,005,001)
preference_prob     <- 0.10     # choose 0.10, 0.50, 0.90
environment_type <- "gradient"  #(choose gradient or homogeneous)

config <- list( base_seed = 1223L,
  
  environment_file = file.path("data", "environments",
    sprintf("envi%d.RData", initial_isolation)),
  
  mating_script = file.path("R", "reproduction.R"),
  
  output_dir = file.path( "output",
    sprintf("isolation_%d_%s_%s_%s_%s", initial_isolation, epistasis_type,
      epistasis_magnitude,sprintf("%02d", round(preference_prob * 10)), environment_type)),
  
  n_iterations = 10L,
  preference_prob = preference_prob,
  epistasis_type = epistasis_type,
  epistasis_magnitude = epistasis_magnitude,
  environment_type = environment_type,
  cores_to_leave_free = 5L
)


# ------------------------------------------------------------------------------
# 2. Epistatic-effect vectors
# ------------------------------------------------------------------------------

EPI_EFFECTS <- list(
  "01" = c(-0.020733651, -0.018625038, -0.006291431, -0.004034495,
    -0.008043853, -0.005442798, -0.003021580, -0.015429728,
    -0.003002657, -0.011406430, -0.010789055, -0.015517526,
    -0.012286533, -0.008870254, -0.017678590),
  
  "005" = c(-0.0008622194, -0.0009184874, -0.0007796196, -0.0008385601,
    -0.0013560065, -0.0002860874, -0.0012618780, -0.0003284061,
    -0.0007399826, -0.0009307868, -0.0005277535, -0.0012485223,
    -0.0001434993, -0.0008188440, -0.0018042037),
  
  "001" = c(-0.004008850, -0.003938429, -0.004865292, -0.005568119,
    -0.005254795, -0.004558156, -0.004817784, -0.004825235,
    -0.005426624, -0.005845740, -0.004434927, -0.005463951,
    -0.005877341, -0.005674150, -0.005675997))

# ------------------------------------------------------------------------------
# 3. Package and input validation
# ------------------------------------------------------------------------------

check_packages <- function() {
  required_packages <- c("AlphaSimR", "future", "future.apply")
  
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  
  if (length(missing_packages) > 0L) {
    stop("Install the following packages before running the simulation: ",
      paste(missing_packages, collapse = ", "), call. = FALSE )
  }
  
  invisible(TRUE)
}

library(AlphaSimR)
library(future)
library(future.apply)


load_model_environment <- function(environment_path, mating_script) {
  if (!file.exists(environment_path)) {
    stop("Environment file not found: ",
      normalizePath(environment_path, mustWork = FALSE),
      call. = FALSE)
  }
  
  if (!file.exists(mating_script)) {
    stop("Mating script not found: ",
      normalizePath(mating_script, mustWork = FALSE),
      call. = FALSE)
  }
  
  # Source the mating code into the same environment as the model objects.
  # This preserves access to any environment objects used by mating3().
  model_env <- new.env(parent = globalenv())
  loaded_objects <- load(environment_path, envir = model_env)
  sys.source(mating_script, envir = model_env)
  
  required_objects <- c(
    "breedA", "breedB", "grid_size", "n_gen", "migration_rate",
    "nProgeny", "n_ind", "pairsE", "add_effA", "dom_effA",
    "add_effB", "dom_effB", "threshold", "mating3"
  )
  
  missing_objects <- required_objects[
    !vapply(required_objects, exists, logical(1), envir = model_env,
      inherits = FALSE)]
  
  sim_param_name <- if (exists("SP", envir = model_env, inherits = FALSE)) {
    "SP"
  } else if (exists("SimParam", envir = model_env, inherits = FALSE)) {
    "SimParam"
  } else {
    NA_character_
  }
  
  if (is.na(sim_param_name)) {
    missing_objects <- c(missing_objects, "SP (or SimParam)")
  }
  
  if (length(missing_objects) > 0L) {
    stop("The environment file is missing: ", paste(missing_objects, collapse = ", "),
      call. = FALSE)
  }
  
  model <- mget(required_objects, envir = model_env, inherits = FALSE)
  model$sim_param <- get(sim_param_name, envir = model_env, inherits = FALSE)
  model$loaded_objects <- loaded_objects
  
  model
}


validate_inputs <- function(config, model, epi_effect) {
  if (!config$environment_type %in% c("homogeneous", "gradient")) {
    stop("`environment_type` must be either 'homogeneous' or 'gradient'.",
      call. = FALSE)
  }
  
  if (!config$epistasis_type %in% c("bdm", "pos")) {
    stop("`epistasis_type` must be either 'bdm' or 'pos'.",call. = FALSE)
  }
  
  if (
    length(config$preference_prob) != 1L ||
    is.na(config$preference_prob) ||
    config$preference_prob < 0 ||
    config$preference_prob > 1
  ) {
    stop("`preference_prob` must be between 0 and 1.", call. = FALSE)
  }
  
  integer_parameters <- c(
    grid_size = model$grid_size,
    n_gen = model$n_gen,
    nProgeny = model$nProgeny,
    n_ind = model$n_ind,
    n_iterations = config$n_iterations
  )
  
  invalid_integer <- !is.finite(integer_parameters) |
    integer_parameters < 1 |
    integer_parameters %% 1 != 0
  
  if (any(invalid_integer)) {
    stop("These parameters must be positive integers: ",
      paste(names(integer_parameters)[invalid_integer], collapse = ", "),
      call. = FALSE)
  }
  
  if (model$grid_size < 2L) {
    stop("`grid_size` must be at least 2.", call. = FALSE)
  }
  
  if (length(model$migration_rate) != 1L ||
    is.na(model$migration_rate) ||
    model$migration_rate < 0 ||
    model$migration_rate > 1
  ) {
    stop("`migration_rate` must be between 0 and 1.", call. = FALSE)
  }
  
  if (!is.matrix(model$pairsE) || ncol(model$pairsE) != 2L) {
    stop("`pairsE` must be a two-column matrix.", call. = FALSE)
  }
  
  if (nrow(model$pairsE) != length(epi_effect)) {
    stop("The number of rows in `pairsE` must equal the number of ",
      "epistatic-effect values.", call. = FALSE)
  }
  
  if (!is.function(model$mating3)) {
    stop("`mating3` exists but is not a function.", call. = FALSE)
  }
  
  if (
    length(model$add_effA) != 100L ||
    length(model$dom_effA) != 100L ||
    length(model$add_effB) != 100L ||
    length(model$dom_effB) != 100L
  ) {
    stop("Each additive and dominance effect vector must contain 100 values.",
      call. = FALSE)
  }
  invisible(TRUE)
}


select_epistatic_effect <- function(type, magnitude) {
  if (!magnitude %in% names(EPI_EFFECTS)) {
    stop("Unknown epistasis magnitude. Available values: ",
      paste(names(EPI_EFFECTS), collapse = ", "),
      call. = FALSE)
  }
  
  base_effect <- EPI_EFFECTS[[magnitude]]
  
  if (identical(type, "bdm")) {
    abs(base_effect)
  } else {
        base_effect
  }
}


# ------------------------------------------------------------------------------
# 4. Population helpers
# ------------------------------------------------------------------------------

is_empty_population <- function(population) {
  is.null(population) || population@nInd < 1L
}


merge_nonempty_populations <- function(...) {
  populations <- list(...)
  
  keep <- vapply( populations,
    function(population) !is_empty_population(population),
    logical(1))
  
  populations <- populations[keep]
  
  if (length(populations) == 0L) {
    return(NULL)
  }
  
  if (length(populations) == 1L) {
    return(populations[[1L]])
  }
  
  AlphaSimR::mergePops(populations)
}


initialise_grid <- function(model) {
  grid_size <- as.integer(model$grid_size)
  
  demes <- matrix(vector("list", grid_size * grid_size),
    nrow = grid_size,  ncol = grid_size )
  
  for (x in seq_len(grid_size)) {
    for (y in seq_len(grid_size)) {
      demes[x, y] <- list(model$breedB)
    }
  }
  
  # For an even grid this gives the central 2 x 2 cells.
  # For an odd grid it gives the single central cell.
  centre_cells <- unique(c(
    floor((grid_size + 1) / 2),
    ceiling((grid_size + 1) / 2)
  ))
  
  for (x in centre_cells) {
    for (y in centre_cells) {
      demes[x, y] <- list(model$breedA)
    }
  }
  
  demes
}


# ------------------------------------------------------------------------------
# 5. Migration
# ------------------------------------------------------------------------------

get_neighbours <- function(x, y, grid_size) {
  neighbours <- list()
  
  if (x > 1L) {
    neighbours[[length(neighbours) + 1L]] <- c(x - 1L, y)
  }
  
  if (x < grid_size) {
    neighbours[[length(neighbours) + 1L]] <- c(x + 1L, y)
  }
  
  if (y > 1L) {
    neighbours[[length(neighbours) + 1L]] <- c(x, y - 1L)
  }
  
  if (y < grid_size) {
    neighbours[[length(neighbours) + 1L]] <- c(x, y + 1L)
  }
  
  neighbours
}


migrate_generation <- function(demes, grid_size, migration_rate) {
  new_demes <- matrix(
    vector("list", grid_size * grid_size),
    nrow = grid_size,
    ncol = grid_size
  )
  
  for (x in seq_len(grid_size)) {
    for (y in seq_len(grid_size)) {
      deme <- demes[[x, y]]
      
      if (!is_empty_population(deme) && deme@nInd > 15L) {
        n_migrants <- min(floor(deme@nInd * migration_rate), deme@nInd - 2L)
        
        if (n_migrants > 0L) {
          migrant_indices <- sample.int(deme@nInd, size = n_migrants,
            replace = FALSE)
          
          migrants <- deme[migrant_indices]
          migrants <- migrants[sample.int(migrants@nInd)]
          
          neighbours <- get_neighbours(x, y, grid_size)
          
          if (length(neighbours) > 0L) {
            migrants_per_neighbour <- ceiling(n_migrants / length(neighbours))
            
            migrant_start <- 1L
            
            for (neighbour_id in sample(seq_along(neighbours))) {
              if (migrant_start > n_migrants) {
                break
              }
              
              migrant_end <- min(migrant_start + migrants_per_neighbour - 1L,
                n_migrants)
              
              migrant_subset <- migrants[migrant_start:migrant_end]
              
              nx <- neighbours[[neighbour_id]][1L]
              ny <- neighbours[[neighbour_id]][2L]
              
              new_demes[nx, ny] <- list(merge_nonempty_populations(
                  new_demes[[nx, ny]], migrant_subset))
              
              migrant_start <- migrant_end + 1L
            }
          }
          
          remaining_indices <- setdiff( seq_len(deme@nInd), migrant_indices)
          
          deme <- if (length(remaining_indices) > 2L) {
            deme[remaining_indices]
          } else {
            NULL
          }
        }
      }
      
      new_demes[x, y] <- list(
        merge_nonempty_populations(new_demes[[x, y]], deme)
      )
    }
  }
  
  new_demes
}


# ------------------------------------------------------------------------------
# 6. Reproduction and selection
# ------------------------------------------------------------------------------

reproduce_population <- function(deme,preference_prob,n_progeny,threshold,
    mating_function,sim_param) {
  if (is_empty_population(deme) || deme@nInd <= 2L) {
    return(NULL)
  }
  use_assortative_mating <- runif(1) <= preference_prob
  
  if (use_assortative_mating) { 
    mating_function(deme, threshold, n_progeny, sim_param )
  } else {
    AlphaSimR::randCross(deme, nCrosses = floor(deme@nInd / 2), nProgeny = n_progeny,
      simParam = sim_param)
  }
}


calculate_fitness <- function(deme, x, grid_size, environment_type, pairsE,
    epi_effect,  add_effA, dom_effA, add_effB, dom_effB, sim_param) {
  genotypes <- AlphaSimR::pullSegSiteGeno(deme, simParam = sim_param)
  
  if (ncol(genotypes) < 230L) {
    stop("At least 230 segregating-site columns are required; found ",
      ncol(genotypes),".", call. = FALSE)
  }
  
  genoB <- genotypes[, 1:100, drop = FALSE]
  genoA <- genotypes[, 101:200, drop = FALSE]
  genoE <- genotypes[, 201:230, drop = FALSE] - 1L
  
  p1 <- pairsE[, 1L]
  p2 <- pairsE[, 2L]
  
  if (
    any(p1 < 1L) ||
    any(p2 < 1L) ||
    any(p1 > ncol(genoE)) ||
    any(p2 > ncol(genoE))
  ) {
    stop("`pairsE` contains an index outside the 30 epistatic loci.",
      call. = FALSE)
  }
  
  pair_products <-
    genoE[, p1, drop = FALSE] *  genoE[, p2, drop = FALSE]
  
  epistatic_value <- drop(pair_products %*% epi_effect)
  
  parental_value_A <- drop((genoA - 1L) %*% add_effA + (genoA == 1L) %*% dom_effA)
  
  parental_value_B <- drop((genoB - 1L) %*% add_effB + (genoB == 1L) %*% dom_effB)
  
  if (identical(environment_type, "homogeneous")) {
    return(parental_value_B + parental_value_A + epistatic_value)
  }
  
  gradient_weight <- (x - 1L) / (grid_size - 1L)
  
  ((1 - gradient_weight) * parental_value_B) +
    (gradient_weight * parental_value_A) + epistatic_value
}


select_population <- function(deme, fitness, maximum_survivors) {
  if (is_empty_population(deme)) {
    return(NULL)
  }
  n_keep <- max(
    1L, min(deme@nInd, as.integer(floor(maximum_survivors))))
  
  selected_indices <- order(fitness, decreasing = TRUE)[seq_len(n_keep)]
  deme[selected_indices]
}


# ------------------------------------------------------------------------------
# 7. One simulation iteration
# ------------------------------------------------------------------------------

run_one_iteration <- function(iteration, model, config, epi_effect, output_dir) {
  had_global_sp <- exists("SP", envir = .GlobalEnv, inherits = FALSE)
  had_global_simparam <- exists("SimParam", envir = .GlobalEnv, inherits = FALSE)
  
  if (had_global_sp) {
    previous_global_sp <- get("SP", envir = .GlobalEnv, inherits = FALSE)
  }
  
  if (had_global_simparam) {previous_global_simparam <- get("SimParam",
      envir = .GlobalEnv, inherits = FALSE)
  }
  
  assign("SP", model$sim_param, envir = .GlobalEnv)
  assign("SimParam", model$sim_param, envir = .GlobalEnv)
  
  on.exit(
    {
      if (had_global_sp) {
        assign("SP", previous_global_sp, envir = .GlobalEnv)
      } else if (exists("SP", envir = .GlobalEnv, inherits = FALSE)) {
        rm("SP", envir = .GlobalEnv)
      }
      
      if (had_global_simparam) {
        assign("SimParam", previous_global_simparam, envir = .GlobalEnv)
      } else if (exists("SimParam", envir = .GlobalEnv, inherits = FALSE)) {
        rm("SimParam", envir = .GlobalEnv)
      }
    },
    add = TRUE
  )
  iteration_dir <- file.path(output_dir, sprintf("iter_%03d", iteration)
  )
  
  dir.create(iteration_dir, showWarnings = FALSE, recursive = TRUE)
  
  grid_size <- as.integer(model$grid_size)
  n_generations <- as.integer(model$n_gen)
  
  demes <- initialise_grid(model)
  
  for (generation in seq_len(n_generations)) {
    message(sprintf("Iteration %d/%d | Generation %d/%d", iteration,
        config$n_iterations, generation, n_generations))
    
    demes <- migrate_generation(demes = demes, grid_size = grid_size,
      migration_rate = model$migration_rate)
    
    demes_next <- matrix(vector("list", grid_size * grid_size),
      nrow = grid_size, ncol = grid_size)
    
    haplotypes <- matrix(vector("list", grid_size * grid_size), nrow = grid_size,
      ncol = grid_size)
    
    for (x in seq_len(grid_size)) {
      for (y in seq_len(grid_size)) {
        deme <- demes[[x, y]]
        
        offspring <- reproduce_population(
          deme = deme,
          preference_prob = config$preference_prob,
          n_progeny = model$nProgeny,
          threshold = model$threshold,
          mating_function = model$mating3,
          sim_param = model$sim_param
        )
        
        if (is_empty_population(offspring)) {
          demes_next[x, y] <- list(NULL)
          haplotypes[x, y] <- list(NULL)
          next
        }
        
        fitness <- calculate_fitness(
          deme = offspring,
          x = x,
          grid_size = grid_size,
          environment_type = config$environment_type,
          pairsE = model$pairsE,
          epi_effect = epi_effect,
          add_effA = model$add_effA,
          dom_effA = model$dom_effA,
          add_effB = model$add_effB,
          dom_effB = model$dom_effB,
          sim_param = model$sim_param
        )
        
        selected_deme <- select_population(deme = offspring, fitness = fitness,
          maximum_survivors = model$n_ind / 2)
        
        haplotypes[x, y] <- list(AlphaSimR::pullSegSiteHaplo(selected_deme,
            simParam = model$sim_param))
        
        demes_next[x, y] <- list(selected_deme)
      }
    }
    
    generation_file <- file.path(iteration_dir, sprintf("gen_%04d.rds", generation))
    
    saveRDS(haplotypes, generation_file, compress = "gzip")
    
    demes <- demes_next
    
    rm(haplotypes)
    
    if (generation %% 5L == 0L) {
      gc(verbose = FALSE)
    }
  }
  
  list(
    iteration = iteration,
    directory = iteration_dir
  )
}


# ------------------------------------------------------------------------------
# 8. Parallel execution
# ------------------------------------------------------------------------------

run_parallel_simulation <- function(model, config, epi_effect, output_dir) {
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  
  thread_variables <- c(
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS"
  )
  
  old_thread_values <- Sys.getenv(
    thread_variables,
    unset = NA_character_
  )
  
  on.exit(
    {
      for (i in seq_along(thread_variables)) {
        variable <- thread_variables[[i]]
        old_value <- old_thread_values[[i]]
        
        if (is.na(old_value)) {
          Sys.unsetenv(variable)
        } else {
          do.call(Sys.setenv, setNames(list(old_value), variable))
        }
      }
    },
    add = TRUE
  )
  
  Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
  
  available_cores <- parallel::detectCores(logical = TRUE)
  
  if (is.na(available_cores)) {
    available_cores <- 1L
  }
  
  n_workers <- max(1L, available_cores - as.integer(config$cores_to_leave_free))
  
  message("Parallel workers: ", n_workers)
  
  future::plan(future::multisession, workers = n_workers)
  
  future.apply::future_lapply(
    X = seq_len(config$n_iterations),
    FUN = run_one_iteration,
    model = model,
    config = config,
    epi_effect = epi_effect,
    output_dir = output_dir,
    future.seed = config$base_seed
  )
}


# ------------------------------------------------------------------------------
# 9. Main entry point
# ------------------------------------------------------------------------------

run_hybgrid_simulation <- function(config) {
  check_packages()
  
  RNGkind(kind = "L'Ecuyer-CMRG", normal.kind = "Inversion",
          sample.kind = "Rejection")
  
  set.seed(config$base_seed)
  
  model <- load_model_environment(environment_path = config$environment_file,
    mating_script = config$mating_script)
  
  epi_effect <- select_epistatic_effect(type = config$epistasis_type,
    magnitude = config$epistasis_magnitude)
  
  validate_inputs(config = config, model = model, epi_effect = epi_effect)
  
  output_dir <- normalizePath(config$output_dir, mustWork = FALSE)
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  metadata <- list(created_at = Sys.time(), configuration = config,
    environment_file = normalizePath(config$environment_file, mustWork = TRUE),
    mating_script = normalizePath(config$mating_script, mustWork = TRUE),
    loaded_environment_objects = model$loaded_objects,
    epistatic_effect = epi_effect,
    model_dimensions = list(grid_size = model$grid_size,
      n_generations = model$n_gen, initial_population_size = model$n_ind),
    session_info = utils::sessionInfo() )
  
  saveRDS(metadata, file.path(output_dir, "simulation_metadata.rds"))
  
  results <- run_parallel_simulation(model = model, config = config,
    epi_effect = epi_effect, output_dir = output_dir)
  
  saveRDS(results, file.path(output_dir, "iteration_index.rds"))
  
  invisible(results)
}


# Run only when this file is executed directly.
if (sys.nframe() == 0L) {
  simulation_results <- run_hybgrid_simulation(config)
}
