# ==============================================================================
# HETEROZYGOSITY
# ==============================================================================

## SETTINGS --------------------------------------------------------------------

project_dir <- ".../HybGrid-Evo"

initial_isolation   <- 10L        # choose 10, 100, 300, 500
epistasis_type      <- "bdm"      # choose "bdm" or "pos"
epistasis_magnitude <- "001"      # choose "01", "005", "001"
preference_prob     <- 0.10       # choose 0.10, 0.50, 0.90
environment_type    <- "gradient" # choose "gradient" or "homogeneous"

n_iterations <- 10L
n_generations <- 50L
scenario_name <- sprintf("isolation_%d_%s_%s_%s_%s", initial_isolation, epistasis_type,
                         epistasis_magnitude,sprintf("%02d", round(preference_prob * 10)), environment_type)
source_folder <- file.path(project_dir, "output",scenario_name)
save_folder   <- source_folder

## LOCI -------------------------------------------------------------------------

loci_qtl   <- 1:200
loci_class <- 201:280
loci_neu   <- 281:300

pop_levels <- c("A", "B", "H")
needed_cols <- sort(unique(c(loci_qtl,loci_class, loci_neu)))
idx_qtl   <- match(loci_qtl, needed_cols)
idx_class <- match(loci_class, needed_cols)
idx_neu   <- match(loci_neu, needed_cols)

## OUTPUT MATRICES --------------------------------------------------------------

qtl_A <- qtl_B <- qtl_H <- matrix(NA_real_, nrow = n_iterations, ncol = n_generations)
neu_A <- neu_B <- neu_H <- matrix(NA_real_, nrow = n_iterations, ncol = n_generations)

## CALCULATE HETEROZYGOSITY -----------------------------------------------------

for (i in seq_len(n_iterations)) {
  
  cat("Iteration", i, "/", n_iterations, "\n")
  
  for (j in seq_len(n_generations)) {
    
    file <- file.path(source_folder,sprintf("iter_%03d", i),sprintf("gen_%04d.rds", j))
    
    if (!file.exists(file)) next
    haps <- readRDS(file)
    ## Store individual heterozygosities from all occupied cells
    qtl_list <- list(A = numeric(0), B = numeric(0), H = numeric(0))
    
    neu_list <- list(A = numeric(0), B = numeric(0), H = numeric(0))
    
    for (hap in haps) {
      if (is.null(hap) || nrow(hap) < 2L) next
      if (ncol(hap) < max(needed_cols)) {
        stop(sprintf("Only %d loci found in iteration %d, generation %d; %d required.",
          ncol(hap),i,j,max(needed_cols)))
      }
      ## Protect against incomplete final haplotype pair
      odd  <- seq.int(1L, nrow(hap) - 1L, by = 2L)
      even <- odd + 1L
      ## Only columns actually needed
      gens <- hap[odd, needed_cols, drop = FALSE] +  hap[even, needed_cols, drop = FALSE]
      
      ## Population classification ---------------------------------------------
      
      cls <- rowMeans(gens[, idx_class, drop = FALSE]) / 2
      ids <- list(
        A = cls >= 0.9,
        B = cls <= 0.1,
        H = cls > 0.1 & cls < 0.9)
      
      ## Individual heterozygosity ---------------------------------------------
      
      het_qtl <- rowMeans(gens[, idx_qtl, drop = FALSE] == 1, na.rm = TRUE)
      het_neu <- rowMeans(gens[, idx_neu, drop = FALSE] == 1, na.rm = TRUE)
      
      ## Pool individuals by population ----------------------------------------
      
      for (pop in pop_levels) {
        
        id <- ids[[pop]]
        if (!any(id)) next
        qtl_list[[pop]] <- c(qtl_list[[pop]], het_qtl[id])
        neu_list[[pop]] <- c(neu_list[[pop]], het_neu[id])
      }
    }
    
    ## Mean heterozygosity across individuals ----------------------------------
    
    if (length(qtl_list$A))
      qtl_A[i, j] <- mean(qtl_list$A, na.rm = TRUE)
    if (length(qtl_list$B))
      qtl_B[i, j] <- mean(qtl_list$B, na.rm = TRUE)
    if (length(qtl_list$H))
      qtl_H[i, j] <- mean(qtl_list$H, na.rm = TRUE)
    if (length(neu_list$A))
      neu_A[i, j] <- mean(neu_list$A, na.rm = TRUE)
    if (length(neu_list$B))
      neu_B[i, j] <- mean(neu_list$B, na.rm = TRUE)
    if (length(neu_list$H))
      neu_H[i, j] <- mean(neu_list$H, na.rm = TRUE)
  }
}

## SAVE -------------------------------------------------------------------------
H <- list(QTL = list( A = qtl_A, B = qtl_B, H = qtl_H), Neutral = list(A = neu_A,
    B = neu_B, H = neu_H))

output_file <- file.path(save_folder, "heterozygosity.rds")
saveRDS(H, output_file, compress = TRUE)
