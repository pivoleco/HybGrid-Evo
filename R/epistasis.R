# ==============================================================================
# EPISTATIC LOCUS METRICS
# ==============================================================================

## SETTINGS --------------------------------------------------------------------

project_dir <- ".../HybGrid-Evo"

initial_isolation   <- 10L        # choose 10, 100, 300, 500
epistasis_type      <- "bdm"      # choose "bdm" or "pos"
epistasis_magnitude <- "001"      # choose "01", "005", "001"
preference_prob     <- 0.10       # choose 0.10, 0.50, 0.90
environment_type    <- "gradient" # choose "gradient" or "homogeneous"

n_iterations  <- 10L
n_generations <- 50L
grid_size     <- 20L
max_ind       <- 50L
scenario_name <- sprintf("isolation_%d_%s_%s_%s_%s", initial_isolation, epistasis_type,
                         epistasis_magnitude,sprintf("%02d", round(preference_prob * 10)), environment_type)
folder <- file.path(project_dir, "output", scenario_name)

## LOCI -------------------------------------------------------------------------

loci_epi   <- 201:230
loci_class <- 201:280

needed_cols <- sort(unique(c(
  loci_epi,
  loci_class
)))

idx_epi   <- match(loci_epi, needed_cols)
idx_class <- match(loci_class, needed_cols)

## OUTPUT ARRAY ----------------------------------------------------------------

epi_gen_pop <- array(
  NA_real_,
  dim = c(
    n_iterations,
    n_generations,
    grid_size,
    grid_size,
    max_ind,
    3,
    3
  ),
  dimnames = list(
    Iteration = NULL,
    Generation = NULL,
    X = NULL,
    Y = NULL,
    Individual = NULL,
    Metric = c(
      "epi_mean",
      "prop_22",
      "prop_00"
    ),
    Pop = c("A", "B", "H")
  )
)

## CALCULATE -------------------------------------------------------------------

for (i in seq_len(n_iterations)) {
  
  message(
    scenario_name,
    " | iteration ",
    i
  )
  
  for (j in seq_len(n_generations)) {
    
    file <- file.path(
      folder,
      sprintf("iter_%03d", i),
      sprintf("gen_%04d.rds", j)
    )
    
    if (!file.exists(file)) next
    
    haps <- readRDS(file)
    
    for (x in seq_len(grid_size)) {
      for (y in seq_len(grid_size)) {
        
        hap <- haps[[x, y]]
        
        if (is.null(hap) || nrow(hap) < 2L) next
        
        if (ncol(hap) < max(needed_cols)) {
          stop(sprintf(
            "Only %d loci found in iter %d, gen %d, deme [%d,%d]; %d required.",
            ncol(hap),
            i,
            j,
            x,
            y,
            max(needed_cols)
          ))
        }
        
        ## Protect against incomplete final haplotype pair
        odd <- seq.int(
          1L,
          nrow(hap) - 1L,
          by = 2L
        )
        
        even <- odd + 1L
        
        ## Only required loci
        gens <- hap[
          odd,
          needed_cols,
          drop = FALSE
        ] +
          hap[
            even,
            needed_cols,
            drop = FALSE
          ]
        
        
        ## Population classification -------------------------------------------
        
        cls <- rowMeans(
          gens[, idx_class, drop = FALSE]
        ) / 2
        
        ids <- list(
          A = cls >= 0.9,
          B = cls <= 0.1,
          H = cls > 0.1 & cls < 0.9
        )
        
        
        ## Epistatic-locus metrics ---------------------------------------------
        
        geno_epi <- gens[
          ,
          idx_epi,
          drop = FALSE
        ]
        
        for (pop in names(ids)) {
          
          idx <- ids[[pop]]
          n_ind <- sum(idx)
          
          if (n_ind == 0L) next
          
          geno_sub <- geno_epi[
            idx,
            ,
            drop = FALSE
          ]
          
          vals <- cbind(
            
            # Mean allele dosage at epistatic loci:
            # 0 = all genotype 0
            # 0.5 = average genotype 1
            # 1 = all genotype 2
            epi_mean = rowMeans(
              geno_sub,
              na.rm = TRUE
            ) / 2,
            
            # Proportion of epistatic loci homozygous "2"
            prop_22 = rowMeans(
              geno_sub == 2,
              na.rm = TRUE
            ),
            
            # Proportion of epistatic loci homozygous "0"
            prop_00 = rowMeans(
              geno_sub == 0,
              na.rm = TRUE
            )
          )
          
          epi_gen_pop[
            i,
            j,
            x,
            y,
            seq_len(n_ind),
            ,
            pop
          ] <- vals
        }
      }
    }
  }
}


## SAVE ------------------------------------------------------------------------

output_file <- file.path( folder, "epistatic_metrics.rds")
saveRDS(epi_gen_pop,output_file,compress = TRUE)

