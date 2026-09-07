# HybGrid-Evo

HybGrid-Evo is an R-based framework for spatially explicit, individual-based simulations of hybridization following secondary contact between genetically differentiated native and introduced populations.

The model was developed to investigate how genetic architecture, pre-contact divergence, assortative mating, and environmental structure influence demographic and evolutionary outcomes following biological invasions.

The model includes:

- migration across a two-dimensional spatial grid;
- random and assortative mating;
- additive, dominance, and pairwise epistatic genetic effects;
- advantageous and Bateson–Dobzhansky–Muller-like (BDM) epistatic interactions;
- homogeneous or spatially heterogeneous environments;
- alternative durations of pre-contact divergence;
- storage of haplotypes through time;
- post-processing of genotypic values, genetic variance, heterozygosity, and epistatic loci.

## Repository structure

```text
HybGrid-Evo/
├── R/
│   ├── initial_populations.R
│   ├── reproduction.R
│   ├── simulation.R
│   ├── genotypic_values.R
│   ├── variances.R
│   ├── heterozygosity.R
│   └── epistasis.R
│
├── data/
│   └── environments/
│       ├── envi10.RData
│       ├── envi100.RData
│       ├── envi300.RData
│       └── envi500.RData
│
└── output/
```

## Simulation scripts

### `R/initial_populations.R`

Creates the initial parental populations and genetic architecture, performs the period of pre-contact divergence, and saves prepared AlphaSimR environments used as starting conditions for the spatial simulations.

### `R/reproduction.R`

Defines the assortative-mating function. Individuals are randomly paired and compatible pairs are identified using genetic distance across population-informative loci.

### `R/simulation.R`

Runs the spatial hybridization simulation, including migration, reproduction, fitness calculation, selection, parallel iterations, and haplotype output.

## Post-processing scripts

### `R/genotypic_values.R`

Calculates genotypic values and their additive, dominance, and epistatic components from stored simulation output.

### `R/variances.R`

Calculates genetic variance components from simulated populations, including additive, dominance, and total genotypic variance.

### `R/heterozygosity.R`

Calculates population-level heterozygosity through time for native, introduced, and hybrid populations.

### `R/epistasis.R`

Calculates metrics associated with epistatic loci and their representation within native, introduced, and hybrid populations through time.

## Prepared environments

The `data/environments/` directory contains four prepared starting environments corresponding to different initial isolation durations:

| File | Initial isolation duration |
|---|---:|
| `envi10.RData` | 10 generations |
| `envi100.RData` | 100 generations |
| `envi300.RData` | 300 generations |
| `envi500.RData` | 500 generations |

## Requirements

The simulations require R and the following packages:

```r
install.packages(c(
  "future",
  "future.apply"
))
```

AlphaSimR must also be installed.

## Running a simulation

Open `R/simulation.R` and choose the initial isolation duration:

```r
initial_isolation <- 10L
```

Available values are:

```r
c(10L, 100L, 300L, 500L)
```

The simulation output is written to the `output/` directory.

## Output

Each simulation iteration is stored in a separate folder. One compressed RDS file is created for each generation, containing the grid of saved haplotypes.

Generated simulation outputs are not tracked by Git because of their size.

## Status

This repository contains research code associated with an ongoing manuscript. The repository structure and documentation may be updated during manuscript preparation.
