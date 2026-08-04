# HybGrid-Evo

HybGrid-Evo is an R-based framework for spatially explicit simulations of hybridization between two differentiated parental populations.

The model includes:

- migration across a two-dimensional grid;
- random and assortative mating;
- additive, dominance, and pairwise epistatic effects;
- homogeneous or gradient environment;
- alternative initial isolation durations;
- storage of haplotypes for downstream analyses.

## Repository structure

```text
HybGrid-Evo/
├── R/
│   ├── initial_populations.R
│   ├── reproduction.R
│   └── simulation.R
├── data/
│   └── environments/
│       ├── envi10.RData
│       ├── envi100.RData
│       ├── envi300.RData
│       └── envi500.RData
└── output/
```

## Scripts

### `R/initial_populations.R`

Creates the initial parental populations, defines the genetic architecture, performs burn-in selection, and saves prepared AlphaSimR environments.

### `R/reproduction.R`

Defines the assortative-mating function. Individuals are randomly paired and compatible pairs are identified using genetic distance across population-informative loci.

### `R/simulation.R`

Runs the spatial hybridization simulation, including migration, reproduction, fitness calculation, selection, parallel iterations, and haplotype output.

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
