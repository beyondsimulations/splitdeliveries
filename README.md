# Optimizing SKU-Warehouse Allocations to minimize Split Parcels in E-Commerce Environments

## About the Project
This repository contains the implementation of novel heuristics for minimizing split deliveries in e-commerce warehouse allocation problems, as presented in our research article. We introduce two new approaches:
- QMK heuristic
- CHI heuristic

Additionally, we implement four state-of-the-art competing heuristics for comparison. This codebase enables full reproduction of our research results and can be adapted for any split-delivery minimization problem.

## Built With
* [Julia](https://julialang.org/) - Primary programming language
* [GAMS](https://www.gams.com) - Mathematical optimization system

## Implemented Heuristics

### Novel Approaches
* **QMKOPT**: QMK heuristic optimization model (CPLEX solver)
* **QMK**: QMK heuristic optimization model (BONMIN solver)
* **CHISOL**: CHI heuristic without local search
* **CHI**: CHI heuristic with local search

### Benchmark Heuristics
* **KLINK**: K-LINK heuristic ([Zhu et al., 2021](https://doi.org/10.1016/j.ejor.2020.08.024))
* **GP**: GREEDY PAIRS heuristic ([Catalan & Fisher, 2012](https://doi.org/10.2139/ssrn.2166687))
* **GS**: GREEDY SEEDS heuristic ([Catalan & Fisher, 2012](https://doi.org/10.2139/ssrn.2166687))
* **BS**: BESTSELLERS heuristic ([Catalan & Fisher, 2012](https://doi.org/10.2139/ssrn.2166687))
* **OPT**: CPLEX optimization model based on Catalan & Fisher (2012) and Zhu et al. (2021)
* **RND**: Random allocation baseline

Configuration of heuristics and parameters can be adjusted in `main_benchmark_settings.jl`.

## Required Solvers
* [CPLEX](https://www.ibm.com/analytics/cplex-optimizer) - For QMKOPT and OPT
* [SBB](https://www.gams.com/latest/docs/S_SBB.html) - GAMS solver

## Getting Started

### Prerequisites
- Julia installation
- GAMS installation
- Valid GAMS license (required for QMKOPT, QMK, and OPT heuristics)

### Installation
1. Install Julia and GAMS
2. Clone this repository:
   ```bash
   git clone [repository-url]
   ```
3. Install required packages:
   ```bash
   julia install_packages.jl
   ```

### Running Benchmarks
1. Configure parameters in `main_benchmark_settings.jl`
2. Execute the benchmark:
   ```bash
   julia main_benchmark_settings.jl
   ```
3. Results will be saved to the `results` directory

## Dependencies
The following Julia packages are required:
* [CSV.jl](https://github.com/JuliaData/CSV.jl) - CSV file handling
* [DataFrames.jl](https://github.com/JuliaData/DataFrames.jl) - Data manipulation
* [Combinatorics.jl](https://github.com/JuliaMath/Combinatorics.jl) - Combinatorial algorithms
* [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) - Probability distributions
* [GAMS.jl](https://github.com/JuliaData/GAMS.jl) - GAMS interface
* [JuMP.jl](https://github.com/jump-dev/JuMP.jl) - Mathematical optimization
* [Octavian.jl](https://github.com/JuliaLinearAlgebra/Octavian.jl) - Matrix operations
* [LoopVectorization.jl](https://github.com/JuliaSIMD/LoopVectorization.jl) - Loop optimization

## License
Distributed under the MIT License. See `LICENSE.txt` for more information.

## References
* Zhu, S., Hu, X., Huang, K. et al. (2021). Optimization of product category allocation in multiple warehouses to minimize splitting of online supermarket customer orders. European Journal of Operational Research. https://doi.org/10.1016/j.ejor.2020.08.024
* Catalan, A., Fisher, M. (2012). Assortment Allocation to Distribution Centers to Minimize Split Customer Orders. SSRN. https://doi.org/10.2139/ssrn.2166687
