# Optimizing SKU-Warehouse Allocations to minimize Split Parcels in E-Commerce Environments

## About the Project
This repository contains the implementation of novel heuristics for minimizing split deliveries in e-commerce warehouse allocation problems, as presented in our research article. We introduce two new approaches:
- QMK heuristic
- CHI heuristic

Additionally, we implement four state-of-the-art competing heuristics for comparison. This codebase enables full reproduction of our research results and can be adapted for any split-delivery minimization problem.

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
* [Gurobi](https://www.gurobi.com) - Mathematical optimization solver

## Getting Started

### Prerequisites
- Julia installation
- Valid Gurobi license (required for QMKOPT, QMK, and OPT heuristics)

### Installation
1. Install Julia from [julialang.org](https://julialang.org/downloads/)
2. Install Gurobi:
   - Download and install Gurobi from [gurobi.com](https://www.gurobi.com/downloads/)
   - Set up your license 
3. Clone this repository:
   ```bash
   git clone https://github.com/beyondsimulations/splitdeliveries
   cd splitdeliveries
   ```
4. Install required packages:
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
* [Gurobi.jl](https://github.com/jump-dev/Gurobi.jl) - Mathematical Solver
* [JuMP.jl](https://github.com/jump-dev/JuMP.jl) - Mathematical optimization
* [Octavian.jl](https://github.com/JuliaLinearAlgebra/Octavian.jl) - Matrix operations
* [LoopVectorization.jl](https://github.com/JuliaSIMD/LoopVectorization.jl) - Loop optimization

## License
Distributed under the MIT License. See `LICENSE.txt` for more information.

## Contact
- Project Link: [https://github.com/beyondsimulations/splitdeliveries](https://github.com/beyondsimulations/splitdeliveries)
- For issues and feature requests, please use the [GitHub Issues](https://github.com/beyondsimulations/splitdeliveries/issues) page

## References
* Zhu, S., Hu, X., Huang, K. et al. (2021). Optimization of product category allocation in multiple warehouses to minimize splitting of online supermarket customer orders. European Journal of Operational Research. https://doi.org/10.1016/j.ejor.2020.08.024
* Catalan, A., Fisher, M. (2012). Assortment Allocation to Distribution Centers to Minimize Split Customer Orders. SSRN. https://doi.org/10.2139/ssrn.2166687
