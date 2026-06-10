# Optimizing SKU-Warehouse Allocations to minimize Split Parcels in E-Commerce Environments

## About the Project
This repository contains the implementation of novel heuristics for minimizing split deliveries in e-commerce warehouse allocation problems, as presented in our research article. We introduce two new approaches:
- QMK heuristic
- CHI heuristic

Additionally, we implement six state-of-the-art competing heuristics for comparison. This codebase enables full reproduction of our research results and can be adapted for any split-delivery minimization problem.

## Implemented Heuristics

### Novel Approaches
* **QMK**: QMK heuristic optimization model (Gurobi solver)
* **QMKS**: QMK heuristic optimization model (SCIP solver)
* **QMKJ**: QMK heuristic optimization model (Juniper solver)
* **CHI**: Chi-square heuristic

### Benchmark Heuristics
* **KL**: K-LINK heuristic ([Zhu et al., 2021](https://doi.org/10.1016/j.ejor.2020.08.024))
* **GO**: GREEDY ORDERS heuristic ([Catalan & Fisher, 2012](https://doi.org/10.2139/ssrn.2166687))
* **GP**: GREEDY PAIRS heuristic ([Catalan & Fisher, 2012](https://doi.org/10.2139/ssrn.2166687))
* **GS**: GREEDY SEEDS heuristic ([Catalan & Fisher, 2012](https://doi.org/10.2139/ssrn.2166687))
* **BS**: BESTSELLERS heuristic ([Catalan & Fisher, 2012](https://doi.org/10.2139/ssrn.2166687))
* **EMCI**: Extended MCI heuristic for D warehouses ([Lin et al., 2025](https://doi.org/10.1111/poms.14114))
* **IIH/IIHS**: Iterative Improvement Heuristic for 2 overlapping warehouses ([Lin et al., 2025](https://doi.org/10.1111/poms.14114)), with Gurobi or SCIP solver
* **OPT**: Exact optimization model (Gurobi)
* **RND**: Random allocation baseline

Configuration of heuristics and parameters can be adjusted in `main_benchmark_settings.jl`.

## Required Solvers
* [Gurobi](https://www.gurobi.com) - Required for QMK, KLQ, IIH, and OPT modes
* [SCIP](https://www.scipopt.org) - Required for QMKS and IIHS modes (installed via Julia package)
* Juniper (with HiGHS/Ipopt) - Required for QMKJ mode (no license needed)

## Getting Started

### Prerequisites
- Julia installation
- Valid Gurobi license (required for QMK, KLQ, IIH, and OPT heuristics)

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

Alternatively, run all six dependency variants in parallel via tmux:
```bash
bash run_benchmarks.sh
```

## Dependencies
The following Julia packages are required:
* [CSV.jl](https://github.com/JuliaData/CSV.jl) - CSV file handling
* [DataFrames.jl](https://github.com/JuliaData/DataFrames.jl) - Data manipulation
* [Combinatorics.jl](https://github.com/JuliaMath/Combinatorics.jl) - Combinatorial algorithms
* [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) - Probability distributions
* [JuMP.jl](https://github.com/jump-dev/JuMP.jl) - Mathematical optimization
* [Gurobi.jl](https://github.com/jump-dev/Gurobi.jl) - Gurobi solver interface
* [SCIP.jl](https://github.com/scipopt/SCIP.jl) - SCIP solver interface
* [HiGHS.jl](https://github.com/jump-dev/HiGHS.jl) - HiGHS solver interface
* [Ipopt.jl](https://github.com/jump-dev/Ipopt.jl) - Ipopt solver interface
* [Juniper.jl](https://github.com/lanl-ansi/Juniper.jl) - Nonlinear integer programming
* [LoopVectorization.jl](https://github.com/JuliaSIMD/LoopVectorization.jl) - Loop optimization

## License
Distributed under the MIT License. See `LICENSE.txt` for more information.

## Contact
- Project Link: [https://github.com/beyondsimulations/splitdeliveries](https://github.com/beyondsimulations/splitdeliveries)
- For issues and feature requests, please use the [GitHub Issues](https://github.com/beyondsimulations/splitdeliveries/issues) page

## References
* Lin, Y.-H., Zhu, S., Wang, R. (2025). Multi-Warehouse Assortment Selection: Minimizing Order Splitting in E-Commerce Logistics. Production and Operations Management. https://doi.org/10.1111/poms.14114
* Zhu, S., Hu, X., Huang, K. et al. (2021). Optimization of product category allocation in multiple warehouses to minimize splitting of online supermarket customer orders. European Journal of Operational Research. https://doi.org/10.1016/j.ejor.2020.08.024
* Catalan, A., Fisher, M. (2012). Assortment Allocation to Distribution Centers to Minimize Split Customer Orders. SSRN. https://doi.org/10.2139/ssrn.2166687
