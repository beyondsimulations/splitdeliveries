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
* **EMCI**: Extended MCI heuristic for multiple warehouses ([Lin et al., 2025](https://doi.org/10.1177/10591478251365581))
* **IIH/IIHS**: Iterative Improvement Heuristic for 2 overlapping warehouses ([Lin et al., 2025](https://doi.org/10.1177/10591478251365581)), with Gurobi or SCIP solver
* **OPT**: Exact optimization model (Gurobi)
* **RND**: Random allocation baseline

Configuration of heuristics and parameters can be adjusted in `main_benchmark_settings.jl`.

## Repository Structure
* `main_benchmark.jl` / `main_benchmark_settings.jl` - benchmark framework and its configuration
* `heuristics/`, `functions/` - algorithm implementations and shared routines
* `benchmark_generator.jl` - synthetic transaction generator (dependency structure + transaction sampling)
* `dependency/`, `capacity/` - dataset and warehouse-constellation definitions of the numerical study
* `results/` - benchmark outputs of the main runs; `results_chigate/` - CHI re-run with the independence gate enabled
* `aggregate_benchmarks.jl` - merges all runs into `results/overall_results.csv`
* `evaluation/` - scripts that generate the result tables of the article
* `helpers/` - small utilities (generator verification, test data)

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
`main_benchmark_settings.jl` holds all parameters and expects three variables to be defined before it is included: the dataset configurations, the SKU scale, and the storage-requirement mode. A single run looks like this:
```julia
dependencies = ["HD-VF"]        # any of ID/MD/HD x SF/VF
experiment   = "1000"           # SKU scale: "100", "1000", "10000", "100000"
weight_mode  = :uniform         # :uniform, :frequency, or :random
include("main_benchmark_settings.jl")
```
Results are saved to the `results` directory as `<experiment>_benchmark_<dependency>[_<weight_mode>].csv`.

To run all six dataset configurations in parallel via tmux (one session per SKU scale, edit the `experiments` array in the script to select scales):
```bash
bash run_benchmarks.sh              # uniform storage requirements
bash run_benchmarks.sh frequency    # frequency-proportional
bash run_benchmarks.sh random       # random
```

### Reproducing the Results of the Article
1. Run the full grid (all six dataset configurations at all four SKU scales) for each of the three weight modes as described above. `results_chigate/` contains the CHI results of the same grid with the aggregate independence gate enabled, which is the production configuration reported in the article.
2. Aggregate all runs into a single file:
   ```bash
   julia aggregate_benchmarks.jl
   ```
   This writes `results/overall_results.csv` (CHI rows taken from the gated run) and `results/overall_results_chiungated.csv` (ungated CHI rows, used for the gate ablation).
3. Generate the article's tables from the aggregated results:
   ```bash
   julia evaluation/generate_detailed_table.jl        # computation times
   julia evaluation/generate_split_ratio_table.jl     # split ratio by scale and order density
   julia evaluation/generate_structure_table.jl       # split ratio by dataset structure
   julia evaluation/generate_warehouse_table_complete.jl
   julia evaluation/generate_weighted_table.jl        # heterogeneous storage requirements
   julia evaluation/generate_gate_ablation_table.jl
   ```
Instance generation is seeded deterministically from the scenario parameters, so re-runs produce identical instances. The `only_pairs` keyword of the `BENCHMARK` function additionally allows re-running a subset of scenarios in isolation while preserving these seeds.

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
* Lin, Y.-H., Zhu, S., Wang, R. (2025). Multi-Warehouse Assortment Selection: Minimizing Order Splitting in E-Commerce Logistics. Production and Operations Management. https://doi.org/10.1177/10591478251365581
* Zhu, S., Hu, X., Huang, K. et al. (2021). Optimization of product category allocation in multiple warehouses to minimize splitting of online supermarket customer orders. European Journal of Operational Research. https://doi.org/10.1016/j.ejor.2020.08.024
* Catalan, A., Fisher, M. (2012). Assortment Allocation to Distribution Centers to Minimize Split Customer Orders. SSRN. https://doi.org/10.2139/ssrn.2166687
