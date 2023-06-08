# Optimising SKU-Warehouse Allocations to minimise Split Parcels in E-Commerce Environments
## About the project
This repository contains the code used in a research article. It contains two new heuristics, the QMK heuristic and the CHI heuristic, as well as four competing heuristics that represent the state-of-the-art in split delivery minimisation. With this repository all heuristics and optimisations applied in the article can be reproduced. Furthermore, the code can be applied to any other split-delivery minimisation problem.

## Built with
* [Julia](https://github.com/JuliaLang)
* [GAMS](https://www.gams.com)

## Heuristics
This repository contains the following heuristics:
* QMKOPT: optimisation model from the QMK heuristic solved with the solver CPLEX
* QMK: optimisation model from the QMK heuristic solved with the solver BONMIN
* CHISOL: CHI heuristic without local search
* CHI: CHI heuristic with local search
* KLINK: K-LINK heuristic by [S. Zhu, X. Hu and K. Huang et al. (2021)](https://doi.org/10.1016/j.ejor.2020.08.024), replicated based on the article
* GP: GREEDY PAIRS heuristic by [A. Catalan and M. Fisher (2012)](https://doi.org/10.2139/ssrn.2166687), replicated based on the article
* GS: GREEDY SEEDS heuristic by [A. Catalan and M. Fisher (2012)](https://doi.org/10.2139/ssrn.2166687), replicated based on the article
* BS: BESTSELLERS heuristic by [A. Catalan and M. Fisher (2012)](https://doi.org/10.2139/ssrn.2166687), replicated based on the article
* OPT: optimisation model to solve the split-delivery minimisation with CPLEX based on the models of [A. Catalan and M. Fisher (2012)](https://doi.org/10.2139/ssrn.2166687) and [S. Zhu, X. Hu and K. Huang et al. (2021)](https://doi.org/10.1016/j.ejor.2020.08.024)
* RND: random allocation of SKUs to warehouses.

The choice of heuristics as well as parameters can be controlled in the file "main_benchmark_seetings.jl".

## Solvers
* [CPLEX](https://www.ibm.com/analytics/cplex-optimizer)
* [SBB](https://www.gams.com/latest/docs/S_SBB.html)

## Getting started
### Prerequisites
Julia and GAMS have to be installed on the machine executing the code in this repository. Furthermore, a valid GAMS license is neccessary to execute the heuristics QMKOPT, QMK and OPT.

### Installation
1. Install Julia and GAMS (use a valid GAMS license if QMKOPT, QMK or OPT should be executed)
1. Clone the repo
2. Execute the file "install_packages.jl" to install all neccessary packages (also listed under "Associated Repositories")

### Reproduce benchmarks from the article
1. Adjust the parameters in the file "main_benchmark_settings.jl" (more details within the comments of the file)
2. Execute the file "main_benchmark_settings.jl" to start the benchmark
3. The results will be saved in the folder "results"

## License
Distributed under the MIT License. See `LICENSE.txt` for more information.

## Associated Repositories
* [CSV.jl](https://github.com/JuliaData/CSV.jl)
* [DataFrames.jl](https://github.com/JuliaData/DataFrames.jl)
* [Combinatorics.jl](https://github.com/JuliaMath/Combinatorics.jl)
* [Distributions.jl](https://github.com/JuliaStats/Distributions.jl)
* [GAMS.jl](https://github.com/JuliaMath/Combinatorics.jl)
* [JuMP.jl](https://github.com/jump-dev/JuMP.jl)
* [Octavian.jl](https://github.com/JuliaLinearAlgebra/Octavian.jl)
* [LoopVectorization.jl](https://github.com/JuliaSIMD/LoopVectorization.jl)
* [Combinatorics.jl](https://github.com/JuliaMath/Combinatorics.jl)

## Acknowledgments

* [S. Zhu, X. Hu and K. Huang et al. (2021)](https://doi.org/10.1016/j.ejor.2020.08.024)
* [A. Catalan and M. Fisher (2012)](https://doi.org/10.2139/ssrn.2166687)
