# Split-Delivery Minimisation in E-Commerce Environments
## About the project
This repository contains the code used in a currently anonymous article. It contains two new heuristics, the QMK heuristic and the CHI heuristic, as well as four competing heuristics that represent the state-of-the-art in split delivery minimisation. With this repository all heuristics and optimisations applied in the article can be reproduced. Furthermore, the code can be applied to any other split-delivery minimisation problem.

## Built with
* [Julia](https://github.com/JuliaLang)
* [GAMS](https://www.gams.com)

## Heuristics
This repository contains the following heuristics:
* QMKOPT: optimisation model from the QMK heuristic solved with the solver CPLEX
* QMK: optimisation model from the QMK heuristic solved with the solver BONMIN
* CHISOL: CHI heuristic without local search
* CHI: CHI heuristic with local search
* KLINK: KLINK heuristic by [Zhang, W.-H. Lin, M. Huang and X. Hu (2021)](https://doi.org/10.1016/j.ejor.2019.07.004), replicated based on the article
* GP: greedy pairs heuristic by [A. Catalan and M. Fisher (2012)](https://doi.org/10.2139/ssrn.2166687), replicated based on the article
* GS: greedy seeds heuristic by [A. Catalan and M. Fisher (2012)](https://doi.org/10.2139/ssrn.2166687), replicated based on the article
* BS: bestselling heuristic by [A. Catalan and M. Fisher (2012)](https://doi.org/10.2139/ssrn.2166687), replicated based on the article
* OPT: optimisation model to solve the split-delivery minimisation with CPLEX based on the models of [A. Catalan and M. Fisher (2012)](https://doi.org/10.2139/ssrn.2166687) and [Zhang, W.-H. Lin, M. Huang and X. Hu (2021)](https://doi.org/10.1016/j.ejor.2019.07.004)
* RND: random allocation of SKUs to warehouses
The choice of heuristics as well as parameters can be controlled in the file "main_benchmark_seetings.jl".

## Solvers
* [CPLEX](https://www.ibm.com/analytics/cplex-optimizer)
* [BONMIN](https://github.com/coin-or/Bonmin)

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

### Apply the benchmark on own datasets
To apply the benchmark on your own dataset you have to allocate two new files:
1. capacity_"your_experiment_name".csv in the folder "capacity"
2. transactions_"your_experiment_name".csv in the folder "transactions"
The capacity file has to contain a matrix. Each new column indicates a new warehouse with the number representing the available space. Each new row indicates a new capacity constellation that will be tested during the benchmark. For example [5 5 0; 6 4 2] represents two capacity constellations that will be benchmarked. The first constellation consists of two warehouses with a capacity of 5 each while the second represents 3 warehouses with a capacity of 6, 4 and 2. For further details take a look at the uploaded files used in our article in the folder "capacity".
The transaction file has to contain a matrix as well. Each column represents a new SKU while each new row is an order. The dataset has to be binary indicating for each order whether the corresponding SKU in the column was part of the order. For further details take a look at the uploaded files used in our article in the folder "transactions".
Execute the follwing steps two start the benchmark after both files are allocated correctly:
1. Rename the variable "experiment" in "main_benchmark_settings.jl" after "your_experiment_name" used in both files
2. Further parameters can be adjusted as explained in "main_benchmark_settings.jl"
3. Execute the file "main_benchmark_settings.jl" to start the benchmark
4. The results will be saved in the folder "results"

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

* [Zhang, W.-H. Lin, M. Huang and X. Hu (2021)](https://doi.org/10.1016/j.ejor.2019.07.004)
* [A. Catalan and M. Fisher (2012)](https://doi.org/10.2139/ssrn.2166687)
