module SplitDeliveries

using DataFrames
using DelimitedFiles
using CSV
using Combinatorics
using Distributions
using JuMP
using Ipopt
using Juniper
using HiGHS
using Gurobi
using SCIP
using LoopVectorization
using Plots
using StatsPlots
using Dates
using Measures

using LinearAlgebra
using Random
using Statistics
using SparseArrays

const PROJDIR = dirname(dirname(@__DIR__))

include(joinpath(PROJDIR, "functions/functions_basic.jl"))
include(joinpath(PROJDIR, "functions/functions_catalan.jl"))
include(joinpath(PROJDIR, "functions/functions_k-links.jl"))
include(joinpath(PROJDIR, "functions/functions_chisquare.jl"))
include(joinpath(PROJDIR, "functions/functions_localsearch.jl"))
include(joinpath(PROJDIR, "functions/new_transactions.jl"))
include(joinpath(PROJDIR, "functions/functions_lin.jl"))

include(joinpath(PROJDIR, "heuristics/heuristic_qmkp.jl"))
include(joinpath(PROJDIR, "heuristics/optimisation_equalcap.jl"))
include(joinpath(PROJDIR, "heuristics/optimisation_buffercap.jl"))
include(joinpath(PROJDIR, "heuristics/heuristic_klinks.jl"))
include(joinpath(PROJDIR, "heuristics/heuristic_greedyseeds.jl"))
include(joinpath(PROJDIR, "heuristics/heuristic_greedypairs.jl"))
include(joinpath(PROJDIR, "heuristics/heuristic_greedyorders.jl"))
include(joinpath(PROJDIR, "heuristics/heuristic_bestselling.jl"))
include(joinpath(PROJDIR, "heuristics/heuristic_chisquare.jl"))
include(joinpath(PROJDIR, "heuristics/heuristic_mci.jl"))
include(joinpath(PROJDIR, "heuristics/heuristic_iih.jl"))

include(joinpath(PROJDIR, "functions/functions_benchmark_helpers.jl"))
include(joinpath(PROJDIR, "main_benchmark.jl"))

end
