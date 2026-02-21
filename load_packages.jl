# activate necessary environment
import Pkg
Pkg.activate("splitdeliveries")

# load all necessary packages (repositories)
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

# load all necessary packages from the standard library module
using LinearAlgebra
using Random
using Statistics
using SparseArrays

## import all basic functions
include("functions/functions_basic.jl")
include("functions/functions_catalan.jl")
include("functions/functions_k-links.jl")
include("functions/functions_chisquare.jl")
include("functions/new_transactions.jl")
include("functions/functions_lin.jl")

## import the heuristic functions
include("heuristics/heuristic_qmkp.jl")
include("heuristics/optimisation_equalcap.jl")
include("heuristics/optimisation_buffercap.jl")
include("heuristics/heuristic_klinks.jl")
include("heuristics/heuristic_greedyseeds.jl")
include("heuristics/heuristic_greedypairs.jl")
include("heuristics/heuristic_greedyorders.jl")
include("heuristics/heuristic_bestselling.jl")
include("heuristics/heuristic_chisquare.jl")
include("heuristics/heuristic_mci.jl")
include("heuristics/heuristic_iih.jl")

## import the main comparison function
include("main_benchmark.jl")
