# Install all necessary packages upon the first use of this repository
using Pkg
Pkg.activate("splitdeliveries")

# Install packages
Pkg.add("DataFrames")
Pkg.add("CSV")
Pkg.add("Combinatorics")
Pkg.add("Distributions")
Pkg.add("JuMP")
Pkg.add("HiGHS")
Pkg.add("Gurobi")
Pkg.add("AmplNLWriter")
Pkg.add("Bonmin_jll")
Pkg.add("Octavian")
Pkg.add("LoopVectorization")
Pkg.add("Plots")
Pkg.add("StatsPlots")
Pkg.add("Measures")

# Precompile all packages to reduce future startup times
Pkg.precompile()

println("Environment setup complete!")
