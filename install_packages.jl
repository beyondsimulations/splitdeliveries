# Install all necessary packages upon the first use of this repository
using Pkg
Pkg.activate("splitdeliveries")
Pkg.add("DataFrames")
Pkg.add("CSV")
Pkg.add("Combinatorics")
Pkg.add("Distributions")
Pkg.add("JuMP")
Pkg.add("GAMS")
Pkg.add("Octavian")
Pkg.add("LoopVectorization")
Pkg.add("Plots")
Pkg.add("StatsPlots")
Pkg.add("Measures")
