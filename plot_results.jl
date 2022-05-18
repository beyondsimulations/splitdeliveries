# Activate necessary environment
import Pkg
Pkg.activate("splitdeliveries")

# Activate necessary packages
using StatsPlots
using LaTeXStrings
using CSV
using DataFrames
using Statistics

# Specify the experiment to evaluate
experiment = "2_sensitivity"
dependency = "IND"

# Load the experiment data
parcels = CSV.read("results/$(experiment)_a_parcels_sent_$dependency.csv", DataFrame)
time    = CSV.read("results/$(experiment)_b_duration_$dependency.csv", DataFrame)
split   = CSV.read("results/$(experiment)_e_split_reduction_$dependency.csv", DataFrame)

# Define Plot Style
plot_font = "Computer Modern";
default(fontfamily=plot_font, size = (600,400), framestyle=:box, label=nothing, grid=false, tickfontsize=10, legend = :topleft)

# Plot A: Parcels sent
plt = time
plt_a = plt[plt[!,:buffer].==0,:]
plt_b = plt[plt[!,:buffer].> 0,:]

plt_a = groupby(plt_a,:capacity)
plt_a = combine(plt_a, Symbol.(names(plt_a)[4:end-1]) .=> mean .=> Symbol.(names(plt_a)[4:end-1]))
plt_a = sort(plt_a, :capacity)
display(@df plt_a plot(:capacity, cols(2:ncol(plt_a))))

plt_b = groupby(plt_b,:capacity)
plt_b = combine(plt_b, Symbol.(names(plt_b)[4:end-1]) .=> mean .=> Symbol.(names(plt_b)[4:end-1]))
plt_b = sort(plt_b, :capacity)
display(@df plt_b plot(:capacity, cols(2:ncol(plt_b))))







plot(t -> p(t,1,0,0.4),0,1500,  labels= L"b = 0.4, p(t) = t^b", size = (600,400), linestyle = :solid, linewidth = 1.0, legend = :bottomleft, 
palette = :Dark2_5, xlabel= L"selling time $t$", ylabel=L"profit contribution $p(t)$")
plot!(t -> p(t,1,0.01,0.4),0,1500,  labels= L"b = 0.4, p(t) = t^b-0.01t",linestyle = :solid, linewidth = 1.0)
plot!(t -> p(t,1,0,0.3),0,1500, labels= L"b = 0.3, p(t) = t^b", linestyle = :solid, linewidth = 1.0)
plot!(t -> p(t,1,0.01,0.3),0,1500,  labels= L"b = 0.3, p(t) = t^b-0.01t",linestyle = :solid, linewidth = 1.0)

savefig("profit_full.pdf")

pd(t,pr,co,el) = pr * el * t^(el-1) - co
plot(t -> pd(t,1,0,0.4),0,1500,  labels= L"b = 0.4, p(t) = t^b",linestyle = :solid, linewidth = 1.0, legend = :topright, 
palette = :Dark2_5, xlabel= L"selling time $t$", ylabel=L"profit contribution $p(t)$")
plot!(t -> pd(t,1,0.01,0.4),0,1500,  labels= L"b = 0.4, p(t) = t^b-0.01t",linestyle = :solid, linewidth = 1.0)
plot!(t -> pd(t,1,0,0.3),0,1500, size = (450,400), labels= L"b = 0.3, p(t) = t^b", linestyle = :solid, linewidth = 1.0)
plot!(t -> pd(t,1,0.01,0.3),0,1500,  labels= L"b = 0.3, p(t) = t^b-0.01t",linestyle = :solid, linewidth = 1.0)

savefig("profit_der.pdf")