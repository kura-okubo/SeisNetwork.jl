module SeisNetwork

using SeisIO, Triangle, LinearAlgebra, Statistics, Geodesy, JLD2, ORCA, PlotlyJS

using LightXML, Dates, HDF5

# import pre and post processing tools
include("triangulation.jl")
include("stationkml2jld.jl")
include("plot_dvv.jl")



end # module
