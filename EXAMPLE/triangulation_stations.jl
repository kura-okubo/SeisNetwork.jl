include("/Users/kurama/.julia/dev/SeisNetwork/src/SeisNetwork.jl")
using .SeisNetwork
#---Input parameters---#
finame = "./src/input_stationlist.jld2"
fogmtname = "./output/station_triangluration.jld2"
fodvvname = "./output/station_dvvtriangluration.jld2"
maxdistance = 60e3 #[m] thresholding triangle mesh with this distance
mindistance = 1e3 #[m] thresholding triangle mesh with this distance
maxangle = 160.0 #[deg] thresholding triangle mesh with this angle
minangle = 20.0 #[deg] thresholding triangle mesh with this angle
studyarea = [37, -119.5, 35, -122]
plotfig = true
figname = "stationloc.png"
#-----------------------#

InputDict = Dict("finame" => finame,
                 "fogmtname" => fogmtname,
                 "fodvvname" => fodvvname,
                 "maxdistance" => maxdistance,
                 "mindistance" => mindistance,
                 "maxangle" => maxangle,
                 "minangle" => minangle,
                 "studyarea" => studyarea,
                 "plotfig"  =>plotfig,
                 "figname"  =>figname
                 )
triangulation(InputDict)
