export stationdensity_selection

"""
stationdensity_selection(InputDict::Dict)
return stationlist which is density-averaged over space.
"""
function stationdensity_selection(InputDict::Dict)

    starttime = InputDict["starttime"]
    endtime = InputDict["endtime"]
    kmlfile = InputDict["kmlfile"]
    studyarea = InputDict["studyarea"]
    includenet = InputDict["includenet"]
    priority = InputDict["priority"]
    TotalNumofStation = InputDict["TotalNumofStation"]
    nx = InputDict["nx"]
    ny = InputDict["ny"]
    plotfig = InputDict["plotfig"]

    # parse kml station file
    println("---Parse KML file---")
    xdoc = parse_file(kmlfile);
    xroot = root(xdoc);
    ces = collect(child_elements(xroot));

    itmp = 1
    stid = 0
    while true
        itmp
        if !isempty(ces[itmp]["name"])
            stid = itmp
            break;
        end
        itmp += 1
        if itmp > 100
            error("Station network is more than 100. abort.")
        end
    end

    s1 = ces[stid:end]

    stationnum = length(s1)

    stationlist = []
    # 1. read station info
    for i = 1:stationnum
        stationname = content(s1[i]["name"][1])
        stationloc = content(s1[i]["Point"][1]["coordinates"][1])
        coord_lon = parse.(Float64, split(stationloc, ","))[1]
        coord_lat = parse.(Float64, split(stationloc, ","))[2]

        # find operating range
        txt = string(s1[i]["Style"][1]["BalloonStyle"][1]["text"])
        tl  =split(string(txt), ";")
        operatingtime =tl[findfirst(occursin.("Operating Range", tl))]
        startoptime =  DateTime(split(operatingtime)[3][1:10])
        endtopime =  DateTime(split(operatingtime)[5][1:10])

        datacentertxt =tl[findfirst(occursin.("Data Center", tl))]

        if startoptime <= starttime && endtime <= endtopime
            # take only stations which is operated over study time.
            st = Dict("stationname" => stationname,
                    "coord_lon" => coord_lon,
                    "coord_lat" => coord_lat)
            push!(stationlist, st)
        end
    end

    # 2. plot location
    if plotfig

        fodir = "./fig"
        if !ispath(fodir) mkpath(fodir); end

        px = [stationlist[x]["coord_lon"] for x in 1:length(stationlist)]
        py = [stationlist[x]["coord_lat"] for x in 1:length(stationlist)]
        text = [stationlist[x]["stationname"] for x in 1:length(stationlist)]
        p1 = scatter(x=px, y=py, mode="markers", text=text,
                       textposition="top center")

        layout = Layout(;width=1200, height=800, title="All station network",
                        xaxis_range=[studyarea[4], studyarea[2]],
                       yaxis_range=[studyarea[3], studyarea[1]],
                       xaxis_title="lon",
                       yaxis_title="lat",
                       showlegend=false)

         p = PlotlyJS.plot(p1, layout)
         PlotlyJS.savefig(p, fodir*"/"*InputDict["figname"])

    end

    # 3. compute station density
    println("---Compute station density---")

    dx = (studyarea[2] - studyarea[4])/nx
    dy = (studyarea[1] - studyarea[3])/ny
    NG = nx * ny
    grid = []
    for i=1:nx
        for j=1:ny
            xl = (i-1) * dx + studyarea[4]
            xr = (i) * dx + studyarea[4]
            yb = (j-1) * dy + studyarea[3]
            yt = (j) * dy + studyarea[3]
            xm = (xl+xr)/2
            ym = (yb+yt)/2
            gridtmp = Dict("xm"=>xm, "ym"=>ym, "xl"=>xl, "xr"=>xr, "yb"=>yb, "yt"=>yt)
            push!(grid, gridtmp)
        end
    end

    stationidlist = []
    for gid = 1:NG
        xl = grid[gid]["xl"]
        xr = grid[gid]["xr"]
        yb = grid[gid]["yb"]
        yt = grid[gid]["yt"]
        numstationpergrid=0
        stationid=Int[]
        for i = 1:length(stationlist)
            x1 = stationlist[i]["coord_lon"]
            y1 = stationlist[i]["coord_lat"]

            if insquare(x1, y1, xl, xr, yb, yt, edgeinclude=true)
                numstationpergrid += 1
                push!(stationid, i)
            end
        end
        grid[gid]["numstationpergrid"] = numstationpergrid
        push!(stationidlist, stationid)
    end

    # 4. reduce the station density
    #   - sort from higher num of station
    #   - remove one of pairs with closest distance and priority
    println("---Averaging station density---")

    r_grid= deepcopy(grid)
    r_stationlist = deepcopy(stationlist)
    r_stationidlist = deepcopy(stationidlist)

    removestations!(TotalNumofStation, r_grid, r_stationlist, r_stationidlist, includenet, priority)

    r_stationlist = [stationlist[x] for x in collect(Iterators.flatten(r_stationidlist))]

    if plotfig

        fodir = "./fig"

        px = [r_stationlist[x]["coord_lon"] for x in 1:length(r_stationlist)]
        py = [r_stationlist[x]["coord_lat"] for x in 1:length(r_stationlist)]
        text = [r_stationlist[x]["stationname"] for x in 1:length(r_stationlist)]
        p1 = scatter(x=px, y=py, mode="markers", text=text,
                       textposition="top center")


       layout = Layout(;width=1200, height=800, title="Density averaged station network",
                       xaxis_range=[studyarea[4], studyarea[2]],
                      yaxis_range=[studyarea[3], studyarea[1]],
                      xaxis_title="lon",
                      yaxis_title="lat",
                      showlegend=false)

        p = PlotlyJS.plot(p1, layout)
        PlotlyJS.savefig(p, fodir*"/"*"densityaveraged_"*InputDict["figname"])

    end

    #===
    save for SeisDownload input file
    ===#

    fofiledir = "./files"
    if !ispath(fofiledir) mkpath(fofiledir); end
    stationlistfoname = fofiledir*"/input_stationlist.jld2"
    jldopen(stationlistfoname, "w") do file
        file["stationlist"] = stationlist;
    end

    #===
    # save stationlist for gmt
    ===#

    # gather network name
    networks = Dict()
    println("---Station networks---")
    for i = 1:length(stationlist)
        networkname =  split(stationlist[i]["stationname"], ".")[1]
        if !haskey(networks, networkname)
            println(networkname)
            #add this network to dict
            networks["$(networkname)"] = Array[]
        end
    end


    for i = 1:length(stationlist)
        stationname = stationlist[i]["stationname"]
        coord_lon = stationlist[i]["coord_lon"]
        coord_lat = stationlist[i]["coord_lat"]

        # Here, elevation is not included because no elevetion in IRIS kml file.
        # However, you can find elevation when analysing dv/v as CorrData has the elevation of stations.

        for key = keys(networks)
            if occursin(stationname[1:3], "$(key).", )
                push!(networks["$(key)"], [coord_lon, coord_lat])
            end
        end
    end

    #save as jld2
    jldopen(fofiledir*"/gmt_stationcoords.jld2", "w") do file
        for key = keys(networks)
            file["$(key)"] = networks["$(key)"]
        end
    end

    return nothing

end

"""
insquare(x1, y1, xl, xr, yb, yt; edgeinclude::Bool=true)
return if x1, y1 is in square
"""
function insquare(x1, y1, xl, xr, yb, yt; edgeinclude::Bool=true)
    if edgeinclude
        if x1 >= xl && x1 <= xr && y1 >= yb && y1 <= yt
            return true
        else
            return false
        end
    else
        if x1 > xl && x1 < xr && y1 > yb && y1 < yt
            return true
        else
            return false
        end
    end
end


"""
removestations!(TotalNumofStation::Int, r_grid::Array{Any,1}, r_stationlist::Array{Any,1},
     r_stationidlist::Array{Any,1}, includenet::Array{String,1}, priority::Array{String,1})
remove stations with averaging density in grid
"""
function removestations!(TotalNumofStation::Int, r_grid::Array{Any,1}, r_stationlist::Array{Any,1},
     r_stationidlist::Array{Any,1}, includenet::Array{String,1}, priority::Array{String,1})

     NG = length(r_grid)

     r_numstlistsort = [r_grid[x]["numstationpergrid"] for x in 1:NG]
     r_sortid = sortperm(r_numstlistsort, rev=true)

     icount = 1
     while true

         if sum(r_numstlistsort) <= TotalNumofStation
             break
         elseif sum(r_numstlistsort) == 0
             error("All station is removed from list. abort.")
         elseif icount > 10000
             error("iteration number exceeds 10000. too many station removal. abort.")
         end

         r_gid = first(r_sortid) # remove from this grid
         st = r_stationidlist[r_gid]
         # compute distance with eachother
         rdist = 1e12
         r_ii = -1
         r_jj = -1
         r_netii = ""
         r_netjj = ""


         for i = 1:length(st)-1
             for j = i+1:length(st)

                 ii = st[i]
                 jj = st[j]

                 id1 = r_stationlist[ii]["stationname"]
                 id2 = r_stationlist[jj]["stationname"]
                 net1 =  string(split(id1, ".")[1])
                 net2 =  string(split(id2, ".")[1])
                 x1 = r_stationlist[ii]["coord_lon"]
                 y1 = r_stationlist[ii]["coord_lat"]
                 x2 = r_stationlist[jj]["coord_lon"]
                 y2 = r_stationlist[jj]["coord_lat"]
                 trial_dist = norm([(x2-x1), (y2-y1)])

                 # decide if we remove this
                 if  @isdefined includenet
                     if  any(occursin.(includenet, net1)) || any(occursin.(includenet, net2))
                         continue;
                     end
                 end

                 if trial_dist < rdist
                     #update potential removal station
                     r_ii = ii
                     r_jj = jj
                     r_netii = net1
                     r_netjj = net2
                     rdist = trial_dist
                 end
             end
         end

         if r_ii == -1 || r_jj == -1
             error("no station is chosen. abort.")
         end

         #evaluate removing either r_ii or r_jj
         if  !@isdefined priority
             priority = []
         end

         ri1 = findfirst(x -> x==r_netii, priority)
         ri2 = findfirst(x -> x==r_netjj, priority)

         if isnothing(ri1) || isnothing(ri2)
             removedid = r_ii #no priority
         elseif ri1 >= ri2
             removedid = r_ii
         else
             removedid = r_jj
         end

         # pop out removedid from list
         r_grid[r_gid]["numstationpergrid"] = r_numstlistsort[r_gid] - 1
         filter!(e -> e != removedid, r_stationidlist[r_gid])

         r_numstlistsort = [r_grid[x]["numstationpergrid"] for x in 1:NG]
         r_sortid = sortperm(r_numstlistsort, rev=true)

         icount += 1
     end
end
