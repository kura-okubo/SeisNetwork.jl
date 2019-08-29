export triangulation

function get_angle_from_sides(a::Float64, b::Float64, c::Float64)
    return acosd((a ^ 2 + b ^ 2 - c ^ 2) / (2 * a * b))
end

function checkangle(a::Float64, b::Float64, c::Float64)
    angle = []
    try
        push!(angle, get_angle_from_sides(a, b, c))
        push!(angle, get_angle_from_sides(c, a, b))
        push!(angle, get_angle_from_sides(b, c, a))
        return sort(angle)
    catch
        return [Nan, Nan, Nan]
    end
end

function triangulation(InputDict::Dict)

    finame = InputDict["finame"]
    fogmtname = InputDict["fogmtname"]
    fodvvname = InputDict["fodvvname"]
    maxdistance = InputDict["maxdistance"]
    mindistance = InputDict["mindistance"]
    maxangle = InputDict["maxangle"]
    minangle = InputDict["minangle"]
    studyarea = InputDict["studyarea"]
    plotfig = InputDict["plotfig"]

    # for convert from LLA coordinate to ENU coordinate
    origin_lla = LLA((studyarea[1]+studyarea[3])/2, (studyarea[2]+studyarea[4])/2, 0.0)
    trans = ENUfromLLA(origin_lla, wgs84)

    t = jldopen(finame)["stationlist"]

    dist = []
    stationnum=length(t)

    # compute Delaunay triangulation of stations
    sta_lonlat = Array{Float64, 2}(undef, stationnum, 2)
    for i = 1:stationnum
        sta_lonlat[i, :] = Float64[t[i]["coord_lon"] t[i]["coord_lat"]]
    end

    println("compute Delaunay.")
    tri_temp = Triangle.basic_triangulation_vertices(sta_lonlat)

    traceall = GenericTrace{Dict{Symbol,Any}}[]

    tri = []
    dist = []

    for i = 1:length(tri_temp)
        # compute distances of all pairs
        p1 = [tri_temp[i][1 ,1], tri_temp[i][1,2]] #[lon, lat]
        p2 = [tri_temp[i][2 ,1], tri_temp[i][2,2]]
        p3 = [tri_temp[i][3 ,1], tri_temp[i][3,2]]

        pp1 = LLA(p1[2], p1[1]) #LLA(lat, lon, alt = 0.0)
        pp2 = LLA(p2[2], p2[1])
        pp3 = LLA(p3[2], p3[1])

        # pe1 = trans(pp1)
        # pe2 = trans(pp2)
        # pe3 = trans(pp3)
        #
        # v1 = diff(hcat([pe1[1], pe2[1]], [pe1[2], pe2[2]]), dims=1)
        # v2 = diff(hcat([pe2[1], pe3[1]], [pe2[2], pe3[2]]), dims=1)
        # v3 = diff(hcat([pe3[1], pe1[1]], [pe3[2], pe1[2]]), dims=1)
        #
        # dist1 = norm(v1)
        # dist2 = norm(v2)
        # dist3 = norm(v3)

        #check
        dist1= distance(pp1, pp2);
        dist2= distance(pp2, pp3);
        dist3= distance(pp3, pp1);
        #
        # if dist11 == dist1 && dist22 == dist2 && dist33 == dist3
        #     println("Pass Distance test")
        # else
        #     println([dist11, dist1, dist22, dist2, dist33, dist3])
        #     @warn "distance wrong"
        # end

        angle = checkangle(dist1, dist2, dist3)

        # println([dist1, dist2, dist3])
        # println(angle)
        # println( any(x -> x>maxdistance, [dist1, dist2, dist3]))
        # println( any(x -> x<mindistance, [dist1, dist2, dist3]))
        # println( any(x -> x<minangle, angle))
        # println( any(x -> x>maxangle, angle))

        if any(x -> x>maxdistance, [dist1, dist2, dist3]) ||
            any(x -> x<mindistance, [dist1, dist2, dist3]) ||
            any(x -> x<minangle, angle) ||
            any(x -> x>maxangle, angle)
            continue
        else

            push!(tri, tri_temp[i])
            push!(dist, [dist1, dist2, dist3])

            trace1 = scatter(;x=tri_temp[i][[1,2],1], y=tri_temp[i][[1,2],2], mode="lines+markers")
            trace2 = scatter(;x=tri_temp[i][[2,3],1], y=tri_temp[i][[2,3],2], mode="lines+markers")
            trace3 = scatter(;x=tri_temp[i][[3,1],1], y=tri_temp[i][[3,1],2], mode="lines+markers")

            push!(traceall, trace1)
            push!(traceall, trace2)
            push!(traceall, trace3)
        end
    end

    if plotfig
        fodir = "./fig"
        if !ispath(fodir) mkpath(fodir); end

        mul = 400
        height= abs(diff([studyarea[1], studyarea[3]])[1]) * mul
        width = abs(diff([studyarea[2], studyarea[4]])[1]) * mul
        layout = Layout(;width=width, height=height, title="Triangulation of station network",
                       xaxis_range=[studyarea[4], studyarea[2]],
                       yaxis_range=[studyarea[3], studyarea[1]],
                       xaxis_title="lon",
                       yaxis_title="lat",
                       showlegend=false)

        p = PlotlyJS.plot(traceall, layout)
        PlotlyJS.savefig(p, fodir*"/"*InputDict["figname"])
        #display(p)
        #println("type return...")
        #readline()
    end

    # compute max, min, average distance
    maxedge= maximum(maximum(dist))
    minedge= minimum(minimum(dist))
    meanedge = mean(mean(dist))
    # paring with station id
    triall = []
    for i = 1:length(tri)
        edgeinfo = []
        for (j, k) in [(1,2), (2,3), (3,1)]
            edge = [tri[i][j, 1] tri[i][j, 2]; tri[i][k, 1] tri[i][k, 2]]
            sta1 = findfirst(x -> x == edge[1], sta_lonlat)
            sta2 = findfirst(x -> x == edge[2], sta_lonlat)
            sta1id = t[sta1]["stationname"]
            sta2id = t[sta2]["stationname"]
            staloc1 =  edge[1, 1:2]
            staloc2 =  edge[2, 1:2]
            midpoint = mean(edge, dims=1)
            push!(edgeinfo, (sta1id, sta2id, staloc1, staloc2, midpoint))
        end
        centerxy= mean(tri[i], dims=1)
        push!(triall, Dict("edgeinfo" => edgeinfo, "centerxy" => centerxy))
    end

    #---output for dv/v---#
    jldopen(fodvvname, "w") do file
        file["tri"] = triall
        file["maxedgesize"] = maxedge
        file["minedgesize"] = minedge
        file["meanedgesize"] = meanedge
    end

    #---output for gmt---#
    jldopen(fogmtname, "w") do file
        edges = []
        for i = 1:length(tri)
            for (j, k) in [(1,2), (2,3), (3,1)]
                edge = [tri[i][j, 1] tri[i][j, 2]; tri[i][k, 1] tri[i][k, 2]]
                if !any(x -> x == edge, edges)
                    #to avoid overlap of line
                    push!(edges, edge)
                end
            end
        end
        file["edges"] = edges
    end

    println("analysis has been successfully done.")
    return nothing

end
