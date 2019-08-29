export plot_dvv

function plot_dvv(InputDict::Dict)

    basefiname      = InputDict["basefiname"]
    triangulationfiname = InputDict["triangulationfiname"]
    outputformat = InputDict["outputformat"]
    plotfig = InputDict["plotfig"]
    figname = InputDict["figname"]

    #1 load data
    t = jldopen(triangulationfiname)
    tris = t["tri"]

    # load existin dvv histories
    basefidir = basefiname[1:end-length(basename(basefiname))]
    dvvpath = ls(basefidir)
    dvvstalist = Array{String, 1}(undef,0)

    for path in dvvpath
        bn = split(basename(path)[length(basename(basefiname))+1:end], ".")
        try
            s1 = join(bn[1:2], ".")
            s2 = join(bn[5:6], ".")
            dvvstalist = vcat(dvvstalist, (s1, s2))
        catch
            #println(path)
            passid = findfirst(x-> x==path, dvvpath)
            dvvpath = dvvpath[setdiff(1:end, passid), :]
        end
    end

    # loop triangles

    loc_and_dvv = []
    # synchronize time to take mean

    stday = Inf
    etday = 0

    for tri in tris
        centerlonlat = tri["centerxy"]
        edgeinfo = tri["edgeinfo"]

        numofedge = 0
        dvvsum = []
        Tsum = []

        for edge in edgeinfo
            stn1 = edge[1]
            stn2 = edge[2]

            #println((stn1, stn2))
            # search corresponding dvv history
            dvvid = findfirst(x->x == (stn1, stn2), dvvstalist)
            dvvidrev = findfirst(x->x == (stn2, stn1), dvvstalist)

            if !isnothing(dvvid)
                fipath = dvvpath[dvvid]
            elseif !isnothing(dvvidrev) # reverse
                fipath = dvvpath[dvvidrev]
            else
                #println("debug: no data on this edge.")
            end

            try
                fiedge = jldopen(fipath, "r")

                push!(Tsum, fiedge["T"])
                push!(dvvsum, fiedge["dvv"])

                close(fiedge)
            catch
                #println("debug: load error edge loading.")
            end
        end

        avgdvv = []
        avgT = []

        for tsumid = 1:length(Tsum)
            for ii = 1:length(Tsum[tsumid])
                # try to average at this time
                testtavg = Tsum[tsumid][ii]
                dvvtemp = []

                if isnothing(findfirst(x -> x==testtavg, avgT))
                    # this time step is already averaged
                    overlapid = zeros(length(Tsum))
                    edgecount = 0
                    for jj = 1:length(Tsum)
                        oid = findfirst(x -> x==testtavg, Tsum[jj])
                        if !isnothing(oid)
                            push!(dvvtemp, dvvsum[tsumid][oid])
                            edgecount += 1
                        end
                    end

                    dvvtemp = filter(!isnan, Array{Float64, 1}(dvvtemp))
                    if !isempty(dvvtemp) && edgecount > 0
                        #no data at this time
                        push!(avgT, testtavg)
                        push!(avgdvv, mean(dvvtemp))
                    else
                        push!(avgT, testtavg)
                        push!(avgdvv, 0)
                    end
                end
            end
        end

        sortid = sortperm(avgT)
        avgT = avgT[sortid]
        avgdvv = avgdvv[sortid]
        avgT = convert(Array{Int64,1}, avgT)
        avgvv =  convert(Array{Float64,1}, avgdvv)
        # println(sortid)


        # smoothing dv/v
        if InputDict["smoothing"]
            trace1 = PlotlyJS.scatter(x=timestamp.(avgT), y=avgdvv)
            SeisNoise.smooth!(avgdvv; half_win=5)
            println(avgdvv)
            trace2 = PlotlyJS.scatter(x=timestamp.(avgT), y=avgdvv, label="after smooth")
            p = PlotlyJS.plot([trace1, trace2])
            display(p)
            readline()
        end

        if !isempty(avgT)
            stday_test = avgT[1]
            etday_test = avgT[end]
            if stday_test < stday stday = stday_test; end
            if etday_test > etday etday = etday_test; end
        end


        push!(loc_and_dvv, (centerlonlat, edgeinfo, avgT, avgdvv))

    end

    dvvdata = Dict( "startday" => timestamp(stday),
                    "endday" => timestamp(etday),
                    "loc_and_dvv" => loc_and_dvv)

    # save data
    tstamp = timestamp(stday)

    if outputformat == "hdf5"

        tstamp =  timestamp(stday)
        foname = basefidir*"/../spatialdvv_$(tstamp).h5"

        if ispath(foname) rm(foname); end

        h5open(foname, "w") do file

            gorigin = "spatiotemporal_dvv" #original group name

            write(file, gorigin*"/startday",  dvvdata["startday"])
            write(file, gorigin*"/endday", dvvdata["endday"])

            loc_and_dvv = dvvdata["loc_and_dvv"]

            #create group on each element
            for i = 1:length(loc_and_dvv)
                centerlonlat  = loc_and_dvv[i][1]
                edgeinfo      = loc_and_dvv[i][2]
                T             = loc_and_dvv[i][3]
                dvvsum        = loc_and_dvv[i][4]

                groupname = join([gorigin, "ele$i"], "/")

                write(file, groupname*"/dvv", Array{Float64, 1}(dvvsum))
                write(file, groupname*"/centerlon", centerlonlat[1])
                write(file, groupname*"/centerlat", centerlonlat[2])
                write(file, groupname*"/T", Array{Int64, 1}(T))

                for edgeid = 1:length(edgeinfo)
                    edge = edgeinfo[edgeid]
                    sta1 = edge[1]
                    sta2 = edge[2]

                    edgegroupname = join([groupname, "edges$(edgeid)"], "/")

                    write(file, edgegroupname*"/sta1", sta1)
                    write(file, edgegroupname*"/sta1edge", edge[3])
                    write(file, edgegroupname*"/sta2", sta2)
                    write(file, edgegroupname*"/sta2edge", edge[4])

                end
            end
        end


    elseif outputformat == "jld2"

        foname = basefidir*"/../spatialdvv_$(tstamp).jld2"
        if ispath(foname) rm(foname); end

        outfile = jldopen(foname, "w")
        outfile["dvvdata"] = dvvdata
        close(outfile)

    else
        println("output format not exist. No output has been done.")
    end

    if plotfig
        fodir = "./fig"
        if !ispath(fodir) mkpath(fodir); end

        println("Plot figure is is working in progress. Please use Matlab or Python to plot.")
        #
        # mul = 400
        # height= abs(diff([studyarea[1], studyarea[3]])[1]) * mul
        # width = abs(diff([studyarea[2], studyarea[4]])[1]) * mul
        # layout = Layout(;width=width, height=height, title="Triangulation of station network",
        #                xaxis_range=[studyarea[4], studyarea[2]],
        #                yaxis_range=[studyarea[3], studyarea[1]],
        #                xaxis_title="lon",
        #                yaxis_title="lat",
        #                showlegend=false)
        #
        # p = PlotlyJS.plot(traceall, layout)
        # PlotlyJS.savefig(p, fodir*"/"*InputDict["figname"])
        #display(p)
        #println("type return...")
        #readline()
    end

    println("plot dvv has been successfully done.")
    return nothing

end
