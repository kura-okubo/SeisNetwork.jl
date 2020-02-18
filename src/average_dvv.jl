export average_dvv

function average_dvv(InputDict::Dict)


	starttime	= InputDict["starttime"]
	endtime     = InputDict["endtime"]
	timeunit    = InputDict["timeunit"]
	basefidir	= InputDict["basefidir"]
	IsPlotFigure 	= InputDict["IsPlotFigure"]
	freqband 	= InputDict["freqband"]
	Nfreqband 	= length(freqband) - 1

	startyear 	= Dates.year(starttime)
	endyear 	= Dates.year(endtime)

	# generate time window for computing averaged dvv
	datetime_unittime = Dates.Second(timeunit)

	# initialize starttime
	twin_left  = starttime
	twin_right = twin_left + datetime_unittime

	icount = 0
	maxcount = 1e5

	timewindows = []
	tvec 		= []
	while true
		# check if the window boundary is within start and end datetime.
		if twin_right > endtime
			# time window is outside of computed cc time
			break;
		end

		push!(timewindows, (twin_left, twin_right))
		push!(tvec, (datetime2unix(twin_left) + datetime2unix(twin_right))/2)

		# slide time window with overlap
		twin_left  = twin_right
		twin_right = twin_left + datetime_unittime

		icount += 1
		if icount >= maxcount
			error("manipulate loop reaches limit $(maxcount). Please check overlapperstack, unitnumperstack, etc.")
		end
	end
	InputDict["timewindows"] = timewindows
	# println(timewindows)
	# println(tvec)

	dvvfiname = []
	for iyear = startyear:endyear
		basedvvdir = basefidir*"/$(iyear)/dvv"
		# read dvv filepaths
		try
			append!(dvvfiname, basedvvdir*"/".*readdir(basedvvdir))
		catch y
			println(y)
			continue
		end
	end

	InputDict["dvvfiname"] = dvvfiname

	# load and store data to 3D array (stationpair, time, ifreq)
	println("---Read and store dvv data---")

	t1 = @elapsed dvv_dict_all = pmap(x -> readdvvfile(x), dvvfiname)

	@show t1

	# compute averaging dvv
	println("---Compute averaged dvv---")
	A = pmap(x->map_average_dvv(x, InputDict, dvv_dict_all), timewindows)


	dvv_mean_all  = Array{Float32, 2}(undef, 0, Nfreqband)
	dvv_std_all   = Array{Float32, 2}(undef, 0, Nfreqband)
	dvv_count_all = Array{Float32, 2}(undef, 0, Nfreqband)
	dvv_cc_all 	  = Array{Any, 1}(undef, length(timewindows))

	for i = 1:length(timewindows)
		dvv_mean_all 	= cat(dvv_mean_all, reshape(A[i][1], 1, Nfreqband), dims=1)
		dvv_std_all 	= cat(dvv_std_all, reshape(A[i][2], 1, Nfreqband), dims=1)
		dvv_count_all 	= cat(dvv_count_all, reshape(A[i][3], 1, Nfreqband), dims=1)
		dvv_cc_all[i]	= A[i][4]
	end

	# @show dvv_mean_all
	# @show dvv_std_all
	# @show dvv_count_all
	println("---Save and plot dvv---")

	if IsPlotFigure; Plots.pyplot(); end

	figdir = InputDict["figdir"]
	if ispath(figdir); rm(figdir, recursive=true) end
	mkpath(figdir);

	for ifreq = 1:Nfreqband
		freqmin = freqband[ifreq]
		freqmax = freqband[ifreq+1]

		# store tvec and dvv
		strfreq = @sprintf("%4.2f-%4.2f", round(freqmin, digits=2), round(freqmax, digits=2))

		jldopen(figdir*"/dvv_averaged_$(strfreq)Hz.jld2", "a+") do fo
			fo["timewindow"] = timewindows
			fo["tvec"] = tvec
			fo["dvv_mean_all"] = dvv_mean_all
			fo["dvv_std_all"] = dvv_std_all
			fo["dvv_count_all"] = dvv_count_all
		end

		if IsPlotFigure
			# plot dvv curve
			dvv_mean_ifreq 	= dvv_mean_all[:, ifreq]
			dvv_std_ifreq	= dvv_std_all[:, ifreq]
			dvv_count_ifreq = dvv_count_all[:, ifreq]

			tvec_stamp = timestamp.(tvec)
			xticklabels = map(x -> x[1:10], tvec_stamp)

			p1 = Plots.plot(xticklabels, dvv_mean_ifreq,
							linewidth = 1.0,
							linecolor = :black,
							ribbon=dvv_std_ifreq,
							fillcolor=:darkblue,
							ylabel="dv/v [%]",
							title=strfreq*"Hz",
							label="",
							fillalpha=0.2,
							# xlims=xlimit,
							xrotation = 45)

			# plot count histogram
			p2 = Plots.bar(xticklabels,  dvv_count_ifreq,
			 				fillcolor=:gray,
							 xrotation = 45,
							 label="",
							 ylabel="Count of cc stretch > threshold")

			p12 = Plots.plot(p1, p2, layout = (2,1), size=(1200,800), link=:x)
			figname = figdir*"/dvv_averaged_$(strfreq)Hz.png"
			Plots.savefig(p12, figname)
		end
	end

	if InputDict["Isoutputccstats"]
		statsdir = figdir*"/../fig_ccstats"
		if ispath(statsdir); rm(statsdir, recursive=true); end
		mkdir(statsdir)
		println("---Compute cc stats---")

		for (it, ts) = enumerate(tvec[1:InputDict["statsplotspan"]:end])

			tvec_stamp = timestamp.(ts)
			stamplabel = tvec_stamp[1:10]

			dvv_cc_t = dvv_cc_all[it]

			jldopen(statsdir*"/dvv_averaged_stats_$(stamplabel).jld2", "a+") do fo
				fo["freqband"]  = freqband
				fo["timestamp"] = ts
				fo["dvv_cc_t"]  = dvv_cc_t
			end


			for ifreq = 1:Nfreqband
				freqmin = freqband[ifreq]
				freqmax = freqband[ifreq+1]
				dvv_cc = dvv_cc_t[:, ifreq]
				stationpairnum = length(dvv_cc)
				# store tvec and dvv
				strfreq = @sprintf("%4.2f-%4.2f", round(freqmin, digits=2), round(freqmax, digits=2))

				# plot histogram of unittime cc for all station pairs
				ph = Plots.histogram(dvv_cc, bins = 42,
				 				normalize=false,
								label = "",
								xlabel = "Correlation coefficient with stretching",
								xlims = (0.0, 1.0),
								ylabel = "Count",
								fillcolor = :blue,
								fillalpha = 0.5,
								title = stamplabel*":"*strfreq*"Hz"*" $(stationpairnum) pairs",
								size=(800,600))
				figname = statsdir*"/dvv_averaged_cchist_$(stamplabel)_$(strfreq)Hz.png"
				Plots.savefig(ph, figname)
			end
		end
	end




    println("Averaging dvv has been successfully done.")
    return nothing

end

function readdvvfile(dvvfiname::String)

	fi = try
		jldopen(dvvfiname, "r")
	catch
		return nothing;
	end

	C = fi["xcorr"]
	Cname			= C.name
	Ccomp			= C.comp
	T 			= unix2datetime.(fi["T"])
	dvv 		= fi["dvv"]
	dvvmode	= fi["dvvmode"]

	dvvdict = Dict("name" 	=> Cname,
				   "comp"	=> Ccomp,
				   "T" 		=> T,
				   "dvv" 	=> dvv,
				   "dvvmode" => dvvmode,
				   "ccdist" => fi["ccdist"])

	if fi["dvvmode"] == "Stretching"
		dvvdict["cc"] = fi["cc"]
	end
	close(fi)

	println("read $(dvvfiname) done.")

	return dvvdict
end

function map_average_dvv(timewindow::Tuple, InputDict::Dict, dvv_dict_all::AbstractArray)

	t1 = now()
	# read files and gather dvv and coherence
	twin_left = timewindow[1]
	twin_right = timewindow[2]
	buffer_twin = Dates.Second(1.0) # [s]: shift time window backward to avoid including timestamp at right edge

	freqband = InputDict["freqband"]
	Nfreqband = length(freqband) - 1

	dvvfiname	= InputDict["dvvfiname"]

	dvv_all = Array{Float32, 2}(undef, 0, Nfreqband)

	dvv_cc_all = Array{Float32, 2}(undef, 0, Nfreqband)

	for dvv_dict in dvv_dict_all

		if isnothing(dvv_dict)
			continue
		end

		#C 	= dvv_dict["xcorr"]
		Cname = dvv_dict["name"]
		Ccomp = dvv_dict["comp"]

		# parse metadata
		stn1 = join(split(Cname, ".")[1:4], ".")
		stn2 = join(split(Cname, ".")[5:8], ".")
		ct = get_corrtype([stn1, stn2])
		comp = Ccomp

		if any(occursin.(ct, InputDict["corrtype"])) && any(occursin.(comp, InputDict["components"]))
			# this station pair has selected corrtype and component pair.

			tind = findall(x -> (x >= twin_left-buffer_twin && x <= twin_right - buffer_twin), dvv_dict["T"])
			dvv = dvv_dict["dvv"]
			ccdist 	= dvv_dict["ccdist"]

			if dvv_dict["dvvmode"] == "Stretching"
				cc 	= dvv_dict["cc"]
			end

			for tt in tind

				dvv_temp = dvv[tt, :]

				# apply thresholding criteria
				if dvv_dict["dvvmode"] == "Stretching"
					cc_temp  = cc[tt, :]
					ccdist_temp = ccdist[tt, :]
					# append to cc all
					dvv_cc_all = cat(dvv_cc_all, reshape(cc_temp, 1, Nfreqband), dims = 1)

					for ifreq = 1:Nfreqband
						#if cc_temp[ifreq] < InputDict["ccthreshold"]
						if (cc_temp[ifreq] < InputDict["ccthreshold"]) || (ccdist_temp > InputDict["ccdistance_threshold"])
							# filter out with small cc
							dvv_temp[ifreq] = NaN
							# DEBUG:
							#println(Cname*Ccomp)

						end
					end
				end

				# append to overall dvv
				dvv_all = cat(dvv_all, reshape(dvv_temp, 1, Nfreqband), dims=1)
			end
		end
	end


	# compute mean, std and count of non-NaN dvv
	dvv_mean_all = zeros(Float32, Nfreqband)
	dvv_std_all = zeros(Float32, Nfreqband)
	dvv_count_all = zeros(Float32, Nfreqband)

	for ifreq = 1:Nfreqband

		dvv_ifreq = dvv_all[:, ifreq]
		#count non-NaN value
		dvv_count_all[ifreq] 	= count(x -> !isnan(x), dvv_ifreq)
		dvv_nonan_ind 			= findall(x -> !isnan(x), dvv_ifreq)

		if isempty(dvv_nonan_ind) ||  dvv_count_all[ifreq] < InputDict["minimumcount"]
			# there is no value with this frequency band at this time
			dvv_mean_all[ifreq] = NaN
			dvv_std_all[ifreq] 	= NaN
		else
			dvv_mean_all[ifreq] = mean(dvv_ifreq[dvv_nonan_ind])
			dvv_std_all[ifreq] 	= std(dvv_ifreq[dvv_nonan_ind])
		end
	end

	t2 = now()
	elapsedtime =  (t2 - t1).value / 1e3
	println(string(twin_left)*"-"*string(twin_right)*" is averaged with $(elapsedtime) second.")

	return (dvv_mean_all, dvv_std_all, dvv_count_all, dvv_cc_all)
end


"""

    corrtype(stnPair::Array{String, 1})

Determine correlation type (cross-correlation, auto-correlation, or cross-channel-correlation) using string slicing.

# Arguments
- `stnPair::Array{String, 1},`    : Station pair, e.g. ["BP.SMNB..BP1", "BP.SMNB..BP3"]

# Output
- `corrtype::String`    : correlation type, e.g. "xchancorr"

"""
function get_corrtype(stnpair::Array{String, 1})

    stn1 = join(split(stnpair[1], ".")[1:2], ".")
    cha1 = split(stnpair[1], ".")[4][end]
    stn2 = join(split(stnpair[2], ".")[1:2], ".")
    cha2 = split(stnpair[2], ".")[4][end]

    if (stn1 == stn2) && (cha1 == cha2)
        ct = "auto-achan"
    elseif (stn1 == stn2) && (cha1 != cha2)
        ct = "auto-xchan"
    elseif (stn1 != stn2) && (cha1 == cha2)
        ct = "cross-achan"
    elseif (stn1 != stn2) && (cha1 != cha2)
        ct = "cross-xchan"
    else
        warning("corrtype error with $(stn1).$(cha1)-$(stn2).$(cha2). at get_corrtype()")
    end

    return ct
#     if stnpair[1] == stnpair[2]
#         ct = "acorr"
#     # same station, different channel
# elseif (stnpair[1][end-3:end] != stnpair[2][end-3:end]) && (stnpair[1][1:end-3] == stnpair[2][1:end-3])
#         ct = "xchancorr"
#     # different station
#     else
#         ct = "xcorr"
#     end
    # return ct
end
