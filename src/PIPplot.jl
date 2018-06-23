module PIPplot
using Reexport

using JuLIP, NBodyIPs, DataFrames, Plots, JLD, StaticArrays
using NBodyIPs.Data: Dat
using Plots; gr()

using PyCall
unshift!(PyVector(pyimport("sys")["path"]), Pkg.dir()*"/NBodyIPs/src")

@pyimport plottools
@pyimport numpy
@pyimport matplotlib.pyplot as plt

export PIPpotplot, PIPenergyplot, PIPforceplot, PIPenergyplot2, PIPforceplot2, PIP3Body3D, PIPfullplot

function unfold(A)
    V = []
    for x in A
        if x === A
            push!(V, x)
        else
            append!(V, unfold(x))
        end
    end
    V
end

function PIPpotplot(potname, IP)
    e0 = string(IP.orders[1].c)
    cutoff2b = string(IP.orders[3].D.rcut)
    cutoff3b = string(IP.orders[4].D.rcut)
    cutoff4b = string(IP.orders[5].D.rcut)
	r0 = 2.89
	rr = linspace(2, 9.0, 200)
	P = plot(rr, IP.(rr); yaxis = ([-0.1, 0.2],), label="2b");
	plot!(P, rr, IP.(rr, rr, rr), label="3b")
	W4 = rr -> IP(SVector(2.9, rr, 2.9, 2.9, 2.9, 2.9))
	title = "Potname=" * potname[3:end-4] * " rc=(" * cutoff2b * "/" * cutoff3b * "/" * cutoff4b * ")"*" e0=" * e0
	plot!(P, rr, W4.(rr), title=title, xlabel="Interatomic distance (Å)", titlefont=8, label="4b", legend=true)
    vline!([2.89], label="r0", color="black")
end

function PIPenergyplot(test_data, IP, s; e_data = 1, returnplot = true)
    d=Dict([(i,test_data[i].config_type) for i in 1:length(test_data)])
    configselection = [[i,d[i]] for i in 1:s:length(test_data)]
    configtypes = unique([configselection[i][2] for i in 1:length(configselection)])
    if e_data == 1
        e_data = hcat([[test_data[configselection[i][1]].E, energy(IP, Atoms(test_data[configselection[i][1]])), configselection[i][2]] for i in 1:length(configselection)]...)
    end
    RMSE_E = sqrt(mean((e_data[1,:]-e_data[2,:]).^2))/54
    if returnplot == true
        xmin = findmin(convert(Array{Float64,1},unfold(e_data[1,:])))[1]/54
        xmax = findmax(convert(Array{Float64,1},unfold(e_data[1,:])))[1]/54
        x = linspace(xmin, xmax, 1000)
        f(x) = x#(abs.(ymin - ymax)/abs.(xmin-xmax))x
        p = plot(x, f(x), color="black", legend=false)
        colors = ["red", "green", "blue", "yellow", "purple", "brown", "cyan", "violet"]
        for i in 1:length(configtypes)
            config = configtypes[i]
            indices = find(e_data[3, 1:end] .== config)
            target = convert(Array{Float64,1}, unfold(e_data[1,findmin(indices)[1]:findmax(indices)[1]]))/54
            test = convert(Array{Float64,1}, unfold(e_data[2,findmin(indices)[1]:findmax(indices)[1]]))/54
            data = transpose(numpy.array([target, test], numpy.float64))
            thin_data, thin_weights = plottools.thin_points(data, r=0.001)
            p = scatter!(p, thin_data[:,1], thin_data[:,2], markersize=4*(thin_weights).^(0.5), markeralpha=0.1, label=config, legend=:bottomright, xlabel="Target Energy per atom (eV) RMSE = "*string(RMSE_E)[1:8], ylabel="Predicted Energy per atom (eV)", legendfont=7, color=colors[i])  ##yaxis = (:log10, (1,100000)), markerstrokealpha=1, markerstrokecolor="black",
        end
        return p
    else
        return e_data
    end
end

function PIPforceplot(test_data, IP, s; f_data = 1, returnplot=true)
    d=Dict([(i,test_data[i].config_type) for i in 1:length(test_data)])
    configselection = [[i,d[i]] for i in 1:s:length(test_data)] ###which different ones are there
    configtypes = unique([configselection[i][2] for i in 1:length(configselection)])
    if f_data == 1
        f_data = hcat([[unfold(test_data[configselection[i][1]].F), unfold(forces(IP, Atoms(test_data[configselection[i][1]]))), configselection[i][2]] for i in 1:length(configselection)]...)
    end
    RMSE_F = sqrt(mean((unfold(f_data[1,:]) - unfold(f_data[2,:])).^2))
    if returnplot == true
        xmin = findmin(convert(Array{Float64,1},unfold(f_data[1,:])))[1]
        xmax = findmax(convert(Array{Float64,1},unfold(f_data[1,:])))[1]
        x = linspace(xmin, xmax, 1000)
        f(x) = x
        p = plot(x, f(x), color="black", legend=false)
        colors = ["red", "green", "blue", "yellow", "purple", "brown", "cyan", "violet"]
        for i in 1:length(configtypes)
            config = configtypes[i]
            indices = find(f_data[3, 1:end] .== config)
            target = convert(Array{Float64,1}, unfold(f_data[1,findmin(indices)[1]:findmax(indices)[1]]))
            test = convert(Array{Float64,1}, unfold(f_data[2,findmin(indices)[1]:findmax(indices)[1]]))
            data = transpose(numpy.array([test, target], numpy.float64))
            thin_data, thin_weights = plottools.thin_points(data, r=0.02)
            p = scatter!(p, thin_data[:,1], thin_data[:,2], markersize=4*(thin_weights).^(0.5), markeralpha=0.1, label=config, legend=:bottomright, xlabel="Target Force (ev/Å) RMSE = "*string(RMSE_F)[1:8], ylabel="Predicted Force (ev/Å)", legendfont=7, color=colors[i], legend=:false)  ##markerstrokealpha=1, markerstrokecolor="black",
        end
        return p
    else
        return f_data
    end
end

function PIP3Body3D(potname, IP; p2b = 2.25, p2e = 4.5, p2s = 50, theta = 60)
	V3 = IP
	W3 = (r1, r2, r3) -> V3(SVector(r1, r2, r3))
	elist = []
	R = []
	for r in linspace(p2b,p2e,p2s)
	    for r12 in linspace(p2b,p2e,p2s)
	        x = [r12*cosd(theta), r12*sind(theta)]
	        y = [r, 0]
	        z = x - y
	        push!(R, [norm(x),norm(y),norm(z)])
	    end
	end
	for i in 1:length(R)
		e = W3.(R[i][1], R[i][2], R[i][3])
		push!(elist, e)
	end
	Z = Matrix(0,p2s)
	for i in 1:p2s:(p2s*p2s)
	    v = transpose(vec(elist[i:i+(50-1)]))
	    Z = cat(1,Z,v)
	end
	X, Y = linspace(p2b,p2e,p2s), linspace(p2b,p2e,p2s)
	Z = convert(Array{Float64,2},Z)
	heatmap(X,Y,Z, xlabel="Interatomic distance r12 (Å)", ylabel="Interatomic distance r13 (Å)")
    plot!([X[1], X[end]], [Y[1], Y[end]], color="white", legend=false)
end

function PIPforceplot2(test_data, IP, s; f_data = 1, returnplot=true)
    d=Dict([(i,test_data[i].config_type) for i in 1:length(test_data)])
    configselection = [[i,d[i]] for i in 1:s:length(test_data)] ###which different ones are there
    configtypes = unique([configselection[i][2] for i in 1:length(configselection)])
    if f_data == 1
        f_data = hcat([[unfold(test_data[configselection[i][1]].F), unfold(forces(IP, Atoms(test_data[configselection[i][1]]))), configselection[i][2]] for i in 1:length(configselection)]...)
    end
    RMSE_F = sqrt(mean((unfold(f_data[1,:]) - unfold(f_data[2,:])).^2))
    if returnplot == true
        xmin = findmin(convert(Array{Float64,1},unfold(f_data[1,:])))[1]
        xmax = findmax(convert(Array{Float64,1},unfold(f_data[1,:])))[1]
        ymin = findmin(convert(Array{Float64,1},unfold(f_data[2,:])))[1]
        ymax = findmax(convert(Array{Float64,1},unfold(f_data[2,:])))[1]
        p = plot()
        colors = ["red", "green", "blue", "yellow", "purple", "brown", "cyan", "violet"]
        for i in 1:length(configtypes)
            config = configtypes[i]
            indices = find(f_data[3, 1:end] .== config)
            target = convert(Array{Float64,1}, unfold(f_data[1,findmin(indices)[1]:findmax(indices)[1]]))
            test = convert(Array{Float64,1}, unfold(f_data[2,findmin(indices)[1]:findmax(indices)[1]]))
            error = abs.(test-target)
            data = transpose(numpy.array([target, error], numpy.float64))
            thin_data, thin_weights = plottools.thin_points(data, r=0.02)
            p = scatter!(p, thin_data[:,1], thin_data[:,2]+1e-10, markersize=4*(thin_weights).^(0.5), markeralpha=0.1, label=config, legend=:bottomright, xlabel="Target Force (ev/Å) RMSE = "*string(RMSE_F)[1:8], ylabel="| Predicted Force Error | (ev/Å)", yaxis = (:log10, (0.00005,2)), legendfont=7, color=colors[i], legend=:false)  ##markerstrokealpha=1, markerstrokecolor="black",
        end
        return p
    else
        return f_data
    end
end

function PIPenergyplot2(test_data, IP, s; e_data = 1, returnplot = true)
    d=Dict([(i,test_data[i].config_type) for i in 1:length(test_data)])
    configselection = [[i,d[i]] for i in 1:s:length(test_data)]
    configtypes = unique([configselection[i][2] for i in 1:length(configselection)])
    if e_data == 1
        e_data = hcat([[test_data[configselection[i][1]].E, energy(IP, Atoms(test_data[configselection[i][1]])), configselection[i][2]] for i in 1:length(configselection)]...)
    end
    RMSE_E = sqrt(mean((e_data[1,:]-e_data[2,:]).^2))/54
    if returnplot == true
        p = plot()
        colors = ["red", "green", "blue", "yellow", "purple", "brown", "cyan", "violet"]
        for i in 1:length(configtypes)
            config = configtypes[i]
            indices = find(e_data[3, 1:end] .== config)
            target = convert(Array{Float64,1}, unfold(e_data[1,findmin(indices)[1]:findmax(indices)[1]]))/54
            test = convert(Array{Float64,1}, unfold(e_data[2,findmin(indices)[1]:findmax(indices)[1]]))/54
            error = abs.(test-target)
            data = transpose(numpy.array([target, error], numpy.float64))
            thin_data, thin_weights = plottools.thin_points(data, r=0.01)
            p = scatter!(p, thin_data[:,1], thin_data[:,2]+1e-10, markersize=4*(thin_weights).^(0.5), markeralpha=0.1, label=config, legend=:bottomright, xlabel="Target Energy per atom (eV) RMSE = "*string(RMSE_E)[1:8], yaxis = (:log10, (0.00005,0.04)), ylabel="Predicted Energy per atom (eV)", legendfont=7, color=colors[i], legend=:false)  ##yaxis = (:log10, (1,100000)), markerstrokealpha=1, markerstrokecolor="black",
        end
        return p
    else
        return e_data
    end
end

function PIPfullplot(potname, test_data, IP; s=10)
    if isfile("./data_"*potname[3:end-4]*".jld") == true
        file = load("./data_"*potname[3:end-4]".jld")
        e_data = file["e_data"]
        f_data = file["f_data"]
    else
        e_data = PIPenergyplot(test_data, IP, s, returnplot=false)
        f_data = PIPforceplot(test_data, IP, s, returnplot=false)
        file = jldopen("./data_"*potname[3:end-4]".jld", "w")
        write(file, "e_data", e_data)
        write(file, "f_data", f_data)
        close(file)
    end
    p1 = PIPpotplot(potname, IP)
    p2 = PIP3Body3D(potname, IP)
    p3 = PIPenergyplot(test_data, IP, e_data=e_data, s)
    p4 = PIPforceplot(test_data, IP, f_data=f_data, s)
    p5 = PIPenergyplot2(test_data, IP, e_data=e_data, s)
    p6 = PIPforceplot2(test_data, IP, f_data=f_data, s)
    p7 = plot(p1,p2,p3,p4,p5,p6, layout=(3,2), size=(900,1500))
    #png(p7, "plot")
end

end # module
