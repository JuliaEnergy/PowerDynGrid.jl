begin
    ENV["PLOTS_USE_ATOM_PLOTPANE"] = "false"
    using Pkg
    Pkg.activate(@__DIR__)
    cd(@__DIR__)
end

begin
    using CSV
    using DataFrames
    using DPSA
    using DPSABase: AbstractNodeParameters
    using DPSAGridTools
    using LaTeXStrings
    using Plots
    using SparseArrays
end

begin
    lines_df = CSV.read("IEEE14_lines.csv")
    df_names = [:from, :to, :R, :X, :charging, :tap_ratio]
    names!(getfield(lines_df, :colindex), df_names)

    busses_df = CSV.read("IEEE14_busses.csv")[[2,5,6,7,8,9,10]]
    df_names = [:type, :P_gen, :Q_gen, :P_load, :Q_load, :intertia, :damping]
    names!(getfield(busses_df, :colindex), df_names)
end
#busses_df

begin
    println("converting dataframes to types of DPSA")
    function busdf2nodelist(busses_df)
        node_list = Array{AbstractNodeParameters,1}()
        for bus_index = 1:size(busses_df)[1]
            bus_data = busses_df[bus_index,:]
            if busses_df[bus_index,:type] == "S"
                append!(node_list, [SlackAlgebraic(U=1)])
            elseif busses_df[bus_index,:type] == "G"
                append!(node_list, [SwingEqLVS(
                    H=busses_df[bus_index,:intertia] * 100  ,
                    P=(busses_df[bus_index,:P_gen] - busses_df[bus_index,:P_load]),
                    D=busses_df[bus_index,:damping],
                    Ω=50,
                    Γ= 2,
                    V=1
                )])
            elseif busses_df[bus_index,:type] == "L"
                append!(node_list, [PQAlgebraic(
                    S= -busses_df[bus_index,:P_load] - im*busses_df[bus_index,:Q_load]
                )])
            end
        end
        return node_list
    end

    node_list = busdf2nodelist(busses_df)
end

begin
    function linedf2LY(lines_df, num_nodes)
        #Y = spzeros(Complex, num_nodes, num_nodes)
        branch_list = Array{DPSAGridTools.TwoPort,1}()
        for line_index = 1:size(lines_df)[1]
            from = lines_df[line_index,:from]
            to = lines_df[line_index,:to]
            tap = lines_df[line_index,:tap_ratio]

            if (from > num_nodes) || (to > num_nodes)
                warn("Skipping line $line_index from $from to $(to)!")
                continue
            end

            if isequal(tap, 1)
                append!(branch_list, [line(
                    from,
                    to,
                    lines_df[line_index,:R],
                    lines_df[line_index,:X],
                    lines_df[line_index,:charging]
                )])
            else
                append!(branch_list, [trafo(
                    from,
                    to,
                    lines_df[line_index,:R],
                    lines_df[line_index,:X],
                    tap
                )])
            end

            # admittance = 1/(lines_df[line_index,:R] + im*lines_df[line_index,:X])
            # println("$from --> $to : $admittance")
            # Y[from, to] = - admittance
            # Y[to, from] = - admittance
            # Y[from, from] += admittance # note the +=
            # Y[to, to] += admittance # note the +=
        end
        return build_nodal_admittance(num_nodes, branch_list)
    end

    # admittance laplacian
    LY = linedf2LY(lines_df, length(node_list))
end


################################################################
# plotting the network representing the power grid
# check the two below for plotting graphs
# using LightGraphs
# using GraphPlot
# g = Graph(Array(LY).!=0)
# gplot(g)
################################################################

# create network dynamics object
g = GridDynamics(node_list, LY, skip_LY_check=true)
# search for fixed point

# find the fixed point = normal operation point
fp = operationpoint(g, ones(SystemSize(g)))

begin
    # define the initial condition as a perturbation from the fixed point
    x0 = copy(fp)
    x0[1, :int, 1] += 0.2 # perturbation on the ω of the first node
    #x0[n, :int, i] : access to the i-th internal variables of the n-th node
    #x0[n, :u] : access to the complex voltage of the n-th node
    #x0[n, :v] : access to the magnitude of the voltage of the n-th node
    #x0[n, :φ] : access to the voltage angle of the n-th node
    timespan = (0.0,10.0)
    # solve it
    sol = solve(g, x0, timespan);
end

Plots.scalefontsizes(0.8)
begin
    gr() # comment this out to use a different plotting backend
    swing_indices = findall(busses_df[:type] .== "G")
    ω_colors = reshape(Plots.get_color_palette(:auto, plot_color(:white), 8)[swing_indices], (1,length(swing_indices)))
    ω_labels = reshape([latexstring(string(raw"\omega", "_{$i}","[$(busses_df[i,:type])]")) for i=swing_indices], (1, length(swing_indices)))
    pl_v = plot(sol, :, :v, legend = (0.4, 1.), ylabel=L"V [p.u.]")
    pl_p = plot(sol, :, :p, legend = (0.6, 1.), ylabel=L"P [p.u.]")
    pl_ω = plot(sol, swing_indices, :int, 1, legend = (0.8, 1.), ylabel=L"\omega \left[\frac{rad}{s}\right]", label=ω_labels, color=ω_colors)
    pl = plot(pl_v, pl_p, pl_ω;
        layout=(3,1),
        size = (960, 540),
        lw=3,
        xlabel=L"t[s]"
        )
    savefig(pl, "ieee14.pdf")
    display(pl)
end
