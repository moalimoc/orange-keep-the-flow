module StaticModel

using JuMP, HiGHS

# Load DataLoader as a local module
include("DataLoader.jl")
using .DataLoader

# Load ECMP as a local module (if you need functions from it)
include("ECMP.jl")
using .ECMP

# ------------------------------------------------------------------------------
# Helper: get list of arcs with non-zero r coefficients
# ------------------------------------------------------------------------------

function get_active_arcs(r_arc::Dict)
    return collect(keys(r_arc))
end

# ------------------------------------------------------------------------------
# Helper: get all segments (i,j) that appear in r_arc
# ------------------------------------------------------------------------------

function get_active_segments(r_arc::Dict)
    seg_set = Set{Tuple{Int,Int}}()
    for (arc, segdict) in r_arc
        for (seg, _) in segdict
            push!(seg_set, seg)
        end
    end
    return collect(seg_set)
end

# ------------------------------------------------------------------------------
# Build and solve the static model (time period 0 only)
# ------------------------------------------------------------------------------

function solve_static(net::NetworkGraph, demands::Vector, maxSeg::Int,
                      r_arc::Dict, time_limit::Float64)

    n = nv(net.graph)
    n_demands = length(demands)

    # Get arcs and segments
    arcs = get_active_arcs(r_arc)
    segments = get_active_segments(r_arc)
    seg_index = Dict(seg => idx for (idx, seg) in enumerate(segments))
    n_segments = length(segments)

    if n_segments == 0
        @warn "No active segments found! Check ECMP computation."
        return :INFEASIBLE, 0.0, 0.0, nothing
    end

    # Build model
    model = Model(HiGHS.Optimizer)
    set_time_limit_sec(model, time_limit)
    set_silent(model)

    # --- Variables: x[d, seg_idx] binary ---
    x = Dict{Tuple{Int,Int}, VariableRef}()
    for d in 1:n_demands
        for (idx, (i,j)) in enumerate(segments)
            x[(d, idx)] = @variable(model, binary=true, base_name="x_$(d)_$(i)_$(j)")
        end
    end

    # --- Flow conservation constraints (1) ---
    for d in 1:n_demands
        s = vertex_from_json_id(net, demands[d].s)
        t = vertex_from_json_id(net, demands[d].t)

        for v in vertices(net.graph)
            # Incoming: sum over i where (i,v) exists
            incoming = AffExpr(0.0)
            for i in vertices(net.graph)
                if haskey(seg_index, (i, v))
                    idx = seg_index[(i, v)]
                    add_to_expression!(incoming, 1.0, x[(d, idx)])
                end
            end

            # Outgoing: sum over j where (v,j) exists
            outgoing = AffExpr(0.0)
            for j in vertices(net.graph)
                if haskey(seg_index, (v, j))
                    idx = seg_index[(v, j)]
                    add_to_expression!(outgoing, 1.0, x[(d, idx)])
                end
            end

            rhs = 0.0
            if v == s
                rhs = 1.0
            elseif v == t
                rhs = -1.0
            end

            @constraint(model, incoming - outgoing == rhs)
        end
    end

    # --- Segment limit constraint (2) ---
    for d in 1:n_demands
        expr = AffExpr(0.0)
        for idx in 1:n_segments
            add_to_expression!(expr, 1.0, x[(d, idx)])
        end
        @constraint(model, expr <= maxSeg)
    end

    # --- Load variables ---
    load = Dict{Tuple{Int,Int}, VariableRef}()
    for arc in arcs
        load[arc] = @variable(model, lower_bound=0, base_name="load_$(arc)")
    end

    # --- Load constraints (3) ---
    for arc in arcs
        expr = AffExpr(0.0)
        segdict = r_arc[arc]

        for ((i,j), frac) in segdict
            seg_idx = seg_index[(i,j)]
            for d in 1:n_demands
                vol = demands[d].vol0
                coeff = vol * frac
                add_to_expression!(expr, coeff, x[(d, seg_idx)])
            end
        end

        @constraint(model, load[arc] == expr)
    end

    # --- Objective: minimize maximum load (a surrogate for lexicographic) ---
    L = @variable(model, lower_bound=0, base_name="L")
    for arc in arcs
        @constraint(model, load[arc] <= L)
    end
    @objective(model, Min, L)

    # --- Solve ---
    optimize!(model)

    status = termination_status(model)
    obj_val = 0.0
    if status == OPTIMAL || status == ALMOST_OPTIMAL
        obj_val = value(L)
    end

    solve_time = solve_time(model)

    # Optionally, you can retrieve the solution values for x and load
    # For now, we just return the key results

    return status, obj_val, solve_time, model
end

# ------------------------------------------------------------------------------
# A more advanced lexicographic objective (commented out for simplicity)
# ------------------------------------------------------------------------------

# function solve_lexicographic(...)
#     # Sequential solves: min L1, fix, min L2, etc.
#     # This is left as an extension
# end

end # module