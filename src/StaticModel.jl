module StaticModel

using JuMP, HiGHS, Graphs, Main.DataLoader, Main.ECMP

export solve_static

function solve_static(net, demands, maxSeg, r_arc, time_limit)

    n = nv(net.graph)
    n_demands = length(demands)

    all_arcs = net.all_arcs
    n_arcs = length(all_arcs)

    # Get all segments from r_arc
    seg_set = Set{Tuple{Int,Int}}()
    for (arc, segdict) in r_arc
        for (seg, _) in segdict
            push!(seg_set, seg)
        end
    end
    segments = collect(seg_set)
    seg_index = Dict(seg => idx for (idx, seg) in enumerate(segments))
    n_segments = length(segments)

    if n_segments == 0
        @warn "No active segments found! Check ECMP computation."
        return :INFEASIBLE, 0.0, 0.0, nothing
    end

    model = Model(HiGHS.Optimizer)
    set_time_limit_sec(model, time_limit)
    set_silent(model)

    x = Dict{Tuple{Int,Int}, VariableRef}()
    for d in 1:n_demands
        for (idx, (i,j)) in enumerate(segments)
            x[(d, idx)] = @variable(model, binary=true, base_name="x_$(d)_$(i)_$(j)")
        end
    end

    for d in 1:n_demands
        s = vertex_from_json_id(net, demands[d].s)
        t = vertex_from_json_id(net, demands[d].t)

        for v in vertices(net.graph)
            incoming = AffExpr(0.0)
            for i in vertices(net.graph)
                if haskey(seg_index, (i, v))
                    idx = seg_index[(i, v)]
                    add_to_expression!(incoming, 1.0, x[(d, idx)])
                end
            end

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

    for d in 1:n_demands
        expr = AffExpr(0.0)
        for idx in 1:n_segments
            add_to_expression!(expr, 1.0, x[(d, idx)])
        end
        @constraint(model, expr <= maxSeg)
    end

    load = Dict{Int, VariableRef}()
    for arc_idx in 1:n_arcs
        load[arc_idx] = @variable(model, lower_bound=0, base_name="load_arc_$(arc_idx)")
    end

    for (unique_arc_key, segdict) in r_arc
        arc_indices = get(net.arc_groups, unique_arc_key, Int[])
        isempty(arc_indices) && continue

        expr = AffExpr(0.0)
        for ((i,j), frac) in segdict
            seg_idx = seg_index[(i,j)]
            for d in 1:n_demands
                vol = demands[d].vol0
                coeff = vol * frac
                add_to_expression!(expr, coeff, x[(d, seg_idx)])
            end
        end

        for arc_idx in arc_indices
            @constraint(model, load[arc_idx] == expr)
            u, v, id, metric, cap = all_arcs[arc_idx]
            @constraint(model, load[arc_idx] <= cap)
        end
    end

    L = @variable(model, lower_bound=0, base_name="L")
    for arc_idx in 1:n_arcs
        @constraint(model, load[arc_idx] <= L)
    end
    @objective(model, Min, L)

    optimize!(model)

    status = termination_status(model)
    obj_val = (status == OPTIMAL || status == ALMOST_OPTIMAL) ? value(L) : 0.0
    solve_time = JuMP.solve_time(model)
    return status, obj_val, solve_time, model
end

end # module