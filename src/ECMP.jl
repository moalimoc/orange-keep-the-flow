module ECMP

using Graphs
include("DataLoader.jl")
using .DataLoader

# ------------------------------------------------------------------------------
# 1. All‑pairs shortest paths with path counts
# ------------------------------------------------------------------------------

"""
    all_pairs_shortest_data(g::AbstractGraph, distmx::AbstractMatrix)

Runs Dijkstra from every vertex with `allpaths = true` to obtain:
- `dist[i, j]`: shortest‑path distance from i to j
- `pathcount[i, j]`: number of shortest paths from i to j

Returns a tuple (dist, pathcount).
"""
function all_pairs_shortest_data(g::AbstractGraph, distmx::AbstractMatrix)
    n = nv(g)
    dist = fill(Inf, n, n)
    pathcount = zeros(Float64, n, n)

    for s in vertices(g)
        state = dijkstra_shortest_paths(g, s, distmx; allpaths = true)
        dist[s, :] = state.dists
        pathcount[s, :] = state.pathcounts
    end

    return dist, pathcount
end

# ------------------------------------------------------------------------------
# 2. Build arc lists and capacities
# ------------------------------------------------------------------------------

"""
    arcs_list(net::NetworkGraph) -> Vector{Tuple{Int,Int}}

Returns a vector of all directed arcs (u, v) in the graph.
"""
function arcs_list(net::NetworkGraph)
    g = net.graph
    return [(src(e), dst(e)) for e in edges(g)]
end

"""
    arc_capacity(net::NetworkGraph, u::Int, v::Int) -> Float64

Returns the capacity of arc (u,v). Assumes the arc exists.
"""
function arc_capacity(net::NetworkGraph, u::Int, v::Int)
    return edge_attributes(net, u, v).capacity
end

# ------------------------------------------------------------------------------
# 3. Main function: compute r for a given time period
# ------------------------------------------------------------------------------

"""
    compute_r(net::NetworkGraph, demands::Vector; t::Int=0) -> Dict

For each demand (s, t_dem) and each arc (u, v), compute the ECMP split coefficient
r[(demand_index, u, v)] = fraction of that demand's flow traversing arc (u,v).

The coefficient is computed for a specific time period `t`; for t>0, the graph
would have been modified (e.g., by removing downed links) before calling this function.
This version assumes the graph is already the correct one for that time period.

Returns a dictionary keyed by (demand_idx, u, v) with value Float64.
"""
function compute_r(net::NetworkGraph, demands::Vector; t::Int=0)
    g = net.graph
    n = nv(g)

    # Get metric matrix (weights)
    D = metric_matrix(net)

    # Run all‑pairs Dijkstra once
    dist, pathcount = all_pairs_shortest_data(g, D)

    # Prepare arc list
    arcs = arcs_list(net)
    arc_set = Set(arcs)  # for fast membership test

    # Result dictionary: key = (demand_idx, u, v) -> fraction
    r = Dict{Tuple{Int, Int, Int}, Float64}()

    # For each demand
    for (idx, dem) in enumerate(demands)
        s = vertex_from_json_id(net, dem.s)      # source in Julia vertex numbering
        t_dem = vertex_from_json_id(net, dem.t)  # target

        # If no path from s to t_dem, all r are zero (skip)
        if !isfinite(dist[s, t_dem])
            continue
        end

        sigma_st = pathcount[s, t_dem]
        if sigma_st == 0
            continue
        end

        # For each arc, check if it lies on a shortest path
        for (u, v) in arcs
            # Condition: dist[s,u] + weight(u,v) + dist[v,t] == dist[s,t]
            if isfinite(dist[s, u]) && isfinite(dist[v, t_dem])
                # Use isapprox with a tolerance for floating-point comparison
                if isapprox(dist[s, u] + D[u, v] + dist[v, t_dem], dist[s, t_dem]; atol=1e-9)
                    # Arc is on a shortest path from s to t_dem
                    sigma_su = pathcount[s, u]
                    sigma_vt = pathcount[v, t_dem]
                    coeff = (sigma_su * sigma_vt) / sigma_st
                    # Store the coefficient
                    r[(idx, u, v)] = coeff
                end
            end
        end
    end

    return r
end

# ------------------------------------------------------------------------------
# 4. Convenience: get r as a dense matrix per demand (for use in JuMP)
# ------------------------------------------------------------------------------

"""
    r_as_matrix(r::Dict, n_demands::Int, n::Int) -> Vector{Matrix{Float64}}

Converts the sparse dictionary r[(d,u,v)] into a vector of n×n matrices,
one per demand, where r[d][u,v] is the coefficient (0 if not present).
"""
function r_as_matrix(r::Dict, n_demands::Int, n::Int)
    matrices = [zeros(n, n) for _ in 1:n_demands]
    for ((d, u, v), val) in r
        matrices[d][u, v] = val
    end
    return matrices
end

end # module