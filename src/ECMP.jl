module ECMP

using Graphs
include("DataLoader.jl")
using .DataLoader

# ------------------------------------------------------------------------------
# 1. All-pairs shortest paths with path counts
# ------------------------------------------------------------------------------

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
# 2. Compute r for all pairs and arcs
# ------------------------------------------------------------------------------

function compute_r_all_pairs(net::NetworkGraph)
    g = net.graph
    n = nv(g)
    D = metric_matrix(net)
    dist, pathcount = all_pairs_shortest_data(g, D)

    arcs = [(src(e), dst(e)) for e in edges(g)]
    arc_set = Set(arcs)

    # r_seg[(i,j)] = Dict((u,v) => fraction)
    r_seg = Dict{Tuple{Int,Int}, Dict{Tuple{Int,Int},Float64}}()
    # r_arc[(u,v)] = Dict((i,j) => fraction)
    r_arc = Dict{Tuple{Int,Int}, Dict{Tuple{Int,Int},Float64}}()

    # For each pair (i,j) that is reachable
    for i in vertices(g)
        for j in vertices(g)
            if i == j
                continue
            end
            if !isfinite(dist[i, j])
                continue
            end
            sigma_ij = pathcount[i, j]
            if sigma_ij == 0
                continue
            end

            # For each arc (u,v), check if it lies on a shortest path from i to j
            seg_dict = Dict{Tuple{Int,Int},Float64}()
            for (u,v) in arcs
                if isfinite(dist[i, u]) && isfinite(dist[v, j])
                    if isapprox(dist[i, u] + D[u, v] + dist[v, j], dist[i, j]; atol=1e-9)
                        frac = (pathcount[i, u] * pathcount[v, j]) / sigma_ij
                        if frac > 0
                            seg_dict[(u,v)] = frac
                        end
                    end
                end
            end

            if !isempty(seg_dict)
                r_seg[(i,j)] = seg_dict
                # Also fill r_arc
                for ((u,v), frac) in seg_dict
                    if !haskey(r_arc, (u,v))
                        r_arc[(u,v)] = Dict{Tuple{Int,Int},Float64}()
                    end
                    r_arc[(u,v)][(i,j)] = frac
                end
            end
        end
    end

    return r_seg, r_arc
end

# ------------------------------------------------------------------------------
# 3. Helper to get list of reachable segments
# ------------------------------------------------------------------------------

function reachable_segments(net::NetworkGraph)
    g = net.graph
    pairs = [(i,j) for i in vertices(g) for j in vertices(g) if i != j && has_path(g, i, j)]
    return pairs
end

# ------------------------------------------------------------------------------
# 4. Optional: compute r for specific demands only (older version)
# ------------------------------------------------------------------------------

# ... (keep the old compute_r if needed, but we'll use the new one)

end # module