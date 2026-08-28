module ECMP

using Graphs, Main.DataLoader   # <-- Main.DataLoader

export compute_r_all_pairs

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

function compute_r_all_pairs(net)
    g = net.graph
    n = nv(g)
    D = metric_matrix(net)
    dist, pathcount = all_pairs_shortest_data(g, D)

    unique_arcs = collect(keys(net.unique_edge_data))

    r_seg = Dict{Tuple{Int,Int}, Dict{Tuple{Int,Int},Float64}}()
    r_arc = Dict{Tuple{Int,Int}, Dict{Tuple{Int,Int},Float64}}()

    for i in vertices(g)
        for j in vertices(g)
            i == j && continue
            !isfinite(dist[i, j]) && continue
            sigma_ij = pathcount[i, j]
            sigma_ij == 0 && continue

            seg_dict = Dict{Tuple{Int,Int},Float64}()
            for (u, v) in unique_arcs
                if isfinite(dist[i, u]) && isfinite(dist[v, j])
                    if isapprox(dist[i, u] + D[u, v] + dist[v, j], dist[i, j]; atol=1e-9)
                        frac = (pathcount[i, u] * pathcount[v, j]) / sigma_ij
                        if frac > 0
                            seg_dict[(u, v)] = frac
                        end
                    end
                end
            end

            if !isempty(seg_dict)
                r_seg[(i, j)] = seg_dict
                for ((u, v), frac) in seg_dict
                    if !haskey(r_arc, (u, v))
                        r_arc[(u, v)] = Dict{Tuple{Int,Int},Float64}()
                    end
                    r_arc[(u, v)][(i, j)] = frac
                end
            end
        end
    end

    return r_seg, r_arc
end

end # module