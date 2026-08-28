module DataLoader

using JSON3, Graphs

export NetworkGraph,
       load_network, load_demands, load_scenario, load_instance,
       json_id, node_name, vertex_from_json_id, edge_attributes,
       capacity_matrix, metric_matrix

struct NetworkGraph{G<:AbstractGraph}
    graph::G
    node_data::Vector{NamedTuple{(:json_id, :name), Tuple{Int, String}}}
    json_to_vertex::Dict{Int, Int}
    all_arcs::Vector{Tuple{Int, Int, Int, Float64, Float64}}
    arc_groups::Dict{Tuple{Int,Int}, Vector{Int}}
    unique_edge_data::Dict{Tuple{Int, Int}, NamedTuple{(:id, :metric, :capacity), Tuple{Int, Float64, Float64}}}
end

function load_network(filename::String)
    data = JSON3.read(read(filename, String))
    hasproperty(data, :directed)   || error("Missing JSON field: directed")
    hasproperty(data, :multigraph) || error("Missing JSON field: multigraph")
    hasproperty(data, :nodes)      || error("Missing JSON field: nodes")
    hasproperty(data, :links)      || error("Missing JSON field: links")

    n = length(data.nodes)
    json_to_vertex = Dict{Int, Int}()
    node_data = Vector{NamedTuple{(:json_id, :name), Tuple{Int, String}}}(undef, n)

    for (v, node) in enumerate(data.nodes)
        id = Int(node.id)
        haskey(json_to_vertex, id) && error("Duplicate node id: $id")
        json_to_vertex[id] = v
        node_data[v] = (json_id = id, name = String(node.name))
    end

    is_directed = Bool(data.directed)
    g = is_directed ? SimpleDiGraph(n) : SimpleGraph(n)

    unique_edge_data = Dict{Tuple{Int, Int}, NamedTuple{(:id, :metric, :capacity), Tuple{Int, Float64, Float64}}}()
    all_arcs = Vector{Tuple{Int, Int, Int, Float64, Float64}}()
    arc_groups = Dict{Tuple{Int, Int}, Vector{Int}}()

    for link in data.links
        from_id = Int(link.from)
        to_id   = Int(link.to)

        haskey(json_to_vertex, from_id) || error("Unknown source id: $from_id")
        haskey(json_to_vertex, to_id)   || error("Unknown target id: $to_id")

        u = json_to_vertex[from_id]
        v = json_to_vertex[to_id]

        id = Int(link.id)
        metric = Float64(link.metric)
        capacity = Float64(link.capacity)

        push!(all_arcs, (u, v, id, metric, capacity))

        key = is_directed ? (u, v) : minmax(u, v)
        if !haskey(arc_groups, key)
            arc_groups[key] = Int[]
        end
        push!(arc_groups[key], length(all_arcs))

        if !haskey(unique_edge_data, key)
            add_edge!(g, u, v)
            unique_edge_data[key] = (id = id, metric = metric, capacity = capacity)
        end
    end

    return NetworkGraph(g, node_data, json_to_vertex, all_arcs, arc_groups, unique_edge_data)
end

json_id(net, v) = net.node_data[v].json_id
node_name(net, v) = net.node_data[v].name
vertex_from_json_id(net, id) = net.json_to_vertex[Int(id)]

function edge_attributes(net, u, v)
    key = is_directed(net.graph) ? (Int(u), Int(v)) : minmax(Int(u), Int(v))
    return net.unique_edge_data[key]
end

function capacity_matrix(net)
    n = nv(net.graph)
    C = zeros(n, n)
    for (key, data) in net.unique_edge_data
        u, v = key
        C[u, v] = data.capacity
    end
    return C
end

function metric_matrix(net)
    g = net.graph
    n = nv(g)
    D = fill(Inf, n, n)
    for v in vertices(g)
        D[v, v] = 0.0
    end
    for (key, data) in net.unique_edge_data
        u, v = key
        w = data.metric
        w < 0 && error("Dijkstra requires non-negative metrics.")
        D[u, v] = w
        if !is_directed(g)
            D[v, u] = w
        end
    end
    return D
end

function load_demands(filename::String)
    data = JSON3.read(read(filename, String))
    hasproperty(data, :num_time_slots) || error("Missing num_time_slots field")
    hasproperty(data, :demands)        || error("Missing demands field")

    demands = Vector{NamedTuple{(:s, :t, :vol0), Tuple{Int, Int, Float64}}}()
    for (idx, d) in enumerate(data.demands)
        hasproperty(d, :s) || error("Demand $idx missing 's'")
        hasproperty(d, :t) || error("Demand $idx missing 't'")
        hasproperty(d, :v) || error("Demand $idx missing 'v'")
        v = d.v
        if length(v) < 1
            error("Demand $idx has no volume data for time slot 0")
        end
        push!(demands, (s = Int(d.s), t = Int(d.t), vol0 = Float64(v[1])))
    end
    return demands
end

function load_scenario(filename::String)
    data = JSON3.read(read(filename, String))
    hasproperty(data, :max_segments) || error("Missing max_segments field")
    return (max_segments = Int(data.max_segments),)
end

function load_instance(prefix::String)
    net = load_network("$prefix-net.json")
    demands = load_demands("$prefix-tm.json")
    scenario = load_scenario("$prefix-scenario.json")
    return net, demands, scenario.max_segments
end

end # module