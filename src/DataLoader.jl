module DataLoader

using JSON3, Graphs

# Export the public API
export NetworkGraph,
       load_network, load_demands, load_scenario, load_instance,
       json_id, node_name, vertex_from_json_id, edge_attributes,
       capacity_matrix, metric_matrix

# ------------------------------------------------------------------------------
# 1. Data structures
# ------------------------------------------------------------------------------

"""
    NetworkGraph{G<:AbstractGraph}

Stores the graph topology together with node and edge attributes.
"""
struct NetworkGraph{G<:AbstractGraph}
    graph::G
    node_data::Vector{NamedTuple{(:json_id, :name), Tuple{Int, String}}}
    json_to_vertex::Dict{Int, Int}
    edge_data::Dict{Tuple{Int, Int}, NamedTuple{(:id, :metric, :capacity), Tuple{Int, Float64, Float64}}}
end

# ------------------------------------------------------------------------------
# 2. Network loading
# ------------------------------------------------------------------------------

"""
    load_network(filename::String) -> NetworkGraph

Reads a JSON network file (same format as the ROADEF challenge) and returns
a NetworkGraph object.
"""
function load_network(filename::String)
    data = JSON3.read(read(filename, String))

    # --- Validate required fields ---
    hasproperty(data, :directed)   || error("Missing JSON field: directed")
    hasproperty(data, :multigraph) || error("Missing JSON field: multigraph")
    hasproperty(data, :nodes)      || error("Missing JSON field: nodes")
    hasproperty(data, :links)      || error("Missing JSON field: links")

    data.multigraph && error("This loader does not support multigraph=true.")

    n = length(data.nodes)

    # --- Build node mapping ---
    json_to_vertex = Dict{Int, Int}()
    node_data = Vector{NamedTuple{(:json_id, :name), Tuple{Int, String}}}(undef, n)

    for (v, node) in enumerate(data.nodes)
        id = Int(node.id)
        haskey(json_to_vertex, id) && error("Duplicate node id: $id")
        json_to_vertex[id] = v
        node_data[v] = (json_id = id, name = String(node.name))
    end

    # --- Build graph ---
    g = data.directed ? SimpleDiGraph(n) : SimpleGraph(n)

    edge_data = Dict{Tuple{Int, Int}, NamedTuple{(:id, :metric, :capacity), Tuple{Int, Float64, Float64}}}()

    for link in data.links
        from_id = Int(link.from)
        to_id   = Int(link.to)

        haskey(json_to_vertex, from_id) || error("Unknown source id: $from_id")
        haskey(json_to_vertex, to_id)   || error("Unknown target id: $to_id")

        u = json_to_vertex[from_id]
        v = json_to_vertex[to_id]

        key = data.directed ? (u, v) : minmax(u, v)

        haskey(edge_data, key) && error("Duplicate edge ($from_id,$to_id) when multigraph=false.")

        add_edge!(g, u, v) || error("Could not add edge ($from_id,$to_id).")

        edge_data[key] = (
            id       = Int(link.id),
            metric   = Float64(link.metric),
            capacity = Float64(link.capacity),
        )
    end

    return NetworkGraph(g, node_data, json_to_vertex, edge_data)
end

# ------------------------------------------------------------------------------
# 3. Access helpers
# ------------------------------------------------------------------------------

json_id(net::NetworkGraph, v::Integer) = net.node_data[v].json_id
node_name(net::NetworkGraph, v::Integer) = net.node_data[v].name
vertex_from_json_id(net::NetworkGraph, id::Integer) = net.json_to_vertex[Int(id)]

function edge_attributes(net::NetworkGraph, u::Integer, v::Integer)
    key = is_directed(net.graph) ? (Int(u), Int(v)) : minmax(Int(u), Int(v))
    return net.edge_data[key]
end

function capacity_matrix(net::NetworkGraph)
    n = nv(net.graph)
    C = zeros(n, n)
    for e in edges(net.graph)
        u, v = src(e), dst(e)
        C[u, v] = edge_attributes(net, u, v).capacity
    end
    return C
end

function metric_matrix(net::NetworkGraph)
    g = net.graph
    n = nv(g)
    D = fill(Inf, n, n)
    for v in vertices(g)
        D[v, v] = 0.0
    end
    for e in edges(g)
        u, v = src(e), dst(e)
        w = edge_attributes(net, u, v).metric
        w < 0 && error("Dijkstra requires non-negative metrics.")
        D[u, v] = w
        if !is_directed(g)
            D[v, u] = w
        end
    end
    return D
end

# ------------------------------------------------------------------------------
# 4. Demands loading
# ------------------------------------------------------------------------------

"""
    load_demands(filename::String) -> Vector{NamedTuple}

Loads traffic demands from a JSON file. Returns a vector of tuples:
(s, t, vol0) where s and t are JSON node IDs, and vol0 is the volume at time 0.
"""
function load_demands(filename::String)
    data = JSON3.read(read(filename, String))

    hasproperty(data, :num_time_slots) || error("Missing num_time_slots field")
    hasproperty(data, :demands)        || error("Missing demands field")

    # We only care about time slot 0 for Steps 1-5
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

# ------------------------------------------------------------------------------
# 5. Scenario loading
# ------------------------------------------------------------------------------

"""
    load_scenario(filename::String) -> NamedTuple{(:max_segments,), Tuple{Int}}

Loads the scenario file and returns only the max_segments value.
(Ignores budget and interventions for Steps 1-5.)
"""
function load_scenario(filename::String)
    data = JSON3.read(read(filename, String))

    hasproperty(data, :max_segments) || error("Missing max_segments field")

    return (max_segments = Int(data.max_segments),)
end

# ------------------------------------------------------------------------------
# 6. Convenience loader for a whole instance (time 0 only)
# ------------------------------------------------------------------------------

"""
    load_instance(prefix::String) -> (NetworkGraph, Vector, Int)

Loads the network, demands, and max_segments for a given instance prefix.
Example: load_instance("setA/setA-09")
"""
function load_instance(prefix::String)
    net = load_network("$prefix-net.json")
    demands = load_demands("$prefix-tm.json")
    scenario = load_scenario("$prefix-scenario.json")
    return net, demands, scenario.max_segments
end

end # module