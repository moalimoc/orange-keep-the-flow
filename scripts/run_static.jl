#!/usr/bin/env julia
# run_static.jl – Runs static model (time 0) on a single instance with demand limit.

using JuMP, HiGHS, DataFrames, CSV, JSON3, Graphs

# ------------------------------------------------------------------------------
# 1. Load the local modules by including the source files
# ------------------------------------------------------------------------------
src_dir = joinpath(@__DIR__, "..", "src")

include(joinpath(src_dir, "DataLoader.jl"))
include(joinpath(src_dir, "ECMP.jl"))
include(joinpath(src_dir, "StaticModel.jl"))
include(joinpath(src_dir, "Utils.jl"))

# Now the modules are defined in Main; bring them into scope
using .DataLoader, .ECMP, .StaticModel, .Utils

# ------------------------------------------------------------------------------
# SETTINGS – change these as needed
# ------------------------------------------------------------------------------
INSTANCE = "setA-01"          # which instance to run
MAX_DEMANDS = 200             # number of demands to use 
TIME_LIMIT = 900.0            # 15 minutes
# ------------------------------------------------------------------------------

function main()
    # ----------------------------------------------------------------------
    # 1. Load data
    # ----------------------------------------------------------------------
    @info "Loading instance: $INSTANCE"
    net = load_network(joinpath("setA", "$INSTANCE-net.json"))
    demands = load_demands(joinpath("setA", "$INSTANCE-tm.json"))
    scenario = load_scenario(joinpath("setA", "$INSTANCE-scenario.json"))
    maxSeg = scenario.max_segments

    # Limit demands
    if length(demands) > MAX_DEMANDS
        @info "  - Truncating demands to $MAX_DEMANDS (original: $(length(demands)))"
        demands = demands[1:MAX_DEMANDS]
    else
        @info "  - Using all $(length(demands)) demands"
    end

    n_vertices = nv(net.graph)
    n_edges = ne(net.graph)
    @info "  - Vertices: $n_vertices, Edges: $n_edges, MaxSeg: $maxSeg"

    # ----------------------------------------------------------------------
    # 2. Compute ECMP coefficients
    # ----------------------------------------------------------------------
    @info "  - Computing ECMP coefficients..."
    r_seg, r_arc = compute_r_all_pairs(net)
    @info "  - Active arcs: $(length(r_arc)), Active segments: $(length(keys(r_seg)))"

    # ----------------------------------------------------------------------
    # 3. Solve the model
    # ----------------------------------------------------------------------
    @info "  - Solving MILP (time limit: $TIME_LIMIT seconds)..."
    status, obj, solve_time, model = solve_static(net, demands, maxSeg, r_arc, TIME_LIMIT)

    # ----------------------------------------------------------------------
    # 4. Record results
    # ----------------------------------------------------------------------
    gap = NaN
    if status == JuMP.OPTIMAL
        gap = 0.0
    elseif model !== nothing && has_attribute(model, JuMP.ObjectiveBound())
        best_bound = JuMP.objective_bound(model)
        gap = Utils.compute_gap(status, obj, best_bound)
    end

    results = DataFrame(
        instance = [INSTANCE],
        n_vertices = [n_vertices],
        n_edges = [n_edges],
        n_demands = [length(demands)],
        max_segments = [maxSeg],
        status = [string(status)],
        objective = [obj],
        cpu_time = [solve_time],
        gap = [gap]
    )

    # ----------------------------------------------------------------------
    # 5. Save and print
    # ----------------------------------------------------------------------
    mkpath("results")
    CSV.write("results/static_results.csv", results)
    @info "Results saved to results/static_results.csv"

    println("\n=== Summary ===")
    println(results)

    return results
end

# ------------------------------------------------------------------------------
# Run the script
# ------------------------------------------------------------------------------
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end