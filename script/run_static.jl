#!/usr/bin/env julia
# run_static.jl – Runs the static model (time period 0) on all SetA instances
# with a 15-minute (900 seconds) time limit.

using JuMP, HiGHS, DataFrames, JSON3, Graphs

# Add the src folder to the load path
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

# Load our modules
using DataLoader, ECMP, StaticModel, Utils

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------

DATA_DIR = joinpath(@__DIR__, "..", "setA")
RESULTS_DIR = joinpath(@__DIR__, "..", "results")
RESULTS_FILE = joinpath(RESULTS_DIR, "static_results.csv")

# Ensure results directory exists
mkpath(RESULTS_DIR)

# List of instances to run (SetA-09 to SetA-14)
INSTANCES = ["setA-09", "setA-10", "setA-11", "setA-12", "setA-13", "setA-14"]

# Time limit in seconds (15 minutes)
TIME_LIMIT = 900.0

# ------------------------------------------------------------------------------
# Main experiment loop
# ------------------------------------------------------------------------------

function run_experiment()
    results = DataFrame(
        instance = String[],
        n_vertices = Int[],
        n_edges = Int[],
        n_demands = Int[],
        max_segments = Int[],
        status = String[],
        objective = Float64[],
        cpu_time = Float64[],
        gap = Float64[]
    )

    for instance in INSTANCES
        @info "Running instance: $instance"

        try
            # --- Load data ---
            net_file = joinpath(DATA_DIR, "$instance-net.json")
            tm_file  = joinpath(DATA_DIR, "$instance-tm.json")
            sc_file  = joinpath(DATA_DIR, "$instance-scenario.json")

            net = load_network(net_file)
            demands = load_demands(tm_file)
            scenario = load_scenario(sc_file)
            maxSeg = scenario.max_segments

            n_vertices = nv(net.graph)
            n_edges = ne(net.graph)
            n_demands = length(demands)

            @info "  - Vertices: $n_vertices, Edges: $n_edges, Demands: $n_demands, MaxSeg: $maxSeg"

            # --- Compute ECMP coefficients ---
            @info "  - Computing ECMP coefficients..."
            r_seg, r_arc = compute_r_all_pairs(net)

            n_active_arcs = length(r_arc)
            n_active_segments = length(keys(r_seg))
            @info "  - Active arcs: $n_active_arcs, Active segments: $n_active_segments"

            # --- Solve the model ---
            @info "  - Solving MILP (time limit: $TIME_LIMIT seconds)..."
            start_time = time()
            status, obj, solve_time, model = solve_static(
                net, demands, maxSeg, r_arc, TIME_LIMIT
            )
            total_time = time() - start_time

            # --- Extract gap if available ---
            best_bound = 0.0
            if status != MOI.OPTIMAL && status != MOI.ALMOST_OPTIMAL
                # Try to get the dual bound from the model
                if has_attribute(model, MOI.DualObjectiveValue())
                    best_bound = MOI.get(model, MOI.DualObjectiveValue())
                end
            end
            gap = Utils.compute_gap(status, obj, best_bound)

            # --- Record results ---
            push!(results, (
                instance = instance,
                n_vertices = n_vertices,
                n_edges = n_edges,
                n_demands = n_demands,
                max_segments = maxSeg,
                status = string(status),
                objective = obj,
                cpu_time = solve_time,
                gap = gap
            ))

            @info "  - Done! Status: $status, Objective: $obj, Time: $solve_time"

        catch e
            @error "Failed on instance $instance" exception=e
            # Record failure
            push!(results, (
                instance = instance,
                n_vertices = 0,
                n_edges = 0,
                n_demands = 0,
                max_segments = 0,
                status = "ERROR",
                objective = NaN,
                cpu_time = NaN,
                gap = NaN
            ))
        end
    end

    # --- Save results ---
    Utils.save_results(results, RESULTS_FILE)

    # Print summary
    @info "\n=== Summary ==="
    println(results)

    return results
end

# ------------------------------------------------------------------------------
# Run the experiment
# ------------------------------------------------------------------------------

if abspath(PROGRAM_FILE) == @__FILE__
    run_experiment()
end