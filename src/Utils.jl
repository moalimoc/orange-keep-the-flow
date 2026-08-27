module Utils

using DataFrames, CSV

"""
    compute_gap(status, objective_value, best_bound)

Computes the optimality gap for a solved model.
- If status is OPTIMAL, gap = 0.0
- Otherwise, gap = 1.0 - (best_bound / objective_value)
"""
function compute_gap(status, objective_value, best_bound)
    if status == MOI.OPTIMAL
        return 0.0
    elseif objective_value > 0 && best_bound > 0
        return 1.0 - (best_bound / objective_value)
    else
        return NaN  # Gap undefined
    end
end

"""
    save_results(df::DataFrame, filename::String)

Saves a DataFrame to a CSV file.
"""
function save_results(df::DataFrame, filename::String)
    CSV.write(filename, df)
    @info "Results saved to $filename"
end

end # module