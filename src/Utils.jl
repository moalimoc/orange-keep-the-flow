module Utils

using DataFrames, CSV, JuMP

export compute_gap, save_results

function compute_gap(status, objective_value, best_bound)
    if status == JuMP.OPTIMAL
        return 0.0
    elseif objective_value > 0 && best_bound > 0
        return 1.0 - (best_bound / objective_value)
    else
        return NaN
    end
end

function save_results(df::DataFrame, filename::String)
    CSV.write(filename, df)
    @info "Results saved to $filename"
end

end # module