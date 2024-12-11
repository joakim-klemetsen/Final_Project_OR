## ---
# Define functions that can be used to solve the European Approach
### ---

# Function that solves the welfare maximization problem with the EU approach
function solve_model(predetermined=nothing, subject=nothing)

    # Initialize model
    model = Model(HiGHS.Optimizer)
    set_silent(model)

    # Define variables
    @variable(model, 0 <= x[i in bids] <= 1)
    @variable(model, 0 <= f[i in location_combinations, j in periods] <= cap)
    @variable(model, 0 <= y[i in bids] <= 1)
    @variable(model, 0 <= z[i in parent_bids] <= 1)

    # Define objective
    @objective(model, Max, sum(data[i, "Quantity"] * data[i, "Price"] * x[i] for i in bids) - sum(FC[j] * z[j] for j in parent_bids))

    # Define constraints

    ## Market balance
    for t in periods
        for l in locations
            filtered_data = data[(data.Period .== t) .& (data.Location .== l), :]
            b = filtered_data.BidID
            @constraint(model,
                sum(filtered_data[i, "Quantity"] * x[b[i]] for i in 1:nrow(filtered_data)) ==
                sum(f[(k, l), t] for k in locations if (k, l) in location_combinations) - 
                sum(f[(l, k), t] for k in locations if (l, k) in location_combinations)
            )
        end
    end

    ## Minimum acceptance ratio fulfillment
    for i in bids 
        @constraint(model, x[i] <= y[i])
    end
    for i in bids
        @constraint(model, x[i] >= (data[i, "AR"] + epsilon) * y[i])    
    end

    ## Fixed cost link
    for l in parent_bids
        filtered_data = data[(data.ParentBidID .== l), :]
        b = filtered_data.BidID
        @constraint(model, sum(y[b[i]] for i in 1:nrow(filtered_data)) <= M[l] * z[l])
        @constraint(model, sum(y[b[i]] for i in 1:nrow(filtered_data)) >= z[l])
    end 

    # Apply predetermined fixes if provided
    if predetermined !== nothing
        for i in 1:nrow(predetermined)
            @constraint(model, z[predetermined.ParentBidID[i]] == predetermined.z_value[i])
        end
    end

    # Apply subject fix if provided
    if subject !== nothing
        for i in subject 
            @constraint(model, z[i] == 0)
        end
    end

    # Solve model
    optimize!(model)

    # Extract the solution of the parent bid ids
    df = DataFrame(ParentBidID = parent_bids, z_value = [value(z[i]) for i in parent_bids])
    
    # Return welfare
    return df, objective_value(model)
end

# Function that performs an iterative search
function iterative_search(candidates)
    it_count = 1
    removed_candidates = Int[]
    fixed_parents = DataFrame(ParentBidID = Int[], z_value = Float64[])

    while !isempty(candidates)
        println("Outer loop iteration: ", it_count)
        search = DataFrame(ParentBidID = Int[], ObjectiveValue = Float64[])

        for i in candidates
            println("Removing parent: ", i)
            dt, obj_val = solve_model(fixed_parents, i) 
            push!(search, (ParentBidID = i, ObjectiveValue = obj_val))
        end

        if isempty(search)
            println("No more candidates to evaluate.")
            break
        end

        max_value, max_pos = findmax(search.ObjectiveValue)
        parent_to_remove = search.ParentBidID[max_pos]
        push!(removed_candidates, parent_to_remove)

        if !isnothing(parent_to_remove)
            push!(fixed_parents, (ParentBidID = parent_to_remove, z_value = 0))
            candidates = setdiff(candidates, [parent_to_remove])
        end

        check_df, check_obj = solve_model(fixed_parents)
        if any(x -> 0 < x < 1, check_df.z_value)
            println("No integral solution found by removing: ", removed_candidates)
            println("Continuing search!")
            it_count += 1
        else
            println("Integral solution found on iteration: ", it_count)
            break
        end
    end

    return check_df
end
    

## ---
# Define the basis for the iterative search
## ---

# solve the model with pure continues relaxation; provides upper bound of solution
df, relaxed_model = solve_model()
CSV.write("output/eu_model/pure_unrestricted.csv",df_test)

# provide a data frame with the predetermined fixed values; parents with binary values from the pure relaxed model
initial_fixed_parents = df[(df.z_value .== 0) .| (df.z_value .== 1),:]

# extract the fractionally accepted parents from the pure model
initial_candidates = setdiff(parent_bids, df[(df.z_value .== 0) .| (df.z_value .== 1),"ParentBidID"])

## ---
# Perform interative search
## ---
z_values_search = iterative_search(initial_candidates)
