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
function iterative_search(candidates, fixed_parents)
    true_candidates = Int[]
    false_candidates = Int[]
    while !isempty(candidates)
        search = DataFrame(ParentBidID = Int[], ObjectiveValue = Float64[])
        for i in candidates
            dt, obj_val = solve_model(fixed_parents, i)  # assuming this function is defined elsewhere
            push!(search, (ParentBidID = i, ObjectiveValue = obj_val))
        end
        min_value, min_index = findmin(search[!, "ObjectiveValue"])
        row_index = findfirst(fixed_parents.ParentBidID .== search.ParentBidID[min_index])  # Finding the row index
        if !isnothing(row_index)
            fixed_parents[row_index, "z_value"] = 1  # Correctly updating the DataFrame
        end
        # check if integral
        check_df, check_obj = solve_model(fixed_parents)  # assuming this function is defined elsewhere
        push!(true_candidates, search.ParentBidID[min_index])
        if any(x -> 0 < x < 1, check_df.z_value)
            candidates = setdiff(candidates, search.ParentBidID[min_index])
        else
            push!(false_candidates, candidates)
            break
        end
    end
    return true_candidates, false_candidates
end    

## ---
# Define the basis for the iterative search
## ---

# solve the model with pure continues relaxation; provides upper bound of solution
df, relaxed_model = solve_model()
any(x -> x <= 0 || x >= 1, df.ParentBidID)
# provide a data frame with the predetermined fixed values; parents with binary values from the pure relaxed model
initial_fixed_parents = df[(df.z_value .== 0) .| (df.z_value .== 1),:]

# extract the fractionally accepted parents from the pure model
initial_candidates = setdiff(parent_bids, df[(df.z_value .== 0) .| (df.z_value .== 1),"ParentBidID"])

## ---
# Perform interative search
## ---
t_can, f_can = iterative_search(initial_candidates, initial_fixed_parents)

df1, relaxed_model1 = solve_model(df[(df.z_value .== 0) .| (df.z_value .== 1),:], initial_candidates)

df1[92,"z_value"] = 1

df2, relaxed_model2 = solve_model(df[(df.z_value .== 0) .| (df.z_value .== 1),:])
