# - Preliminary Information ----
# This file contains a solver for the welfare maximization problem using the European approach.
# The solver requires the output of binary variables from the base model and that the base model is run to 
# obtain the packages, sets and parameters used in the model.  
# ------------------------------


# Initialize the model
eu_model = Model(HiGHS.Optimizer)
set_silent(eu_model)

# Define variables
@variable(eu_model, 0 <= x[i in bids] <= 1)
@variable(eu_model, 0 <= f[i in location_combinations, j in periods] <= cap)
@variable(eu_model, 0 <= y[i in bids] <= 1)
@variable(eu_model, 0 <= z[i in parent_bids] <= 1)

# Define objective
@objective(eu_model, Max, sum(data[i, "Quantity"] * data[i, "Price"] * x[i] for i in bids) - sum(FC[j] * z[j] for j in parent_bids))

# Define constraints

## Market balance
market_balance_eu = Dict()
for t in periods
    for l in locations
        filtered_data = data[(data.Period .== t) .& (data.Location .== l), :]
        b = filtered_data.BidID
        market_balance_eu[t, l] = @constraint(eu_model,
            sum(filtered_data[i, "Quantity"] * x[b[i]] for i in 1:nrow(filtered_data)) ==
            sum(f[(k, l), t] for k in locations if (k, l) in location_combinations) - 
            sum(f[(l, k), t] for k in locations if (l, k) in location_combinations)
        )
    end
end

## Minimum acceptance ratio fulfillment
ar_link_cond = Dict()
for i in bids 
    ar_link_cond[i] = @constraint(eu_model, x[i] <= y[i])
end

ar_geq_cond = Dict()
for i in bids
    ar_geq_cond[i] = @constraint(eu_model, x[i] >= (data[i, "AR"] + epsilon) * y[i])    
end

## Fixed cost link
fc_upper = Dict()
fc_lower = Dict()
for l in parent_bids
    filtered_data = data[(data.ParentBidID .== l), :]
    b = filtered_data.BidID
    fc_upper[l] = @constraint(eu_model, sum(y[b[i]] for i in 1:nrow(filtered_data)) <= M[l] * z[l])
    fc_lower[l] = @constraint(eu_model, sum(y[b[i]] for i in 1:nrow(filtered_data)) >= z[l])
end

# solve model
optimize!(eu_model)
objective_value(eu_model)

# output results
test = copy(data)
test.x_solution = [value(x[i]) for i in bids]
test.y_solution = [value(y[i]) for i in bids]
z_solution_map = Dict(j => value(z[j]) for j in parent_bids)
test[!, "z_solution"] = [z_solution_map[test[i, "ParentBidID"]] for i in 1:nrow(test)]
test.cleared_volume = [test[i,"Quantity"]*test[i,"x_solution"] for i in 1:nrow(test)]
CSV.write("output/eu_model/eu_model_interim_output.csv",test)

# output flows
result = DataFrame(Location = zeros(Int, length(periods)*length(locations)),
                   Period = zeros(Int, length(periods)*length(locations)),
                   Price = zeros(Float64, length(periods)*length(locations)),
                   Demand = zeros(Float64, length(periods)*length(locations)),
                   Supply = zeros(Float64, length(periods)*length(locations)),
                   Netflow = zeros(Float64, length(periods)*length(locations)),
                   f_to_1 = zeros(Float64, length(periods)*length(locations)),
                   f_to_2 = zeros(Float64, length(periods)*length(locations)),
                   f_to_3 = zeros(Float64, length(periods)*length(locations)),
                   f_to_4 = zeros(Float64, length(periods)*length(locations)),
)

index = 1
for (loc, period) in [(l, p) for l in locations, p in periods]
    # Filter data for the current location and period
    filtered_data = data[(data.Location .== loc) .& (data.Period .== period), :]
    
    # Fill in the result DataFrame
    result[index, :Location] = loc
    result[index, :Period] = period
    result[index, :Price] = (-1) * dual(market_balance[period, loc])
    
    # Calculate Demand
    positive_quantity_bids = filtered_data[filtered_data.Quantity .>= 0, :]
    result[index, :Demand] = sum(
        [positive_quantity_bids[i, "Quantity"] * value(x[positive_quantity_bids[i, "BidID"]]) 
         for i in 1:nrow(positive_quantity_bids)]
    )
    
    # Calculate Supply
    negative_quantity_bids = filtered_data[filtered_data.Quantity .< 0, :]
    result[index, :Supply] = sum(
        [negative_quantity_bids[i, "Quantity"] * value(x[negative_quantity_bids[i, "BidID"]]) 
         for i in 1:nrow(negative_quantity_bids)]
    ) 

    # Fill in flows
        # Fill in the flow columns
        for to_loc in 1:4  # Assuming 4 locations numbered 1 to 4
            if loc != to_loc
                # Flow is defined
                result[index, Symbol("f_to_$(to_loc)")] = value(f[(loc, to_loc), period])
            else
                # Flow is undefined (same location)
                result[index, Symbol("f_to_$(to_loc)")] = 0
            end
        end
        
        index += 1    
end
result.Netflow = [((result[i,"Demand"] + result[i,"Supply"])) for i in 1:nrow(result)]
CSV.write("output/base_model/result_flow_output_bm.csv", result)
