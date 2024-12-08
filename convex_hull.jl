# - Preliminary Information ----
# This file contains a solver for the welfare maximization problem using the Convex Hull approach.
# The solver requires the base model is run to 
# obtain the packages, sets and parameters used in the model.  
# ------------------------------

# initialize model
ch_model = Model(HiGHS.Optimizer)
set_silent(ch_model)

# variables
@variable(ch_model, 0 <= x[i in bids] <= 1)          
@variable(ch_model, 0 <= f[i in location_combinations, j in periods] <= cap)  
@variable(ch_model, alpha[i in bids])      
@variable(ch_model, beta[i in parent_bids])             

# objective
@objective(ch_model, Max, 
    sum(data[i, "Quantity"] * data[i, "Price"] * x[i] for i in bids) -
    sum(FC[j] * beta[j] for j in parent_bids)
)

# constraints

## - 1. market balance constraints ----
market_balance_ch = Dict()
for t in periods
    for l in locations
        filtered_data = data[(data.Period .== t) .& (data.Location .== l), :]
        b = filtered_data.BidID
        market_balance_ch[t, l] = @constraint(ch_model,
            sum(filtered_data[i, "Quantity"] * x[b[i]] for i in 1:nrow(filtered_data)) ==
            sum(f[(k, l), t] for k in locations if (k, l) in location_combinations) -
            sum(f[(l, k), t] for k in locations if (l, k) in location_combinations)
        )
    end
end

## - 2. convex hull constraints ----
ch_constraints = Dict() 
for i in bids
    @constraint(ch_model, x[i] <= alpha[i])
    @constraint(ch_model, (data[i,"AR"]+epsilon)*alpha[i] <= x[i])
    ch_constraints[i] = @constraint(ch_model, alpha[i] <= beta[data[data.BidID .== i,"ParentBidID"][1]])
end

# solve model
optimize!(ch_model)
objective_value(ch_model)

# output results
test = copy(data)
test.x_solution = [value(x[i]) for i in bids]
test.y_solution = [value(alpha[i]) for i in bids]
test.cleared_volume = [test[i,"Quantity"]*test[i,"x_solution"] for i in 1:nrow(test)]
z_solution_map = Dict(j => value(beta[j]) for j in parent_bids)
test[!, "z_solution"] = [z_solution_map[test[i, "ParentBidID"]] for i in 1:nrow(test)]
#dual_market_balance = Dict((t, l) => (-1)*dual(market_balance_ch[t, l]) for t in periods for l in locations)
#test[!, "pi_star"] = [
#    dual_market_balance[(test[i, "Period"], test[i, "Location"])] for i in 1:nrow(test)
#]
dual_y_fix = Dict(i => (-1)*dual(ch_constraints[i]) for i in bids)
test[!, "delta_star"] = [dual_y_fix[test[i, "BidID"]] for i in 1:nrow(test)]
#test.payment =[test[i,"pi_star"]*abs(test[i,"Quantity"])*test[i,"x_solution"]-test[i,"delta_star"]*test[i,"y_solution"] for i in 1:nrow(test)]
#test.surplus = [abs(test[i,"Price"]-test[i,"pi_star"])*abs(test[i,"Quantity"])*test[i,"x_solution"] for i in 1:nrow(test)]
CSV.write("output/ch_model/ch_model_interim_output.csv", test)

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
    result[index, :Price] = (-1) * dual(market_balance_ch[period, loc])
    
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
CSV.write("output/ch_model/result_flow_output_ch.csv", result)
