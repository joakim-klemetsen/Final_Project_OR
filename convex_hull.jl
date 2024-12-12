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
    ch_constraints[i] =@constraint(ch_model, (data[i,"AR"]+epsilon)*alpha[i] <= x[i])
    @constraint(ch_model, alpha[i] <= beta[data[data.BidID .== i,"ParentBidID"][1]])
end

# solve model
optimize!(ch_model)
objective_value(ch_model)

# output results
result_ch = copy(data)
result_ch.x_solution = [value(x[i]) for i in bids]
result_ch.y_solution = [value(alpha[i]) for i in bids]
result_ch.cleared_volume = [result_ch[i, "Quantity"] * result_ch[i, "x_solution"] for i in 1:nrow(result_ch)]
z_solution_map = Dict(j => value(beta[j]) for j in parent_bids)
result_ch[!, "z_solution"] = [z_solution_map[result_ch[i, "ParentBidID"]] for i in 1:nrow(result_ch)]
dual_market_balance = Dict((t, l) => (-1) * dual(market_balance_ch[t, l]) for t in periods for l in locations)
result_ch[!, "pi_star"] = [
    dual_market_balance[(result_ch[i, "Period"], result_ch[i, "Location"])] for i in 1:nrow(result_ch)
]
dual_y_fix = Dict(i => (-1) * dual(ch_constraints[i]) for i in bids)
result_ch[!, "delta_star"] = [
    (result_ch[i, "y_solution"] == 0 ? 0 : dual_y_fix[result_ch[i, "BidID"]]) for i in 1:nrow(result_ch)
]

# Initialize surplus columns
result_ch[!, "consumer_surplus"] = zeros(nrow(result_ch))
result_ip[!, "producer_surplus"] = zeros(nrow(result_ip))

# Calculate consumer and producer surplus
for i in 1:nrow(result_ch)
    if result_ch[i, "cleared_volume"] >= 0
        result_ch[i, "consumer_surplus"] = ((result_ch[i, "Price"] - result_ch[i, "pi_star"]) * result_ch[i, "Quantity"] * result_ch[i, "x_solution"]) - (result_ch[i, "delta_star"] * result_ch[i, "y_solution"])
    else
        result_ch[i, "producer_surplus"] = ((result_ch[i, "Price"] - result_ch[i, "pi_star"]) * result_ch[i, "Quantity"] * result_ch[i, "x_solution"]) - (result_ch[i, "delta_star"] * result_ch[i, "y_solution"])
    end
end

# Write to CSV
CSV.write("output/ch_model/ch_model_output.csv", result_ch)

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
                   transfer_surplus = zeros(Float64, length(periods)*length(locations))
)
index = 1
for (loc, period) in [(l, p) for l in locations, p in periods]
    filtered_data = data[(data.Location .== loc) .& (data.Period .== period), :]
    result[index, :Location] = loc
    result[index, :Period] = period
    result[index, :Price] = (-1) * dual(market_balance_ip[period, loc])
    
    positive_quantity_bids = filtered_data[filtered_data.Quantity .>= 0, :]
    result[index, :Demand] = sum([positive_quantity_bids[i, "Quantity"] * value(x[positive_quantity_bids[i, "BidID"]]) for i in 1:nrow(positive_quantity_bids)])
    
    negative_quantity_bids = filtered_data[filtered_data.Quantity .< 0, :]
    result[index, :Supply] = sum([negative_quantity_bids[i, "Quantity"] * value(x[negative_quantity_bids[i, "BidID"]]) for i in 1:nrow(negative_quantity_bids)])
    
    # compute transfer surplus
    transfer_surplus = 0  
    for to_loc in 1:4  
        if loc != to_loc
            flow_volume = value(f[(loc, to_loc), period])
            price_difference = result[index, :Price] - (-1) * dual(market_balance_ip[period, to_loc])
            result[index, Symbol("f_to_$(to_loc)")] = flow_volume
            transfer_surplus += abs(price_difference) * flow_volume
        else
            result[index, Symbol("f_to_$(to_loc)")] = 0
        end
    end
    result[index, :transfer_surplus] = transfer_surplus
    index += 1
end
result.Netflow = [((result[i,"Demand"] + result[i,"Supply"])) for i in 1:nrow(result)]
CSV.write("output/ip_model/result_flow_output_ip.csv", result)

# auxiliary computations
## compute fixed cost
fixed_cost_df = unique(select(result_ip, :ParentBidID, :z_solution, :FC), :ParentBidID)
total_fixed_cost = sum(filter(row -> row.z_solution == 1, fixed_cost_df).FC)

# output ip model aggregated results
report_results = DataFrame(consumer_surplus = sum(result_ip.consumer_surplus),
                           producer_surplus = sum(sum(result_ip.producer_surplus)),
                           transfer_surplus = sum(result.transfer_surplus),
                           fixed_cost = total_fixed_cost,
                           compensation = (-1)*sum(result_ip.delta_star)                                                                                
)
CSV.write("output/ip_model/report_ip_aggregated_results.csv",report_results)