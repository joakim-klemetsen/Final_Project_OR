# - Preliminary Information ----
# This file contains a solver for the welfare maximization problem using the IP Pricing approach.
# The solver requires the output of binary variables from the base model and that the base model is run to 
# obtain the packages, sets and parameters used in the model.  
# ------------------------------

# load base model output
base_output = CSV.read("output/base_model/base_model_extensive_output.csv", DataFrame)

# initialize ip pricing model
ip_model = Model(HiGHS.Optimizer)
set_silent(ip_model)

# variables
@variable(ip_model, 0 <= x[i in bids] <= 1)
@variable(ip_model, 0 <= f[i in location_combinations, j in periods] <= cap)
@variable(ip_model, y[i in bids])
@variable(ip_model, z[i in parent_bids])

# objective 
@objective(ip_model, Max, sum(data[i, "Quantity"] * data[i, "Price"] * x[i] for i in bids) - sum(FC[j] * z[j] for j in parent_bids))

# constraints

## - 1. market balance constraints ----
market_balance_ip = Dict()
for t in periods
    for l in locations
        filtered_data = data[(data.Period .== t) .& (data.Location .== l), :]
        b = filtered_data.BidID
        market_balance_ip[t, l] = @constraint(ip_model,
            sum(filtered_data[i, "Quantity"] * x[b[i]] for i in 1:nrow(filtered_data)) ==
            sum(f[(k, l), t] for k in locations if (k, l) in location_combinations) - 
            sum(f[(l, k), t] for k in locations if (l, k) in location_combinations)
        )
    end
end

## - 2. minimum acceptance ratio fulfillment ----

### link between acceptance ratio and binary variable
ar_link_cond = Dict()
for i in bids 
    ar_link_cond[i] = @constraint(ip_model, x[i] <= y[i])
end

### ensures that acceptance ratio atleast meets the min. acc. requirement
ar_geq_cond = Dict()
for i in bids
    ar_geq_cond[i] = @constraint(ip_model, x[i] >= (data[i,"AR"]+epsilon)*y[i])    
end

## - 3. fix binary variables ----

### ensures that y is fixed to the optimal value of y in the base model
fix_y = Dict()
for i in bids
    fix_y[i] = @constraint(ip_model, y[i] == base_output[i,"y_solution"])    
end

### ensures that z is fixed to the optimal value of z in the base model
fix_z = Dict()
for j in parent_bids
    fix_z[j] = @constraint(ip_model, z[j] == base_output[base_output.ParentBidID .== j,"z_solution"][1])    
end

# solve model
optimize!(ip_model)
objective_value(ip_model)

# output results
result_ip = copy(data)
result_ip.x_solution = [value(x[i]) for i in bids]
result_ip.y_solution = [value(y[i]) for i in bids]
result_ip.cleared_volume = [result_ip[i, "Quantity"] * result_ip[i, "x_solution"] for i in 1:nrow(result_ip)]
z_solution_map = Dict(j => value(z[j]) for j in parent_bids)
result_ip[!, "z_solution"] = [z_solution_map[result_ip[i, "ParentBidID"]] for i in 1:nrow(result_ip)]
dual_market_balance = Dict((t, l) => (-1) * dual(market_balance_ip[t, l]) for t in periods for l in locations)
result_ip[!, "pi_star"] = [
    dual_market_balance[(result_ip[i, "Period"], result_ip[i, "Location"])] for i in 1:nrow(result_ip)
]
dual_y_fix = Dict(i => (-1) * dual(fix_y[i]) for i in bids)
result_ip[!, "delta_star"] = [
    (result_ip[i, "y_solution"] == 0 ? 0 : dual_y_fix[result_ip[i, "BidID"]]) for i in 1:nrow(result_ip)
]

# Initialize surplus columns
result_ip[!, "consumer_surplus"] = zeros(nrow(result_ip))
result_ip[!, "producer_surplus"] = zeros(nrow(result_ip))

# Calculate consumer and producer surplus
for i in 1:nrow(result_ip)
    if result_ip[i, "cleared_volume"] >= 0
        result_ip[i, "consumer_surplus"] = ((result_ip[i, "Price"] - result_ip[i, "pi_star"]) * result_ip[i, "Quantity"] * result_ip[i, "x_solution"]) - (result_ip[i, "delta_star"] * result_ip[i, "y_solution"])
    else
        result_ip[i, "producer_surplus"] = ((result_ip[i, "Price"] - result_ip[i, "pi_star"]) * result_ip[i, "Quantity"] * result_ip[i, "x_solution"]) - (result_ip[i, "delta_star"] * result_ip[i, "y_solution"])
    end
end

# Write to CSV
CSV.write("output/ip_model/ip_model_output.csv", result_ip)

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