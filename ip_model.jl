# - Preliminary Information ----
# This file contains a solver for the welfare maximization problem using the IP Pricing approach.
# The solver requires the output of binary variables from the base model and that the base model is run to 
# obtain the packages, sets and parameters used in the model.  
# ------------------------------

# load base model output
Y = CSV.read("output/base_model_y_output.csv", DataFrame)
Z = CSV.read("output/base_model_z_output.csv", DataFrame)

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

# market balance constraints
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
    ar_geq_cond[i] = @constraint(ip_model, x[i] >= data[i,"AR"]*y[i])    
end

# fix binary variables
fix_y = Dict()
for i in bids
    fix_y[i] = @constraint(ip_model, y[i] == Y[i,"Y"])    
end

fix_z = Dict()
for i in parent_bids
    fix_z[i] = @constraint(ip_model, z[i] == Z[i,"Z"])    
end

# solve model
optimize!(ip_model)

# output dual values for z
CSV.write("output/dual_fixed_z_ip.csv",DataFrame(ParentBidID = parent_bids, 
                                                 Z = [Z[i,"Z"] for i in parent_bids],
                                                 Dual = [dual(fix_z[i]) for i in parent_bids])
)

# Add solution values for x, y, and z to the original data
data_with_solutions = DataFrame(data)
data_with_solutions[!, "x_solution"] = [value(x[i]) for i in data.BidID]
data_with_solutions[!, "y_solution"] = [value(y[i]) for i in data.BidID]
CSV.write("output/ip_model_output.csv", data_with_solutions)
# Add z solutions (mapped by ParentBidID)
#z_solution_map = Dict(j => value(z[j]) for j in parent_bids)
#data_with_solutions[!, "z_solution"] = [z_solution_map[data[i, "ParentBidID"]] for i in 1:nrow(data)]

# Add dual values of the market_balance constraint
#dual_market_balance = Dict((t, l) => dual(market_balance[t, l]) for t in periods for l in locations)
#data_with_solutions[!, "dual_market_balance"] = [
#    dual_market_balance[(data[i, "Period"], data[i, "Location"])] for i in 1:nrow(data)
#]

# Write the updated data to a CSV file
#CSV.write("output/extended_base_model_output.csv", data_with_solutions)
dual(market_balance_ip[6,4])