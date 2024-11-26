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