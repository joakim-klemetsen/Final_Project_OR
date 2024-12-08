# loading packages
using CSV, DataFrames, JuMP, HiGHS

# loading data
data = CSV.read("data/ModifiedProjectData.csv", DataFrame)

# initilize model
m = Model(HiGHS.Optimizer)
set_silent(m)

# sets
bids = unique(data.BidID)
parent_bids = unique(data.ParentBidID)
locations = unique(data.Location)
location_combinations = [(l, k) for l in locations, k in locations if l != k]
periods = unique(data.Period)

# parameters
cap = 300
M = Dict(parent_bids .=> [nrow(data[data.ParentBidID .== i, :]) for i in parent_bids])
FC = Dict(parent_bids .=> [data[data.ParentBidID .== i, "FC"][1] for i in parent_bids])
epsilon = 1e-5

# variables
@variable(m, 0 <= x[i in bids] <= 1)
@variable(m, 0 <= f[i in location_combinations, j in periods] <= cap)
@variable(m, y[i in bids], Bin)
@variable(m, z[i in parent_bids], Bin)

# objective
@objective(m, Max, sum(data[i,"Quantity"]*data[i,"Price"]*x[i] for i in bids) - sum(FC[j]*z[j] for j in parent_bids))

# constraints

## - 1. market balance ----
market_balance = Dict()
for t in periods
    for l in locations
        filtered_data = data[(data.Period .== t) .& (data.Location .== l), :]
        b = filtered_data.BidID
        market_balance[t, l] = @constraint(m,
            sum(filtered_data[i, "Quantity"] * x[b[i]] for i in 1:nrow(filtered_data)) ==
            sum(f[(k, l),t] for k in locations if (k, l) in location_combinations) - 
            sum(f[(l, k),t] for k in locations if (l, k) in location_combinations)
        )
    end
end

## - 2. minimum acceptance ratio fulfillment ----

### link between acceptance ratio and binary variable
ar_link_cond = Dict()
for i in bids 
    ar_link_cond[i] = @constraint(m, x[i] <= y[i])
end

### ensures that acceptance ratio atleast meets the min. acc. requirement
ar_geq_cond = Dict()
for i in bids
    ar_geq_cond[i] = @constraint(m, x[i] >= (data[i,"AR"]+epsilon)*y[i])    
end

## - 3. fixed cost link ----

### upper bound: becomes active when \sum(y_i) > 0 and forces z_l = 1
fc_upper = Dict()
for l in parent_bids
    filtered_data = data[(data.ParentBidID .== l),:]
    b = filtered_data.BidID
    fc_upper[l] = @constraint(m, sum(y[b[i]] for i in 1:nrow(filtered_data)) <= M[l]*z[l])
end

### lower bound: becomes active when \sum(y_i) = 0 and forces z_l = 0
fc_lower = Dict()
for l in parent_bids
    filtered_data = data[(data.ParentBidID .== l),:]
    b = filtered_data.BidID
    fc_lower[l] = @constraint(m, sum(y[b[i]] for i in 1:nrow(filtered_data)) >= z[l])
end

# solve model
optimize!(m)

# Output optimal binary variables that
y_fixed = DataFrame(BidID = data.BidID,
                    Y = [value(y[i]) for i in bids]
)

CSV.write("output/base_model/base_model_y_output.csv",y_fixed)

z_fixed = DataFrame(ParentBidID = parent_bids,
                    Z = [value(z[i]) for i in parent_bids]
)

CSV.write("output/base_model/base_model_z_output.csv",z_fixed)

CSV.write("output/base_model/base_model_output.csv", DataFrame(BidID = data.BidID,
                                                    X = [value(x[i]) for i in bids],
                                                    Y = [value(y[i]) for i in bids]))

# output results
test = copy(data)
test.x_solution = [value(x[i]) for i in bids]
test.y_solution = [value(y[i]) for i in bids]
z_solution_map = Dict(j => value(z[j]) for j in parent_bids)
test[!, "z_solution"] = [z_solution_map[test[i, "ParentBidID"]] for i in 1:nrow(test)]
test.cleared_volume = [test[i,"Quantity"]*test[i,"x_solution"] for i in 1:nrow(test)]
CSV.write("output/base_model/base_model_interim_output.csv",test)

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
