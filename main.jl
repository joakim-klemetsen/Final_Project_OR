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
M = Dict()
for i in parent_bids 
   M[i] = nrow(data[data.ParentBidID .== i,:])
end
FC = Dict()
for i in parent_bids
    FC[i] = data[(data.ParentBidID .== i),"FC"][1]
end

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
    ar_geq_cond[i] = @constraint(m, x[i] >= data[i,"AR"]*y[i])    
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

objective_value(m)
for i in bids
    println(value(x[i]))    
end
for i in bids
    println(value(y[i]))    
end
for i in parent_bids
    println(value(z[i]))
end

# explore results
result = copy(data)
result.FAR = [value(x[i]) for i in bids ]
CSV.write("output/test_results.csv", result)

DataFrame(
    loc_comb = location_combinations,
    flow = [value(f[i,2]) for i in location_combinations]
)

# extract binary variable outputs to be used in IP Pricing case
y_fixed = [value(y[i]) for i in bids]
z_fixed = [value(z[j]) for j in parent_bids]

test = copy(data)
test.x = [value(x[i]) for i in bids]
test.y_fixed = [value(y[i]) for i in bids]
test[!, :z_fixed] = [value(z[parent_bid]) for parent_bid in test[!, :ParentBidID]]
CSV.write("output/test_fixed_bin.csv", test)
