# loading packages
using CSV, DataFrames, JuMP, HiGHS

# loading data
data = CSV.read("data/ModifiedProjectData.csv", DataFrame)

# initilize model
m = Model(HiGHS.Optimizer)

# sets
bids = 1:nrow(data)
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

nrow(data[data.ParentBidID .== 1,:])

# variables
@variable(m, 0 <= x[i in bids] <= 1)
@variable(m, 0 <= f[i in pipelines] <= cap)
@variable(m, y[i in bids], Bin)
@variable(m, z[i in parent_bids], Bin)

# objective
@objective(m, Max, sum(data[i,"Quantity"]*data[i,"Price"]*x[i] for i in bids) - sum(data[j,"FC"]*z[j] for j in parent_bids))

# constraints

## - 1. market balance ----
market_balance = Dict() 
for t in periods
    for l in locations
        filtered_data = data[data.Period .== t .& data.Location .== l, :]
        market_balance[l] = @constraint(m,
            sum(row.Quantity * x[row.BidID] for row in eachrow(filtered_data)) ==
            sum(f[(k, l)] for k in locations if (k, l) in location_combinations) - 
            sum(f[(l, k)] for k in locations if (l, k) in location_combinations)
        )
    end
end

## - 2. minimum acceptance ratio fulfillment ----

### link between acceptance ratio and binary variable
for i in bids 
    @constraint(m, x[i] <= y[i])
end

### ensures that acceptance ratio atleast meets the min. acc. requirement
for i in bids
    @constraint(m, x[i] >= data[i,"AR"]*y[i])    
end

## - 3. fixed cost link ----

### upper bound: becomes active when \sum(y_i) > 0 and forces z_l = 1
for l in locations
    @constraint(m, sum(y[i] for i in bids) <= M[l]*z[l])
end

### lower bound: becomes active when \sum(y_i) = 0 and forces z_l = 0
for l in locations
    @constraint(m, sum(y[i] for i in bids) >= z[l])
end


