# loading packages
using CSV, DataFrames, JuMP, HiGHS

# loading data
data = CSV.read("data/ProjectData.csv", DataFrame)



# initilize model
m = Model(HiGHS.Optimizer)

# sets
bids = 1:nrow(data)
parent_bids = unique(data.ParentBidID)
locations = unique(data.Location)
pipelines = [(l, k) for l in locations, k in locations if l != k]

# parameters
cap = 300

# variables
@variable(m, 0 <= x[i in bids] <= 1)
@variable(m, 0 <= f[i in pipelines] <= cap)
@variable(m, y[i in bids], Bin)
@variable(m, z[i in parent_bids], Bin)

# objective
@objective(m, Max, sum(data[i,"Quantity"]*data[i,"Price"]*x[i] for i in bids) - sum(data[j,"FC"]*z[j] for j in parent_bids))

