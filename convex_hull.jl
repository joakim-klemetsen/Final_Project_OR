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
@variable(ch_model, 0 <= alpha[i in bids] <= 1)      
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
convex_hull_constraints = Dict()
for i in bids
    parent = data[data.BidID .== i, "ParentBidID"][1]

    ### ensures that alpha sets an upper limit to x for each i 
    convex_hull_constraints[i] = @constraint(ch_model, x[i] <= alpha[i])

    ### ensures that alpha is > 0 if and only if beta > 0 for each bid i
    @constraint(ch_model, alpha[i] <= beta[parent])
end

## - 3. fixed cost link constraints ----
fixed_cost_constraints = Dict()
for l in parent_bids
    filtered_data = data[data.ParentBidID .== l, :]
    b = filtered_data.BidID

    ### ensures that beta > 0 if and only if at least one child alpha is > 0  
    fixed_cost_constraints[l] = @constraint(ch_model,
        sum(alpha[b[i]] for i in 1:nrow(filtered_data)) >= beta[l]
    )
end

# solve model
optimize!(ch_model)
objective_value(ch_model)