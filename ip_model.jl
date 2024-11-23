# initilize model
m_ip = Model(HiGHS.Optimizer)
set_silent(m_ip)

# variables
@variable(m_ip, 0 <= x[i in bids] <= 1)
@variable(m_ip, 0 <= f[i in location_combinations, j in periods] <= cap)
@variable(m_ip, y[i in bids], Bin)
@variable(m_ip, z[i in parent_bids], Bin)

# objective
@objective(m_ip, Max, sum(data[i,"Quantity"]*data[i,"Price"]*x[i] for i in bids) - sum(FC[j]*z[j] for j in parent_bids))

# constraints

## - 1. market balance ----
market_balance = Dict()
for t in periods
    for l in locations
        filtered_data = data[(data.Period .== t) .& (data.Location .== l), :]
        b = filtered_data.BidID
        market_balance[t, l] = @constraint(m_ip,
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
    ar_link_cond[i] = @constraint(m_ip, x[i] <= y[i])
end

### ensures that acceptance ratio atleast meets the min. acc. requirement
ar_geq_cond = Dict()
for i in bids
    ar_geq_cond[i] = @constraint(m_ip, x[i] >= data[i,"AR"]*y[i])    
end

## - 3. fixed cost link ----

### upper bound: becomes active when \sum(y_i) > 0 and forces z_l = 1
fc_upper = Dict()
for l in parent_bids
    filtered_data = data[(data.ParentBidID .== l),:]
    b = filtered_data.BidID
    fc_upper[l] = @constraint(m_ip, sum(y[b[i]] for i in 1:nrow(filtered_data)) <= M[l]*z[l])
end

### lower bound: becomes active when \sum(y_i) = 0 and forces z_l = 0
fc_lower = Dict()
for l in parent_bids
    filtered_data = data[(data.ParentBidID .== l),:]
    b = filtered_data.BidID
    fc_lower[l] = @constraint(m_ip, sum(y[b[i]] for i in 1:nrow(filtered_data)) >= z[l])
end

## - 4. fix binary variables

### fixing bids binary variables
fix_y = Dict()
for i in bids 
   fix_y[i] = @constraint(m_ip, y[i] == y_fixed[i]) 
end

### fixing parent id binary variable
fix_z = Dict()
for i in parent_bids 
   fix_y[i] = @constraint(m_ip, z[i] == z_fixed[i]) 
end

# solve model
optimize!(m_ip)

objective_value(m_ip)
dual(market_balance[1,1])
dual(fix_z[1])
fix_z[2]