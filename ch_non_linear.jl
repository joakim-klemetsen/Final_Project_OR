using Ipopt

ch_test = Model(Ipopt.Optimizer)
set_silent(ch_test)

@variable(ch_test, 0 <= x[i in bids] <= 1)
@variable(ch_test, 0 <= f[i in location_combinations, j in periods] <= cap)
@variable(ch_test, 0 <= alpha[i in bids] <= 1)
@variable(ch_test, 0 <= beta[j in parent_bids] <= 1)

@objective(ch_test, Max, sum(data[i,"Quantity"]*data[i,"Price"]*x[i] for i in bids) 
                        - sum(FC[j]*beta[j] for j in parent_bids)
)

market_balance_ch_test = Dict()
for t in periods
    for l in locations
        filtered_data = data[(data.Period .== t) .& (data.Location .== l), :]
        b = filtered_data.BidID
        market_balance_ch_test[t, l] = @constraint(ch_test,
            sum(filtered_data[i, "Quantity"] * x[b[i]] for i in 1:nrow(filtered_data)) ==
            sum(f[(k, l), t] for k in locations if (k, l) in location_combinations) - 
            sum(f[(l, k), t] for k in locations if (l, k) in location_combinations)
        )
    end
end

ch_constraint = Dict()
for l in parent_bids
    filtered_data = data[(data.ParentBidID .== l),:]
    b = filtered_data.BidID
    n = nrow(filtered_data)
    ch_constraint[l] = @constraint(ch_test, sum(filtered_data[i,"AR"]*alpha[b[i]] for i in 1:n)*beta[l] <= sum(x[b[i]] for i in 1:n))
end
ch_constraint[1]
optimize!(ch_test)

dual(market_balance_ch_test[2,1])

#bids_dict = Dict()
#for i in bids 
#    bids_dict[i] = @constraint(ch_test, x[i] >= (data[i,"AR"])*alpha[i])
#end

#ch_constraint = Dict()
#for l in parent_bids
#    filtered_data = data[(data.ParentBidID .== l),:]
#    b = filtered_data.BidID
#    n = nrow(filtered_data)
#    ch_constraint[l] = @constraint(ch_test, sum(alpha[b[i]] for i in 1:n)*beta[l] <= sum(x[b[i]] for i in 1:n))
#end
#parent_link_dict = Dict()
#for l in parent_bids
#    filtered_data = data[(data.ParentBidID .== l),:]
#    b = filtered_data.BidID
#    parent_link_dict[l] = @constraint(ch_test, sum(x[b[i]] for i in 1:nrow(filtered_data)) <= beta[l])
#end
