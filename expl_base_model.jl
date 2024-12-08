b_model = Model(HiGHS.Optimizer)

@variable(b_model, 0 <= x[i in bids, l in locations] <= 1)
@variable(b_model, y[i in bids], Bin)
@variable(b_model, z[i in parent_bids], Bin)
@variable(b_model, 0 <= f[i in location_combinations, h in periods] <= cap)

@objective(b_model, Max, sum(data[i,"Price"]*data[i,"Quantity"]*x[i,a] for i in bids, a in locations) 
                        - sum(data[data.ParentBidID .== j,"FC"][1]*z[j] for j in parent_bids)
)

# nb: missing something to link x to f
market_balance_bm = Dict()
for a in locations 
    for h in periods
        df = data[(data.Location .== a) .& (data.Period .== h), :]
        b = df.BidID
        n = nrow(df)
        for l in locations
            if l != a
                @constraint(b_model, f[(a,l),h] == sum(df[i,"Quantity"]*x[b[i],l] for i in 1:n))
            end
        end
        market_balance_bm[a,h] = @constraint(b_model, sum(data[i,"Quantity"]*x[b[i],a] for i in 1:n) == 0)
    end
end


for i in bids
    @constraint(b_model, sum(x[i,a] for a in locations) <= y[i])    
    @constraint(b_model, sum(x[i,a] for a in locations) >= data[i,"AR"]*y[i])
end 

for j in parent_bids
    df = data[data.ParentBidID .== j,:]
    b = df.BidID
    n = nrow(df)
    @constraint(b_model, sum(y[b[i]] for i in 1:n) >= z[j])
    @constraint(b_model, sum(y[b[i]] for i in 1:n) <= n*z[j])
end

optimize!(b_model)
objective_value(b_model)

value(x[2,1])
sum([value(z[j]) for j in parent_bids])
value(f[(3,2),1])