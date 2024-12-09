b_model = Model(HiGHS.Optimizer)

@variable(b_model, 0 <= x[i in bids, l in locations] <= 1)
@variable(b_model, y[i in bids], Bin)
@variable(b_model, z[i in parent_bids], Bin)

@objective(b_model, Max, sum(data[i,"Price"]*data[i,"Quantity"]*x[i,a] for i in bids, a in locations) 
                        - sum(data[data.ParentBidID .== j,"FC"][1]*z[j] for j in parent_bids)
)

market_balance_bm = Dict()
flow = Dict()
for a in locations 
    for h in periods
        df = data[(data.Location .== a) .& (data.Period .== h), :]
        b = df.BidID
        n = nrow(df)
        t = unique(data[data.Location .!= a,"Location"])
        market_balance_bm[a,h] = @constraint(b_model, sum(df[i,"Quantity"]*x[b[i],a] for i in 1:n) == 0)
        for loc in t 
            flow[a,h,loc] = @constraint(b_model, sum(df[i,"Quantity"]*x[b[i],loc] for i in 1:n) <= cap)
        end
    end
end

for i in bids
    @constraint(b_model, sum(x[i,a] for a in locations) <= y[i])    
    @constraint(b_model, sum(x[i,a] for a in locations) >= (data[i,"AR"]+epsilon)*y[i])
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

output_data = copy(data)
for loc in locations
    output_data[!, Symbol("x[i,$(loc)]")] = [value(x[i, loc]) for i in output_data.BidID]
end
output_data."y[i]" = [value(y[i]) for i in bids]
z_solution_map = Dict(j => value(z[j]) for j in parent_bids)
output_data[!, "z[j]"] = [z_solution_map[output_data[i, "ParentBidID"]] for i in 1:nrow(output_data)]
for loc in locations
    output_data[!, Symbol("cleared_volume[i,$(loc)]")] = [data[i,"Quantity"]*value(x[i, loc]) for i in output_data.BidID]
end
CSV.write("output/interim_bm_output.csv", output_data)