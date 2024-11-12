# loading packages
using CSV, DataFrames, JuMP, HiGHS

# loading data
data = CSV.read("data/ProjectData.csv", DataFrame)