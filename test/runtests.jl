include("../sortToys/sortToy.jl")

# As function of Integer width
#sortToy(40,2.^[10;],[40;],[Int8,Int16,Int32,Int64,Int128],"random")
# As function of N, with data sorted in reverse, Tim is best
sortToy(40,2.^[2:10;],[10;],[String],"reverse")
# As function of N, with data sorted in reverse, Radix is best 
sortToy(100,2.^[2:14;],[50;],[Float64],"random")
