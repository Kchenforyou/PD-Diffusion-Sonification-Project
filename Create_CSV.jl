using DataFrames, CSV

#include("ScaffSpind.jl")
#include("Skewer.jl")

#mySST = ScaffSpindType(1000)
#spndl = size(mySST.N)[1]

#Clr = RainbowTypes(spndl)
#Sp = SkewerProcess(mySST, Clr)

function create_CSV(Sp::Vector{Matrix{Float64}}, f::Int64)
    frames = size(Sp)[1]
    num_space = 500
    digits = floor(Int, log(10, num_space)) + 1
    
    for i in 1:num_space # for now, we're doing every every 10 images (up to 100)
        df = DataFrame(Sp[i*2], :auto)
	label = floor(Int, log(10, i))
	n = digits - label
        CSV.write("CSV_files" * "$f" * "/SkewerProcess" * '0'^n * "$i.csv", df)
    end
    return nothing
end
