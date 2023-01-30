function SkewerProcess(mySST::ScaffSpindType,
	Clr::Array{Float64,2}, # Colour is a 2 by n array
	epochs::Int64 = 500)  # how many frames is the animation?

	IPP = Vector{Array{Float64,2}}(undef,epochs) # IPP is a 2 by n array
	# "interval partition process"
	# an array of *epochs* interval partitions (compositions of n, really)
	for i in 1:epochs
		IPP[i] = Array{Float64}(undef,0,
			1 +  # block size (population) for each block (herd)
			3 + # size(Clr)[2] +  # 3(?) RGB color coordinates for each block (herd) alive at that time
			size(mySST.ptype[1])[2]  # 2(?) position (or type) coordinates
			)
	end
	# we make a 500-array vector, with each array being a 1x6 matrix in 2 dimensions.
	
	Y = Float64[] 	# at what levels do we evaluate the skewer?
	if mySST.theta==0  Y = collect(range(0,length=epochs,stop=maximum(mySST.X[:,3])))
	else  Y = collect(range(0,length=epochs,stop=mySST.X[1,2])) # KC: By indexing the 2nd column of the first row, we get the birth time of the left most spindle--the highest.
	end # KC: We subdivide from 0 to the highest birth time into 500 intervals. Between 0 to mySST.X[1,2], these are times where all the spindles will be born and "killed".
	
	dy = Y[2] # The '1st' interval.
	numSpdl = length(mySST.N) # number of spindles.
	L,K=1,1
	spHeights = [1]
	
	for i in 1:numSpdl
		# for each spindle, if it crosses at least one slice level, we
		# slice it at each level in Y and add the slices to the appropriate IPs
		L = 1 + ceil(Int,mySST.X[i,2]/dy) # round float to lowest integer; starting at the first spindle, L = 500. We have L intervals.
			# = lowest index where Y[L] >= mySST.X[i,2]
			# the "1+" is because Y[1] = 0 instead of = dy
		K = min(1 + floor(Int,mySST.X[i,3]/dy) , epochs ) - L # the argument in floor() finds the number of intervals between the highest death time.
			# = highest index where Y[L+K] <= mySST.X[i,3]
			# = 1 less than number of slices where this spindle contributes to the skewer
		if K >= 0
			spHeights = round.(Int, Y[L:(L+K)] .- mySST.X[i,2] ) .+ 1
			for j in 0:K
				IPP[L+j] = vcat(IPP[L+j] ,  # add a new herd to the list of herds alive at time Y[L+j]
						hcat([mySST.N[i][spHeights[j+1]]], # population of new herd
							Clr[i,1:3]' , # color of herd KC: Apostrophe is for transposing from column to row vector.
								# NF: I added in "1:3" for row indices extracted from Clr, since the RainbowTypes function
								# outputs an array with an extraneous, 4th column
							mySST.ptype[i][spHeights[j+1],:]'  # KC to add ptype of herd at that time; # KC: note the colon, we want the entire row of ptype.
							))
			end
		end
	end
	return IPP
end


function SkewerProcess2(mySST::ScaffSpindType,
	Clr::Array{Float64,2}) # Colour is a 2 by n array

	maxHeight = 0 # KC: This is what used to be epochs. This will be used for the dimension of IPP.
	if mySST.theta==0  maxHeight = maximum(mySST.X[:,3])
	else  maxHeight = mySST.X[1,2] # KC: By indexing the 2nd column of the first row, we get the birth time of the left most spindle--the highest.
	end # KC: We subdivide from 0 to the highest birth time into 500 intervals. Between 0 to mySST.X[1,2], these are times where all the spindles will be born and "killed".
	
	IPP = Vector{Array{Float64,2}}(undef,maxHeight) # IPP is a 2 by n array
	# "interval partition process"
	# an array of *epochs* interval partitions (compositions of n, really)
	for i in 1:maxHeight
		IPP[i] = Array{Float64}(undef,0,
			1 +  # block size (population) for each block (herd)
			3 + # size(Clr)[2] +  # 3(?) RGB color coordinates for each block (herd) alive at that time
			size(mySST.ptype[1])[2]  # 2(?) position (or type) coordinates
			)
	end
	# we make a vector of arrays, with each array being a 1x6 matrix in 2 dimensions.
	
	numSpdl = length(mySST.N) # number of spindles.
	m,M=1,1
	
	for i in 1:numSpdl
        m = mySST.X[i,2]
        M = min(mySST.X[i,3] , maxHeight) - 1 # shouldn't these by indented?
		if m < maxHeight
			for j in m:M
				#if j==maxHeight # KC: I added this because we keep indexing IPP[maxHeight + 1], which is out of range.
						# KC: program crashes during first iteration because IPP isn't large enough to be indexed by j+1
						# KC: Why is IPP < M?
				#	break
				#end
				IPP[j+1] = vcat(IPP[j+1] ,  # add a new herd to the list of herds alive at time j
						hcat([mySST.N[i][j-m+1]], # population of new herd
							Clr[i,1:3]' , # color of herd KC: Apostrophe is for transposing from column to row vector.
								# NF: I added in "1:3" for row indices extracted from Clr, since the RainbowTypes function
								# outputs an array with an extraneous, 4th column
							mySST.ptype[i][j-m+1,:]'  # KC: note the colon, we want the entire row of ptype.
							))
			end
		end
	end
	return IPP
end
