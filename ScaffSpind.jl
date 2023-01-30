using Distributions, Random

alpha = .5
theta = .7
InitPop = 10

# As InitPop -> infty, spindle birth rate (1/(avg spindle length - 1)) 
# should decrease like alpha / InitPop. However, we get nicer pictures 
# faster w/ lower InitPop, so we run some code to check the average 
# spindle length. If we set SBR too high, our code has a chance to loop
# infinitely, so we err on the side of setting SBR too low.
# Here are some recorded values:

# SBR = spindle birth rate

#SBR = 1/35.2  # for a=.3, initpop=10
SBR = 1/21.06  # for a=.5, initpop=10; 1/21.02 is probably more accurate
#SBR = 1/17.7  # for a=.6, initpop=10
#SBR = 1/16.3  # for a=.7, initpop=10
#SBR = 1/13.55  # for a=.8, initpop=10


mutable struct ScaffSpind
	X::Array{Int64,2} # nx3 array
		# 1st col is x coord
		# 2nd col is y-coord at bottom of jump
		# 3rd col is y-coord at top of jump
	N::Vector{Vector{Int64}}
		# vec of spindles. Each spindle is a vec of amplitudes
	p::Vector{Int64}
		# index in N of the parent of each spindle;
		# = 0 if spindle has no parent in N
		# if the 10th spindle in N is a child of the 4th spindle in N, then
			# p[10] = 4 ("parent of spindle 10 is spindle 4")
	alpha::Float64
	theta::Float64
	function ScaffSpind(X::Array{Int64,2} = Array{Int64}(undef,0,3),
						N::Vector{Vector{Int64}} = Vector{Int64}[],
						p::Vector{Int64} = Int64[]
						)
		return new(X,N,p,alpha,theta)
	end
end

mutable struct ScaffSpindType
	X::Array{Int64,2} # nx3 array
		# 1st col is x coord
		# 2nd col is y-coord at bottom of jump
		# 3rd col is y-coord at top of jump
	N::Vector{Vector{Int64}}
		# vec of spindles. Each spindle is a vec of amplitudes
	p::Vector{Int64}
		# index in N of the parent of each spindle;
		# = 0 if spindle has no parent in N
	alpha::Float64
	theta::Float64
	#Clr::Array{Float64,2} # nx3 array
		# R-G-B for each spindle
	ptype::Vector{Array{Int64,2}} # genetic marker of each "bison herd" represented by varying distance
		# vec of Lx2 arrays; 1st col is x coords; 2nd is y-coords
		# where size(ptype[j]) = [length(N[j]) , 2]
		# (flexibility to do higher-dimensional positions, if wanted)
		# KC: I want to use rand(Uniform(-1, 1), 1, 2) to generate the initial position of the spindles. 
		# I'm not too clear as how to incorporate the output of this function here...
		# What's a type alias? The output of rand(Uniform(-1, 1), 1, 2)
		
		# KC: Should probably make this discrete (800.*rand(Uniform(-1, 1), 1, 2))) and then applying ceiling function.
		

	# color has been deleted
	function ScaffSpindType(X::Array{Int64,2} = Array{Int64}(undef,0,3),
						N::Array{Array{Int64,1},1} = Array{Int64,1}[],
						p::Vector{Int64} = Int64[],
						ptype::Vector{Array{Int64,2}} = Array{Int64,2}[]
						)
		return new(X,N,p,alpha,theta,ptype)
	end
end

function ScaffSpindType(SS::ScaffSpind,
	ptype::Vector{Array{Int64,2}} = Vector{Array{Int64,2}}[])
	return ScaffSpindType(SS.X,SS.N,SS.p,ptype)
end

mySeed = MersenneTwister(2603725848649352970)

function spindle(CM::Int64 = InitPop)
# approximation of BESQ(-2*alpha) process via a critical Galton-Watson
# branching process w/ variance 1 w/ emmigration rate alpha
# (normally, emmigration means deterministic, but we want fractional rate,
# so we do a Bernoulli(alpha)).
	x = [CM]
	while (x[end] > 0)
		push!(x, 2*rand(mySeed,Binomial(x[end],.50)) -
				rand(mySeed,Bernoulli(alpha)))
	end
	x[end] = 0
	return x
end


function spindleLengthTest(n::Int64)
	# estimates the mean and stddev of the length of a spindle in n trials;
	# depends on value of (global variable) alpha
	mySum,mySumSq,L = 0,0,0
	for i in 1:n
		L = length(spindle())
		mySum += L
		mySumSq += L^2
	end
	return [mySum/n, sqrt((mySumSq/n - (mySum/n)^2))]
end

# If you don't know what spindle birth rate (SBR) to use in your scaffolding based
# on the values of alpha and initpop, run the following line of code:

# SLT = spindleLengthTest(Int64(round(20000000/InitPop)))
# SBR = 1 / (SLT[1]-1 + 3*SLT[2])

# takes a while to run, but not forever...
# length of spindle should go linearly w/ initpop,
# so divide by initpop to control run time

function OCRP(minSize::Int64, massBD::Int64)
	totPop = 1 # actual population
	ePop = 0 # effective population, only includes large enough tables (>minSize)
	K = 1 # total number of tables
	IP = [1]  # list of tables (interval partition / set composition of [n])
	p = [1-alpha,theta,alpha]
	# vector of probability weights of being added to each table / in each new table slot
	# = vcat(IP-alpha , theta , rep(alpha,tables))

	J = 1 # declaring placeholder var

	while ePop < massBD
		J = rand(mySeed,Categorical( p ./ (totPop+theta) ))
		# choose a random number from 1 to 2K+1, w/ vector of probabilities
		# p / (totPop+theta) (dividing to normalize p, which is a vec of probab weights)
		if J <= K
			# add customer to Jth table
			IP[J] += 1
			p[J] += 1
			if IP[J] == minSize  		ePop += IP[J] # add table that crossed min pop threshold
			elseif IP[J] > minSize  	ePop += 1 # add new customer
			end
		else
			# add new table to the left of the (J-K)th (or on far right if J-K=K+1)
			splice!(IP,(J-K):(J-K-1),1) # add in between (J-K-1)st and (J-K)th tbl
			splice!(p,(J-K):(J-K-1),1-alpha) # update vector of insertion probabs for singleton table
			push!(p,alpha) # also, now there is one more insertion site
			K += 1 # update table count
		end
		totPop += 1 # update total population
	end
	return IP[IP .>= minSize] # delete small tables
end


function len(mySST::ScaffSpindType)
	if size(mySST.X)[1] == 0  return 0  end
	return mySST.X[end,1] + mySST.X[end,3]
end

function len_pm(mySST::ScaffSpindType)
	if size(mySST.X)[1] == 0  return 0  end
	return mySST.X[end,1] + mySST.X[end,3] - min(mySST.X[1,1],0)
end

function SSheight(SST::ScaffSpindType)  return maximum(SST.X[:,3])  end

function ScaffShift(mySST::ScaffSpindType, x::Int64 = 0, y::Int64 = 0)
	if (x != 0)  mySST.X[:,1] .+= x  end
	if (y != 0)
		mySST.X[:,2] .+= y
		mySST.X[:,3] .+= y
	end
	return mySST
end

function randwalk(L::Int64)
	# assume we're on an infinite plane... for now.
	#pos = Array{Int64}[]
	#init_c = floor.(rand(Uniform(-400, 400), 1, 2)) # note each herd spawn must start from the position of its parent.
	#push!(pos, init_c)

	# random walk goes here. This returns integers.
	step_x = 2 .* bitrand(mySeed, L) .- 1
	step_y = 2 .* bitrand(mySeed, L) .- 1
	
	x = cumsum(step_x)
	y = cumsum(step_y)
	
	return [x y]
end

function Clade(CM::Int64 = InitPop, InitType::Array{Int64,2} = [0 0]) # we need an optional argument for the range/size of the torus/grid.
	f = spindle(CM) # vector of spindles initialized here in Clade
	
	# randwalk will be called here, and accept 'f' as an argument.
	rw = randwalk(length(f)) .+ InitType
	DX = rand(mySeed,Geometric(SBR))+1  # +1 because in Julia, Geometric RVs can be 0, instead of starting at 1
	# SBR is a global variable.
	if DX > length(f)-1  return ScaffSpindType([0 0 (length(f)-1)],[f],[0],[rw])  end
	
	mySST = ScaffSpindType([0 0 (length(f)-1)],[f],[0],[rw]) # replace with MySST, color deleted.
	
	
	# KC: [f] is a nested list
	# KC: [0] means no parent.
	# KC: Modify this function such that we get a position process. 
	# KC: We'll need to pass a 1-D vector [pos] with an entry that is a 2d array of the spindle position at time 't'.
	# KC: delete colour for now.
	
	# KC: problem for Clade:
	# Assume there is a "Spindle G". When generating a sub-branch, we assume that it is a child of G. 
	# Based on G's position in time (the y-coordinate from attribute X), what is the position (relative to ptype)
	# of the herd when it spawns a child?
	
	# KC's possible solution. 
		# To find time, we'll need the height difference between y1 and y2 (where a new Clade gets generated).
		# Then round. Can't decide on floor vs ceiling
		# Because rate is 1 step on rw/1 unit of time, then use the rounded difference to index the random walk vectors.
	
	curPos = [DX, length(f)-1-DX] # track x and y position from X.
	subBranch = ScaffSpindType()
	while curPos[2] > 0
		# subBranch calls Clade to create children of children.
		BT = curPos[2] # BT is the birth time of our sub-branch
		subBranch = Clade(InitPop, [rw[BT, 1] rw[BT, 2]]) # BT will be used to find the randwalk coordinates @ the time the sub-branch is born.
		# KC: append on to parent array (p) the parent array of my sub-branch.
		# the addition below is to tell program that spindle+1 isn't root!
		append!(mySST.p, vcat([1], subBranch.p[2:end] .+ length(mySST.N)))
		append!(mySST.N,subBranch.N)
		append!(mySST.ptype, subBranch.ptype) 
		mySST.X = vcat(mySST.X,
			subBranch.X .+ [curPos[1] curPos[2] curPos[2]]) # this shifts position of new clade according to where parent is
		DX = rand(mySeed,Geometric(SBR))+1
		curPos += [len(subBranch) + DX , -DX]
	end
	return mySST 
		
end


function SSSlope(SST::ScaffSpindType)
	nSpdl = length(SST.N)
	return mean( (SST.X[1:(nSpdl-1),3] - SST.X[2:nSpdl,2]) ./ (SST.X[2:nSpdl,1] - SST.X[1:(nSpdl-1),1]) )
end


function InitTypes(sidelength::Int64 = 200)
	return ceil.(Int64,rand(Uniform(-1*sidelength/2, sidelength/2), 1, 2)) # random coordinates
end


function ScaffSpindType(Scale::Real)
	# KC: This below is not a recursive call. Rather we're function over-loading.
	# The ScaffSpindType() from the 'class' above outputs an 'empty' data structure.
	mySST = ScaffSpindType()
	if Scale < 1   return mySST   end

	CM = OCRP(InitPop,ceil(Int64,Scale*InitPop))
	curX = 0
	nextClade = ScaffSpindType() # KC: another example of function overloading?

	for j in 1:length(CM)
		nextClade = Clade(CM[j],InitTypes()) # KC: InitTypes generates random x-y pos for ancestor of the clade.
		append!(mySST.p, vcat(0,nextClade.p[2:end] .+ length(mySST.N)))
		append!(mySST.N,nextClade.N)
		append!(mySST.ptype, nextClade.ptype)
		mySST.X = vcat(mySST.X, nextClade.X .+ [curX 0 0])
		curX += len(nextClade)
	end
	if theta != 0
		myImm = SSImmigration(ceil(Int64,maximum(mySST.X[:,3])*1.2))
		mySST.p[mySST.p .!= 0] .+= length(myImm.N)
		append!(myImm.N,mySST.N)
		myImm.X = vcat(myImm.X, mySST.X)
		append!(myImm.p,mySST.p)
		append!(myImm.ptype,mySST.ptype)
		return myImm
	end
	return mySST
end


function SSImmigration(height::Int64)
	SIR = SBR*theta/alpha  # Spindle immigration rate
	DX = rand(mySeed,Geometric(SIR))+1
	if DX > height  return ScaffSpindType()  end

	curPos = [DX, height - DX]
	nextClade = ScaffSpindType()
	mySST = ScaffSpindType()

	while curPos[2] > 0
		# Need KC to add random initial types
		nextClade = Clade(InitPop,InitTypes())
		append!(mySST.p, vcat([0], nextClade.p[2:end] .+ length(mySST.N)))
		append!(mySST.N,nextClade.N)
		mySST.X = vcat(mySST.X,
			nextClade.X .+ [curPos[1] curPos[2] curPos[2]])
		append!(mySST.ptype,nextClade.ptype)
		DX = rand(mySeed,Geometric(SIR))+1
		curPos += [len(nextClade) + DX , -DX]
	end
	mySST.X[:,1] .-= len(mySST)
	return mySST
end


function RWTest(n::Int64, dx::Float64, trials::Int64, darkest::Float64 = .5, brightest::Float64 = 2.3)
	d = zeros(Float64,trials)
	moves = rand(mySeed,Uniform(-sqrt(3)*dx,sqrt(3)*dx),(n*trials,3))
	x0 = Array{Float64}(undef,3)
	x = Array{Float64}(undef,3)

	for i in 1:trials
		x0 = rand(mySeed,3)
		x = x0
		for j in (n*(i-1)+1):(n*i)
			x = abs.(x .+ moves[j,:])
			x[x .>1] = -x[x .>1] .+ 2
			if sum(x) > brightest
				x .*= (2*(brightest/sum(x))-1)
			elseif sum(x) < darkest
				x .*= ((2*darkest/sum(x)) - 1)
			end
		end
		d[i] = sum(abs.(x-x0))
	end
	return d
end


function rainbowRing(x::Float64)
	# cycles through the rainbow, in RGB, as x covers each unit interval
	x = 6*(x-floor(x))
	c = floor(Int64,x) # which side of the rainbox hexagon to be on
	x -= c	# how far to go along that side of the hexagon

	if iseven(c)
		if c==0	  		return Float64[1,x,0]
		elseif c==2  	return Float64[0,1,x]
		else  			return Float64[x,0,1]
		end
	else
		if c==1  		return Float64[1-x,1,0]
		elseif c==3  	return Float64[0,1-x,1]
		else  			return Float64[1,0,1-x]
		end
	end
end


function SSClr(mySS::ScaffSpind, meandiff::Float64 = 1.4, darkest::Float64 = .5, 
	brightest::Float64 = 2.3)
	g = Array{Int64}(undef,length(mySS.p)+1)
	g[1] = 0
	for i in 1:length(mySS.p)
		g[i+1] = g[mySS.p[i]+1] + 1
	end

	dx = min(sqrt(3/(mean(g[g .> 1])-1)) * meandiff * .42 , .5)
	# the 0.42 is from experimentation. Should give
	# avg L1 distance between RGB of descendant and of ancestor
	# approx = meandiff, assuming darkest and brightest are
	# left at default values

	Clr = rand(mySeed,Uniform(-dx,dx),(length(mySS.N),3))
	newClr = zeros(Float64,3)
	curBow = rand(mySeed)
	
	for i in 1:length(mySS.N)
		if mySS.p[i] == 0
			Clr[i,:] = rainbowRing(curBow)
			curBow += rand(mySeed,Uniform(1/6,1/2))
		else
			newClr = abs.(Clr[mySS.p[i],:] .+ Clr[i,:])
			newClr[newClr.>1] = -newClr[newClr.>1] .+ 2

			if sum(newClr) > brightest
				newClr .*= (2*(brightest/sum(newClr)) - 1)
			elseif sum(newClr) < darkest
				newClr .*= ((2*darkest/sum(newClr)) - 1)
			end
			Clr[i,:] = newClr
		end
	end
	return Clr
end


function SSLineage(SS::ScaffSpind, J::Array{Int64,1})
	lineages = Array{Int64,1}[]
	j = 0

	for i in 1:length(J)
		j = J[i]
		push!(lineages,[])
		while j > 0
			push!(lineages[i],j)
			j = SS.p[j]
		end
	end

	return lineages
end


function RainbowTypes(n::Int64)
	# random numbers are generated, then converted into rgb values via rainbow
	Clr = Array{Float64}(undef,n,3)
	Clr[:,3] = rand(n)
	for i in 1:n  Clr[i,1:3] = rainbowRing(Clr[i,3])  end
	return Clr
end


function ScaffMin(SS::ScaffSpind)
	I = filter(i -> SS.p[i]==0 , 1:length(SS.N))
	J = vcat(I[2:length(I)].-1 , length(SS.N))
	return hcat( X[I,1] , X[J,1] , X[I,2] )
end

using FileIO
using JLD


SST = ScaffSpindType(10000)
N = SST.N
X = SST.X
ptype = SST.ptype
Clr = RainbowTypes(length(N))[:,1:3]

#cd("C:\\Users\\Noah\\Dropbox\\Advising\\NF and KC - PD diffusions\\")
@save "mySSn.jld" SST Clr
#@load "NXperp.jld"

#cd("C:\\Users\\Noah\\Dropbox\\Coding\\AD\\Data\\")
#cd("C:\\Users\\Noah Forman\\Dropbox\\Coding\\AD\\")
#@save "mySSn.jld" varsthatiwanttosave
#@load "NXperp.jld"

#N = S.N; X = S.X;
#@rput N X Clr IPP slope
