# Colin Beckingham 2023
# Julia routine to run a Monte Carlo experiment
# to examine the output of an economic network
# where beans can be traded and money in the form
# of promissory notes can be used.

using Graphs, SimpleWeightedGraphs, GraphPlot
using Random, Compose, Colors, Combinatorics
import Cairo, Fontconfig
#
function say(wav::String,str::String)
	cmd = `flite -voice awb -t $str -o $wav`
	run(cmd)
end
#
function indexof(str,vec)
	return findall(x->x==str, vec)[1]
end
#
function validEdges(num_v,num_e)::Bool
	# when dealing with combinations of 2 items from a list
	# the number of possible combinations 
	# resolves to a simple calculation 
	# that avoids longer factorials
	possibleEdges = ((num_v-1) * num_v) / 2
	return num_e <= possibleEdges
end
#
"""
	makeEdges decides randomly which nodes are connected
	inputs: the graph, the number of nodes, the number of edges
	output: vector of edge weights
"""
function makeEdges(G,num_v,num_e)
	# decide randomly which nodes have edges
	edgecombs = 
		sort!(shuffle!(collect(combinations(1:num_v,2)))[1:num_e])
	edgeweights = []
	for i in 1:num_e
		# random weights
		wt = round(max(rand(),0.1),digits=2)
		# build an edge
		(x,y) = edgecombs[i]
		add_edge!(G, x, y, wt)
		push!(edgeweights,wt)
	end
	return edgeweights
end
#
"""
	trial() creates a new network graph
	inputs: number of nodes, number of edges
	outputs: the graph, edgeweights, stocks and node labels
"""
function trial(num_v::Int,num_e::Int)
	@assert validEdges(num_v,num_e) "init of edges wrong"
	G = SimpleWeightedGraph()
	add_vertices!(G, num_v)
	edgeweights = makeEdges(G,num_v,num_e)
	@assert length(edgeweights) == num_e "bad count of edges"
	stocks = [[Int(round(randn()*2,digits=0)),0,0,0] for i in 1:num_v]
	labels = ["$i : "*string(x) for (i,x) in enumerate(stocks)]
	return (G,edgeweights,stocks,labels)
end
#
"""
	paint() plots the network
	inputs: the graph definition, edgeweights, node labels
		and a stage (number that identifies ordinality of one
		of a series)
	output: draws the network as a png image
"""
function paint(G,edgeweights,labels,stage)
	# draw the current network
	num_v = nv(G);num_e = ne(G)
	edgelabels = [G.weights[src(e),dst(e)] for e in edges(G)]
	ran = range(0.45,stop=0.95,length=100)
	nodecolors = [RGB(rand(ran,3)...) for i in 1:num_v]
	edgecolors = [colorant"lightgreen" for i in 1:num_e]
	draw(
		PNG("test$stage.png",20cm,20cm), 
		gplot(G,
			nodelabel=labels, 
			layout=stressmajorize_layout, 
			nodefillc=nodecolors,
			edgelabel=edgelabels,
			edgestrokec=edgecolors
		)
	)
end
#
"""
	get_matches() examines the nodes that need crop seeds
	and tries to match them with other nodes that have surplus.
	Not all deficits can be matched with surplus nodes
	since there may not exist a confidence path that links them.
	
	inputs: the graph, a vector indicating those nodes in deficit,
	and those nodes that are in surplus
	output: a vector of paths linking deficit nodes with surplus
	nodes and the confidence value of the path
"""
function get_matches(G,needs,surps)
	matches = []
	for n in needs
		for s in surps
			if has_path(G,n,s)
				spath = enumerate_paths(
					dijkstra_shortest_paths(G, n),s)
				lspath = length(spath)-1
				wts = []
				for step in 1:lspath
					push!(wts, G.weights[spath[step],spath[step+1]])
				end
				push!(matches,[round(prod(wts),digits=2),[n,s]])
			end
		end
	end
	# 
	sort!(matches,rev=true)
	return matches
end
#
"""
	get_besttrade()
	
	which trades are possible and which
	involves the greatest tradeable quantity
"""
function get_besttrade(G,matches,stocks)
	besttrade = [0,[]]
	for match in matches
		confidence = match[1]
		if confidence >= 0.5
			needer,supplier = match[2]
			needamt = abs(stocks[needer][1])
			suppamt = stocks[supplier][1]
			# quantity that can be exchanged
			cantrade = min(needamt,suppamt)
			if cantrade > besttrade[1]
				besttrade = [cantrade, [needer, supplier]]
			end
		end
	end
	return besttrade
end
#
"""
	newstep() sets up the scene for functions get_matches()
	and get_besttrade()
	
	inputs: the graph, existent stock levels for each node
	output: [0,0] if no matches are possible, or
	the number of matches found and the nodes involved in the 
	best trade
"""
function newstep(G,stocks)
	crops = ["beans","potatoes","tomatoes","squash"]
	needs = [i for i in 1:nv(G) if stocks[i][indexof("beans",crops)] < 0]
	surps = [i for i in 1:nv(G) if stocks[i][indexof("beans",crops)] > 0]
	matches = get_matches(G,needs,surps)
	lm = length(matches)
	if lm == 0
		return [0,0]
	else
		besttrade = get_besttrade(G,matches,stocks)
		return [lm,besttrade]
	end
end
#
"""
	report() runs a visual graph report
	inputs: the output from a monte carlo run, and
	a number that specifies the interval at which 
	a report will be produced. If the monte carlo run
	has 1M outputs, we clearly do not need 1M reports.
"""
function report(monte,interval)
	lm = length(monte)
	for i in 0:lm
		if i % interval == 0
			divn = Int(i/interval)
			stem = "test$divn"
			img = "$stem.png"
			wav = "$stem.wav"
			mon = monte[i+1]
			ntrades = mon[1]
			if ntrades > 0
				saystuff = "in image $divn, $ntrades trades are possible"
				bestq = mon[2][1]
				if bestq > 0
					needer,supplier = mon[2][2]
					saystuff *= " and the best quantity is $bestq"
					saystuff *= " between nodes $needer and $supplier"
				else
					saystuff *= " but no trades qualify"
				end
				println(saystuff)
				say(wav,saystuff)
			else
				println("in image $img, no trades are possible")
				say(wav,"in graph $divn, no trades are possible")
			end
		end
	end
end
#
"""
	flatter() produces a regular matrix from the vector of vectors
	produced from a Monte Carlo run.
"""
function flatter(vec)
	return ( vec[2][1] > 0 ? 
				[vec[1],vec[2][1],vec[2][2][1],vec[2][2][2]] : 
				[vec[1],vec[2][1],0,0] )
end
#
"""
	main() is the entry function for simulation
	inputs are all integers
	-	num_v : number of vertices or nodes
	-	num_e : number of edges
	-	limit : max iterations for simulation
	-	interval : points at which a new plot is produced
			
	returns a matrix containing, by column
	- the number of	trades possible given the stocks and edges found 
	- the number of trades where an exchange could take place given 
	the existing confidence levels
	- the two node numbers involved in that best trade
	
	example:
	
	monte = main(12,8,1000,100)
	
	which signifies:
	create a graph with 12 vertices and 8 random 
	edges with normally distributed quantities of beans in store,
	then repeat the experiment 1000 times and record how many 
	trades were possible given inventory stocks and confidence levels.
	Store the node numbers of the nodes in the best trade and
	at intervals paint a new graphplot PNG image of the network
"""
function main(num_v::Int,num_e::Int,limit::Int,interval::Int)
	@assert interval*10 == limit "check interval for 10+1 graphplots"
	monte = []
	(G,edgeweights,stocks,labels) = trial(num_v,num_e);
	# plot and store the first network
	paint(G,edgeweights,labels,0)
	push!(monte,newstep(G,stocks))
	for i in 1:limit
		(G,edgeweights,stocks,labels) = trial(num_v,num_e);
		if i % interval == 0
			imnum = Int((i/limit)*10)
			paint(G,edgeweights,labels,imnum)
		end
		push!(monte,newstep(G,stocks))
	end
	report(monte,interval)
	return reduce(hcat,flatter.(monte))'
end
#
