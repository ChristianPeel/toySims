function sortToy(Ns::Int,N::Array{Int,1},L::Array{Int,1},dataType::Array{DataType,1},distType)
#sortToy
#  Do Monte-Carlo sims of sorting techniques. N data points of type
#  dataType and width L comprise each set of data to be sorted. Ns
#  different sets are timed, with the average time shown.
#
# Examples
#  # As function of Integer width
#  sortToy(40,2.^[10;],[40;],[Int8,Int16,Int32,Int64,Int128],"random")
#  # As function of N, with data sorted in reverse, Tim is best
#  sortToy(40,2.^[2:10;],[10;],[AbstractString],"reverse")
#  # As function of N, with data sorted in reverse, Radix is best 
#  sortToy(100,2.^[2:14;],[50;],[Int8],"random")

#
# By Christian Peel  (chris.peel@ieee.org)
# Last Modified: Mon 9 Mar 15, 10:09pm by peel

# using PyPlot
# Pkg.add("SortingAlgorithms")
# using CPUTime
    
# import Base.Sort: QuickSort, MergeSort, InsertionSort
# using SortingAlgorithms #Provides HeapSort, RadixSort, TimSort

sorts = ["Heap", "Radix", "Tim", "Quick", "Merge", "Default","Insertion"]
if any(N.>1024)
    sorts = sorts[sorts.!="Insertion"];
end

println("      Ns      M      L   dataType")

# Not using the following 'getAlgs()' method because it doesn't allow
# me to easily plot genie-aided receivers.

out = cell(0)
if length(N)>1
  for s = 1:length(N)
    push!(out,sortSim(Ns,N[s],L[1],dataType[1], distType, sorts));
  end
  plotfun = loglog;
  xval = N;
  xlab = "Set Size";
  tstr = @sprintf("Ns=%d,L=%d,dataType=%s",Ns,L[1],string(dataType[1]));
elseif length(L)>1
  for s = 1:length(L)
    push!(out,sortSim(Ns,N[1],L[s],dataType[1], distType, sorts));
  end
  plotfun = semilogy;
  xval = L;
  xlab = "Max";
  tstr = @sprintf("Ns=%d,N=%d,dataType=%s",Ns,N[1],string(dataType[1]));
elseif length(dataType)>1
  for s = 1:length(dataType)
    push!(out,sortSim(Ns,N[1],L[1],dataType[s], distType, sorts));
  end
  plotfun = semilogy;
  xval = 1:length(dataType)
  xlab = "Max";
  tstr = @sprintf("Ns=%d,N=%d,L=%d",Ns,N[1],L[1]);
end

pColor = {"r>-","bo--","kx-.","gd-","c^--","m*-.",
          "rs--","gp-.","bv-","kh--","c+-.","m.-",};
pIdx   = 1;

Nout = size(out,1);

clf()
times = zeros(Nout);
for a=1:length(out[1][2])
    for k=1:Nout
        times[k] = out[k][1][a];
#        mi(k) = out(k).mi(a);
    end
#    figure(1)
#println(" finished times")
    plotfun(xval,times,pColor[pIdx],label=out[1][2][a]);
    #semilogy(xval,times,pColor[pIdx],label=out[1][2][a]);
    hold(true);
##     figure(2)
##     plot(xval,mi,pColor[pIdx]);  hold(true);
    if pIdx == length(pColor)
         pIdx = 1;
    else
         pIdx = pIdx + 1;
    end
end

hold(false);
##figure(1)
xlabel(xlab);
ylabel("execution time");
legend(loc=2);
grid(which="both",axis="y")
grid(which="major",axis="x")
title(tstr)

return
end

#######################################################################
function sortSim(Ns,N,L,dataType,distType,sorts)
sortDict = Dict("Insertion" => d->sort(d,alg=InsertionSort),
                   # Worst with large N, ok for small
                "Quick"     => d->sort(d,alg=QuickSort),
                   # 
                "Merge"     => d->sort(d,alg=MergeSort),
                "Heap"      => d->sort(d,alg=HeapSort),
                "Radix"     => d->sort(d,alg=RadixSort),
                   # Best for Int8,largeN, worst with small N
                "Tim"       => d->sort(d,alg=TimSort),
                   # Best w sorted||reverse sorted AbstractString, large N
                "Default"   => d->sort(d)); 

if dataType<:Integer
    randfn! = d -> rand!(d, 1:L);
elseif dataType<:AbstractString
    sorts = sorts[sorts.!="Radix"];
    randfn! = d -> (for i = 1:length(d); d[i] = randstring(L); end; d);
elseif dataType<:FloatingPoint
    randfn! = d -> L*rand!(d::AbstractArray{dataType,1});
else
    error("Unknown Type $string(dataType)")
end
distDict = Dict("random"=> randfn!,
                "sorted"=> d -> begin
                        randfn!(d);
                        sort!(d);
                    end,
                "reverse"=> d -> begin
                        randfn!(d);
                        sort!(d);
                        reverse!(d);
                    end,
                "ones"=> d -> ones(d),
                "rand4"=> d -> begin
                        randfn!(d);
                        data4=data[rand(1:4,N)]
                    end,
                "sortedW10rand"=> d -> begin
                        randfn!(d);
                        sort!(d);
                        d[end-9:end] = randfn!(Array(dataType,10));
                        return d
                    end,
                );
randf! = distDict[distType];

Nsort = length(sorts);
times = zeros(Nsort,Ns);
algNames = cell(Nsort)

data = Array(dataType,N);

for ix = 1:Ns
    CPUtic();
    data = randf!(data);
    times[end,ix] = CPUtoq();
    algNames[end] = "randfn";

    #
    # Sort algorithms
    #
    for ax = 1:length(sorts)
        CPUtic();
        sorted = sortDict[sorts[ax]](data);
        times[ax,ix] = CPUtoq();
        algNames[ax] = sorts[ax];
    end

end
@printf("%8d %6d %6d %12s", Ns,  N,  L, string(dataType))

#outtimes = mean(times,2);
outtimes = zeros(Nsort,1);
for ax = 1:Nsort
    stimes = sort(vec(times[ax,:]));
    outtimes[ax] = mean(stimes[1:floor(Ns*2/3)]);
    if ax<8
        @printf("    %7.4f",outtimes[ax]*1000)
    end 
end
@printf("\n")

return outtimes, algNames
end
