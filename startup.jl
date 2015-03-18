
@printf("Entering sortToys startup.jl\n")

using PyPlot
#import PyPlot; const plt = PyPlot

Pkg.add("SortingAlgorithms")
using CPUTime
    
import Base.Sort: QuickSort, MergeSort, InsertionSort
using SortingAlgorithms #Provides HeapSort, RadixSort, TimSort
