
@printf("Entering startup.jl in %s\n",pwd())

Pkg.add("PyPlot")
Pkg.add("CPUTime")
Pkg.add("SortingAlgorithms")

using PyPlot
using CPUTime

# Used in the sortToys subdirectory
import Base.Sort: QuickSort, MergeSort, InsertionSort
using SortingAlgorithms #Provides HeapSort, RadixSort, TimSort
