
@printf("Entering startup.jl in toySims\n")

Pkg.add("PyPlot")
Pkg.add("CPUTime")
Pkg.add("JSON")
Pkg.add("SortingAlgorithms")
Pkg.add("HTTPClient")

using PyPlot
using CPUTime

# Used in the sortToys subdirectory
import Base.Sort: QuickSort, MergeSort, InsertionSort
using SortingAlgorithms #Provides HeapSort, RadixSort, TimSort
