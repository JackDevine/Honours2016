# Honours2016
Honours project on Brownian dynamics, the source code is in the `src` directory and includes unit tests for the important functions in the project, the most important file is called finiteDifferences.jl, here you will find my implementation pf finite differences.

In the `notebooks` directory, you will find some notebooks that walk through examples on how to use the code and some important results from the project.

If you want to run the code on your own machine, then you will need to install [Julia](http://julialang.org/), once you have Julia, there are some packages that you will need to run the code, simply go to the Julia command line and type.

```julia
julia> Pkg.add("PyPlot")  # Used for plotting.
julia> Pkg.add("ForwardDiff")  # Used for automatically finding the derivatives of functions.
julia> Pkg.add("SymPy")  # Used to obtain analytical solutions to the Smoluchowski equation.
```

Once you have these prerequisites, simply clone this repo and start playing. You may need to change some paths in some of the notebooks, but that won't be too hard.
