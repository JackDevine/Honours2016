# Honours2016
Honours project on Brownian dynamics, the source code is in the `src` directory and includes unit tests for the important functions in the project.

In the `notebooks` directory, you will find some notebooks that walk through examples on how to use the code and some important results from the project.

If you want to run the code on your own machine, then you will need to install [Julia](http://julialang.org/), once you have Julia, there are some packages that you will need to run the code.

```julia
  Pkg.add("PyPlot")  # Used for plotting.
  Pkg.add("ForwardDiff")  # Used for automatically finding the derivatives of functions.
  Pkg.add("SymPy")  # Used to obtain analytical solutions to the Smoluchowski equation.
```

Once you have these prerequisites, simply clone this repo and start playing. You may need to change some paths in some of the notebooks, but that won't be too hard.
