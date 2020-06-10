# RigidBodyTools.jl

*Tools for creating, moving, and discretizing rigid bodies*

The purpose of this package is to provide tools for rigid bodies with
point-discretized surfaces. It includes methods for

* a library of surface shape definitions and associated point discretizations
* calculation of geometric properties
* rigid-body motion and transformation of surface points
* collections of multiple rigid bodies

## Installation

This package works on Julia `1.0` and above and is registered in the general Julia registry. To install from the REPL, type
e.g.,
```julia
] add RigidBodyTools
```

Then, in any version, type
```julia
julia> using RigidBodyTools
```

The plots in this documentation are generated using [Plots.jl](http://docs.juliaplots.org/latest/).
You might want to install that, too, to follow the examples.
