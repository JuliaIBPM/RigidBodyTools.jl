## RigidBodyTools.jl

_Tools for creating, moving, and discretizing rigid bodies_

| Documentation | Build Status |
|:---:|:---:|
| [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaIBPM.github.io/RigidBodyTools.jl/stable) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaIBPM.github.io/RigidBodyTools.jl/dev) | [![Build Status](https://github.com/JuliaIBPM/RigidBodyTools.jl/workflows/CI/badge.svg)](https://github.com/JuliaIBPM/RigidBodyTools.jl/actions) [![Coverage](https://codecov.io/gh/JuliaIBPM/RigidBodyTools.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaIBPM/RigidBodyTools.jl) |

## About the package

The purpose of this package is to provide tools for rigid bodies with
point-discretized surfaces. It currently includes methods for

* a library of 2D surface shape definitions and associated point discretizations
* calculation of geometric properties
* collections of multiple rigid bodies
* prescribed kinematics of single or linked systems of rigid and deforming bodies, with a library of joints

These tools support a variety of classes of problems, including those that involve bodies interacting with fluids, in which we would immerse these point-discretized representations into a computational grid.

<!--
Documentation can be found at https://JuliaIBPM.github.io/RigidBodyTools.jl/latest.

**RigidBodyTools.jl** is registered in the general Julia registry. To install, enter the package manager by typing
```julia
] add RigidBodyTools
```

Then, in any version, type
```julia
julia> using RigidBodyTools
```
For examples, consult the documentation or see the example Jupyter notebooks in the Examples folder.
-->
