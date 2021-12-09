# Index

```@meta
DocTestSetup = quote
  using RigidBodyTools
end
```

## Shapes
```@docs
BasicBody
Circle
Ellipse
NACA4
Plate
Rectangle
SplinedBody
Square
```

## Body list functions
```@docs
BodyList
getrange
Base.collect(::BodyList)
Base.sum(::AbstractVector,::BodyList,::Int)
Base.view(::AbstractVector,::BodyList,::Int)
```


## Rigid transformations
```@docs
RigidTransform
RigidTransformList
Base.vec(::RigidTransform)
Base.vec(::RigidTransformList)
```


## Surface functions
```@docs
centraldiff
Base.diff(::Body)
Base.diff(::BodyList)
dlength
dlengthmid
Base.length(::Body)
midpoints
normal
normalmid
```

## Motion functions
```@docs
BasicDirectMotion
RigidBodyMotion
MotionList
motion_state
motion_velocity
surface_velocity!
surface_velocity
```

## Rigid body kinematics types
```@docs
Kinematics
Oscillation
OscillationX
OscillationY
OscillationXY
PitchHeave
Pitchup
RotationalOscillation
SwitchedKinematics
```
