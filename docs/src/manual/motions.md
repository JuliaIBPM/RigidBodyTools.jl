# Motions

```@meta
DocTestSetup = quote
  using RigidBodyTools
end
```


## Motion types
```@docs
BasicDirectMotion
RigidBodyMotion
RigidBodyMotion(::Kinematics)
RigidBodyMotion(::Any,::Any)
RigidAndDirectMotion
RigidAndDirectMotion(::Kinematics,::AbstractDirectlySpecifiedMotion)
RigidAndDirectMotion(::Kinematics,::Any,::Any)
RigidAndDirectMotion(::Any,::Any,::Any,::Any)
```

## Surface velocity functions
```@docs
surface_velocity!
surface_velocity
```

## Motion functions
```@docs
motion_state
update_body!
motion_velocity
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
