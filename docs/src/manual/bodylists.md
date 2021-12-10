# Body, motion, and transform lists

```@meta
DocTestSetup = quote
  using RigidBodyTools
end
```

```@docs
BodyList
getrange
Base.collect(::BodyList)
Base.sum(::AbstractVector,::BodyList,::Int)
Base.view(::AbstractVector,::BodyList,::Int)
MotionList
RigidTransformList
Base.vec(::RigidTransformList)
```
