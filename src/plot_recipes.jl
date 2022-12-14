using RecipesBase
using ColorTypes
using LaTeXStrings
import PlotUtils: cgrad

#using Compat
#using Compat: range

const mygreen = RGBA{Float64}(151/255,180/255,118/255,1)
const mygreen2 = RGBA{Float64}(113/255,161/255,103/255,1)
const myblue = RGBA{Float64}(74/255,144/255,226/255,1)


@recipe function plot(b::Body{N,CS}) where {N,CS}

    x = RigidBodyTools._extend_array(b.x,CS)
    y = RigidBodyTools._extend_array(b.y,CS)

    #x = [b.x; b.x[1]]
    #y = [b.y; b.y[1]]
    linecolor --> mygreen
    if CS == ClosedBody
      fillrange --> 0
      fillcolor --> mygreen
    else
      fill := :false
    end
    aspect_ratio := 1
    legend := :none
    grid := false
    x := x
    y := y
    ()
end

@recipe function f(m::RigidBodyMotion;tmax=10)

    t = 0.0:0.01:tmax
    cx = map(ti -> translational_position(m(ti))[1],t)
    cy = map(ti -> translational_position(m(ti))[2],t)
    α = map(ti -> angular_position(m(ti)),t)
    ux = map(ti -> translational_velocity(m(ti))[1],t)
    uy = map(ti -> translational_velocity(m(ti))[2],t)
    adot = map(ti -> angular_velocity(m(ti)),t)
    xlims --> (0,tmax)
  layout := (2,3)
  grid --> :none
  linewidth --> 1
  legend --> :none
  framestyle --> :frame
  xguide --> L"t"
  xlim = min(minimum(cx),minimum(cy)),max(maximum(cx),maximum(cy))
  alim = extrema(α)
  ulim = min(minimum(ux),minimum(uy)),max(maximum(ux),maximum(uy))
  adlim = extrema(adot)
  @series begin
      subplot := 1
      ylims --> xlim
      yguide --> L"x"
      t, cx
    end
  @series begin
      subplot := 2
      ylims --> xlim
      yguide --> L"y"
      t, cy
    end
  @series begin
      subplot := 3
      ylims --> alim
      yguide --> L"\alpha"
      t, α
  end
  @series begin
      subplot := 4
      ylims --> ulim
      yguide --> L"u"
      t, ux
    end
  @series begin
      subplot := 5
      ylims --> ulim
      yguide --> L"v"
      t, uy
    end
   @series begin
      subplot := 6
      ylims --> adlim
      yguide --> L"\dot{\alpha}"
      t, adot
    end
end

function RecipesBase.RecipesBase.apply_recipe(plotattributes::Dict{Symbol, Any}, bl::BodyList)
    series_list = RecipesBase.RecipeData[]
    for b in bl
        append!(series_list, RecipesBase.RecipesBase.apply_recipe(copy(plotattributes), b) )
    end
    series_list
end
