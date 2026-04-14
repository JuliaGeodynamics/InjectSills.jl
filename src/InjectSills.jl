module InjectSills

using GeometryBasics
export Vec, Vec1, Vec2, Vec1f, Vec2f, Point, Point2, Point3, Point2f, Point3f

using GeoParams
export m, Pa, deg, NoUnits, GeoUnit, C, km, Pas, GEO_units, SI_units, NO_units 
export nondimensionalize, dimensionalize, UnitValue, Value

abstract type AbstractSill{N,_T} <: AbstractMaterialParam end
export AbstractSill

volume(s::AbstractSill) =   error("volume not implemented for type $(typeof(s))")    
area(s::AbstractSill)   =   error("area not implemented for type $(typeof(s))")                   
export volume, area

include("Utils.jl")
include("PennyShapedSills.jl")
include("Mogi_McTigue.jl")
include("FiniteEllipsoidalCavity.jl")

# Extensible entry points — implemented by optional extensions
function inject_sill! end
function plot_sill end
function surface_displacement end
export inject_sill!, plot_sill, surface_displacement

end # module InjectSills
