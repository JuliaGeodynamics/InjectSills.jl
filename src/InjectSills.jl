module InjectSills

using GeometryBasics
export Vec, Vec1, Vec2, Vec1f, Vec2f, Point, Point2, Point3, Point2f, Point3f

using GeoParams
export m, Pa, deg, NoUnits, GeoUnit, C, km, Pas, GEO_units, SI_units, NO_units 
export nondimensionalize, dimensionalize, UnitValue, Value

abstract type AbstractSill{N,_T} <: AbstractMaterialParam end  
export AbstractSill

include("Utils.jl")
include("PennyShapedSills.jl")

end # module InjectSills
