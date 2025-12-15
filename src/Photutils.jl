module Photutils

using Reexport

include("aperture/Apertures.jl")

@reexport using .Apertures

end #module