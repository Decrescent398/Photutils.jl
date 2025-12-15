# src/aperture/Apertures.jl
module Apertures

export BoundingBox,
        shape,
        from_float,
        intersection,
        union
        
export  ApertureMask,
        cutout,
        get_overlap_slices,
        get_values,
        multiply,
        to_image

include("boundingboxes.jl")
include("mask.jl")

end #module