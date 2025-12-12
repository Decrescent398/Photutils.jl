#=
Most of this work is derived from astropy/photutils. The relevant derivations
are considered under the BSD 3-clause license. =#

module mask

export ApertureMask

using BitArrays
using .boundingboxes

struct ApertureMask
    data::Matrix{Float64}
    box::BoundingBox
    mask::BitMatrix

end

function ApertureMask(data::AbstractArray{<:Real, 2}, box::BoundingBox)

    """
    Parameters
    ----------
    data : array_like
        A 2D array representing the fractional overlap of an aperture
        on the pixel grid. This should be the full-sized (i.e., not
        truncated) array that is the direct output of one of the
        low-level `photutils.geometry` functions.

    bbox : `photutils.aperture.boundingboxes`
        The bounding box object defining the aperture minimal bounding
        box.
    """

    data_arr = Array(data)
    if size(data) != shape(box) 
        throw(ArgumentError("Mask data and bounding box must have the same size"))
    end

    mask = (data_arr .> 0)

    return ApertureMask(data_arr,
                        box,
                        mask
    )

end

end