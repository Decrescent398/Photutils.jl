#=
Most of this work is derived from astropy/photutils. The relevant derivations
are considered under the BSD 3-clause license. =#

module boundingboxes

export BoundingBox,
    from_float,
    center, 
    shape,
    overlap_slices,
    extent,
    union,
    intersection

using StaticArrays

"""
    A rectangular bounding box in pixel indices.

    Parameters
    ----------
    ixmin, ixmax, iymin, iymax :: int
        The bounding box pixel indices. Note that the upper values
        (``iymax`` and ``ixmax``) are exclusive as for normal slices
        in Python. The lower values (``ixmin`` and ``iymin``) must not
        be greater than the respective upper values (``ixmax`` and
        ``iymax``).
"""

struct BoundingBox

    ixmin::Int
    ixmax::Int
    iymin::Int
    iymax::Int

    function BoundingBox(ixmin::Int, ixmax::Int, iymin::Int, iymax::Int)
        ixmin > ixmax && error("ixmin must be <= ixmax")
        iymin > iymax && error("iymin must be <= iymax")
        new(ixmin, ixmax, iymin, iymax)
    end

end

end