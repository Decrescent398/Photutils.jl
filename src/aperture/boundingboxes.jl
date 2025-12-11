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

function from_float(xmin::Real, xmax::Real, ymin::Real, ymax::Real)
    """
    Return the smallest bounding box that fully contains a given
        rectangle defined by float coordinate values.

        Following the pixel index convention, an integer index
        corresponds to the center of a pixel and the pixel edges span
        from (index - 0.5) to (index + 0.5). For example, the pixel edge
        spans of the following pixels are:

        * pixel 0: from -0.5 to 0.5
        * pixel 1: from 0.5 to 1.5
        * pixel 2: from 1.5 to 2.5

        In addition, because `BoundingBox` upper limits are exclusive
        (by definition), 1 is added to the upper pixel edges. See
        examples below.

        Parameters
        ----------
        xmin, xmax, ymin, ymax : float
            The floating-point coordinates defining a rectangle. The
            lower values (``xmin`` and ``ymin``) must not be greater
            than the respective upper values (``xmax`` and ``ymax``).

        Returns
        -------
        bbox : `BoundingBox` object
            The minimal ``BoundingBox`` object fully containing the
            input rectangle coordinates.

    """

    ixmin = floor(xmin + 0.5)
    ixmax = ceil(xmax + 0.5)
    iymin = floor(ymin + 0.5)
    iymax = ceil(ymax + 0.5)

    return BoundingBox(ixmin, ixmax, iymin, iymax)

end

end