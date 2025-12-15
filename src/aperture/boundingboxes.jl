#=
Most of this work is derived from astropy/photutils. The relevant derivations
are considered under the BSD 3-clause license. =#

#=
Most of this comes from 
https://github.com/astropy/photutils/blob/main/photutils/aperture/bounding_box.py =#

#= 
Major Improvements:
1. Static typing
2. Inline functions
3. UnitRange allocation
4. No temp objects 
5. Type Stability=#

#TODO: Write Plotting functions with photutils.rectangle and makie

module boundingboxes

export BoundingBox,
    from_float,
    center, 
    shape,
    get_overlap_slices,
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

@inline function center(box::BoundingBox)

    """
    The ``(y, x)`` center of the bounding box as StaticVector
    """

    y_center = 0.5 * (box.iymax - 1 + box.iymin)
    x_center = 0.5 * (box.ixmax - 1 + box.ixmin)

    return @SVector [y_center, x_center]

end

"""
The ``(ny, nx)`` shape of the bounding box as tuple.
"""

@inline shape(box::BoundingBox) = (box.iymax - box.iymin, 
                                    box.ixmax - box.ixmin)

@inline function extent(box::BoundingBox)

    """
    The extent of the mask, defined as the ``(xmin, xmax, ymin,
    ymax)`` bounding box from the bottom-left corner of the lower-
    left pixel to the upper-right corner of the upper-right pixel.

    The upper edges here are the actual pixel positions of the
    edges, i.e., they are not "exclusive" indices used for python
    indexing. The extent is useful for plotting the bounding box.
    """

    return (box.ixmin - 0.5,
            box.ixmax - 0.5,
            box.iymin - 0.5,
            box.iymax - 0.5)
    
end

function get_overlap_slices(box::BoundingBox, ny::Int, nx::Int)

    """
    Get slices for the overlapping part of the bounding box and a 2D
    array.

    Parameters
    ----------
    ny, nx: Shape of Array

    Returns
    -------
    slices_large : tuple of slices or `None`
        A tuple of slice objects for each axis of the large array,
        such that ``large_array[slices_large]`` extracts the region
        of the large array that overlaps with the small array.
        `None` is returned if there is no overlap of the bounding
        box with the given image shape.

    slices_small : tuple of slices or `empty`
        A tuple of slice objects for each axis of an array enclosed
        by the bounding box such that ``small_array[slices_small]``
        extracts the region that is inside the large array. `empty`
        is returned if there is no overlap of the bounding box with
        the given image shape.
    """

    xmin, xmax = box.ixmin, box.ixmax
    ymin, ymax = box.iymin, box.iymax

    #No overlap
    if xmin >= nx || ymin >= ny || xmax <= 0 || ymax <= 0
        return (empty, empty), (empty, empty)
    end
    
    y_slices_large = max(ymin, 0) + 1: min(ymax, ny)
    x_slices_large = max(xmin, 0) + 1: min(xmax, nx)

    y_slices_small = max(1 - ymin, 1) + 1: min(ymax - ymin, ny - ymin)
    x_slices_small = max(1 - xmin, 0) + 1: min(xmax - xmin, nx - xmin)

    return (y_slices_large, x_slices_large), (y_slices_small, x_slices_small)

end

@inline function union(box::BoundingBox, other::BoundingBox)

    """
    Return a `BoundingBox` representing the union of this
    `BoundingBox` with another `BoundingBox`.

    Parameters
    ----------
    other : `BoundingBox`
        The `BoundingBox` to join with this one.

    Returns
    -------
    result : `BoundingBox`
        A `BoundingBox` representing the union of the input
        `BoundingBox` with this one.
    """
    
    return BoundingBox(min(box.ixmin, other.ixmin),
                        max(box.ixmax, other.ixmax),
                        min(box.iymin, other.iymin),
                        max(box.iymax, other.iymax))

end

@inline function intersection(box::BoundingBox, other::BoundingBox)

    """
    Return a `BoundingBox` representing the intersection of this
    `BoundingBox` with another `BoundingBox`.

    Parameters
    ----------
    other : `BoundingBox`
        The `BoundingBox` to intersect with this one.

    Returns
    -------
    result : `BoundingBox`
        A `BoundingBox` representing the intersection of the input
        `BoundingBox` with this one or an empty BoundingBox
    """

    ixmin = max(box.ixmin, other.ixmin)
    ixmax = min(box.ixmax, other.ixmax)
    iymin = max(box.iymin, other.iymin)
    iymax = min(box.iymax, other.iymax)

    if ixmax < ixmin || iymax < iymin
        return BoundingBox(1, 0, 1, 0) #empty range
    end

    return BoundingBox(ixmin, ixmax, iymin, iymax)

end

end #module