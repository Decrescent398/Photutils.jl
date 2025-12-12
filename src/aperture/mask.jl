#=
Most of this work is derived from astropy/photutils. The relevant derivations
are considered under the BSD 3-clause license. =#

module mask

export ApertureMask

using BitArrays
using .boundingboxes: BoundingBox, get_overlap_slices, shape

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

    mask = BitMatrix(data_arr .== 0)

    return ApertureMask(data_arr,
                        box,
                        mask
    )

end


function get_overlap_slices(mask::ApertureMask, ny::Int, nx::Int)

    """
    Get slices for the overlapping part of the bounding box and a 2D
    array.

    Parameters
    ----------
    ny, nx: size of Array

    Returns
    -------
    slices_large : tuple of slices or `None`
        A tuple of slice objects for each axis of the large array,
        such that ``large_array[slices_large]`` extracts the region
        of the large array that overlaps with the small array.
        `None` is returned if there is no overlap of the bounding
        box with the given image size.

    slices_small : tuple of slices or `None`
        A tuple of slice objects for each axis of an array enclosed
        by the bounding box such that ``small_array[slices_small]``
        extracts the region that is inside the large array. `None`
        is returned if there is no overlap of the bounding box with
        the given image size.
    """

    return get_overlap_slices(mask.box, ny, nx)

end

function to_image(mask::ApertureMask, ny::Int, nx::Int; dtype=Float64)

    """
    Return an image of the mask in a 2D array of the given size,
    taking any edge effects into account.

    Parameters
    ----------
    ny, nx: size of Array

    dtype : data-type, optional
        The desired data type for the array. This should be a
        floating data type if the `ApertureMask` was created with
        the "exact" or "subpixel" mode, otherwise the fractional
        mask weights will be altered. An integer data type may be
        used if the `ApertureMask` was created with the "center"
        mode.

    Returns
    -------
    result : A 2D array of the mask.
    """

    slices_large, slices_small = get_overlap_slices(mask.box, ny, nx)
    if slices_small === nothing
        return nothing
    end

    image = zeros(dtype, ny, nx)
    @views image[slices_large] .= mask.data[slices_small]

    return image

end

function cutout(mask::ApertureMask, data::AbstractArray{<:Real,2}; fill_value=0.0, copy=false)

    """
    Create a cutout from the input data over the mask bounding box

    Parameters
    ----------
    data : array_like
        A 2D array on which to apply the aperture mask.

    fill_value : float, optional
        The value used to fill pixels where the aperture mask does
        not overlap with the input ``data``. The default is 0.

    copy : bool, optional
        If `True` then the returned cutout array will always be hold
        a copy of the input ``data``. If `False` and the mask is
        fully within the input ``data``, then the returned cutout
        array will be a view into the input ``data``. In cases where
        the mask partially overlaps or has no overlap with the input
        ``data``, the returned cutout array will always hold a copy
        of the input ``data`` (i.e., this keyword has no effect).

    Returns
    -------
    result : `array` or `None`
        A 2D array cut out from the input ``data`` representing
        the same cutout region as the aperture mask. If there is a
        partial overlap of the aperture mask with the input data,
        pixels outside the data will be assigned to ``fill_value``.
        `None` is returned if there is no overlap of the aperture
        with the input ``data``.
    """
    
    ny, nx = size(data)
    slices_large, slices_small = get_overlap_slices(mask, ny, nx)

    if slices_small === nothing
        return nothing
    end

    cutout_size = (length(first(slices_small)), length(last(slices_small)))

    if cutout_size == size(mask.data)

        cutout = @view data[slices_large]

        if copy
            cutout = copy(cutout)
        end

        return cutout

    end

    dtype = promote_type(typeof(fill_value), eltype(data))
    cutout = zeros(dtype, size(mask.data))
    cutout[:] = fill_value
    cutout[slices_small] = data[slices_large]

    return cutout

end

function multiply(mask::ApertureMask, data::AbstractArray{<:Real,2}; fill_value=0.0)

    """
    Multiply the aperture mask with the input data, taking any edge
    effects into account.

    The result is a mask-weighted cutout from the data.

    Parameters
    ----------
    data : array
        The 2D array to multiply with the aperture mask.

    fill_value : float, optional
        The value is used to fill pixels where the aperture mask
        does not overlap with the input ``data``. The default is 0.

    Returns
    -------
    result : `array` or `None`
        A 2D mask-weighted cutout from the input ``data``. If
        there is a partial overlap of the aperture mask with the
        input data, pixels outside the data will be assigned to
        ``fill_value`` before being multiplied with the mask. `None`
        is returned if there is no overlap of the aperture with the
        input ``data``.
    """

    c = cutout(mask, data; fill_value=fill_value, copy=False)
    
    if c === nothing
        return nothing
    end

    weighted_cutout = c .* mask.data
    weighted_cutout[mask.mask] = fill_value

    return weighted_cutout

end

function _get_overlap_cutouts(mask::ApertureMask, ny::Int, nx::Int; _mask=nothing)

    """
    Get the aperture mask weights, pixel mask, and slice for the
    overlap with the input size.

    If input, the ``_mask`` is included in the output pixel mask
    cutout.

    Parameters
    ----------
    ny, nx : The size of data.

    _mask : array_like (bool), optional
        A boolean mask with the same size as ``size`` where a
        `True` value indicates a masked pixel.

    Returns
    -------
    slices_large : tuple of slices or `None`
        A tuple of slice objects for each axis of the large array
        of given ``size``, such that ``large_array[slices_large]``
        extracts the region of the large array that overlaps with
        the small array. `None` is returned if there is no overlap
        of the bounding box with the given image size.

    aper_weights: 2D float `array`
        The cutout aperture mask weights for the overlap.

    pixel_mask: 2D bool `array`
        The cutout pixel mask for the overlap.

    Notes
    -----
    This method is separate from ``get_values`` to facilitate
    applying the same slices, aper_weights, and pixel_mask to
    multiple associated arrays (e.g., data and error arrays). It is
    used in this way by the `PixelAperture.do_photometry` method.
    """
    
    if (_mask !== nothing) && (size(_mask) != (ny, nx))
        error("Mask and data size must be the same")
    end

    slc_large, slc_small = get_overlap_slices(mask, ny, nx)
    
    if slc_large === nothing
        return nothing, nothing, nothing
    end

    aper_weights = mask.data[slc_small]
    pixel_mask = aper_weights .> 0.0

    if _mask !== nothing 
        pixel_mask .&= .!_mask[slc_large]
    end

    return slc_large, aper_weights, pixel_mask

end

function get_values(mask::ApertureMask, data::AbstractArray{<:Real,2}; _mask=nothing)

    """
    Get the mask-weighted pixel values from the data as a 1D array.

    If the ``ApertureMask`` was created with ``method='center'``,
    (where the mask weights are only 1 or 0), then the returned
    values will simply be pixel values extracted from the data.

    Parameters
    ----------
    data : array
        The 2D array from which to get mask-weighted values.

    mask : array_like (bool), optional
        A boolean mask with the same shape as ``data`` where a
        `True` value indicates the corresponding element of ``data``
        is not returned in the result.

    Returns
    -------
    result : `array`
        A 1D array of mask-weighted pixel values from the input
        ``data``. If there is no overlap of the aperture with the
        input ``data``, the result will be an empty array with shape
        (0,).
    """
    
    ny, nx = size(data)
    slc_large, aper_weights, pixel_mask = _get_overlap_cutouts(mask, ny, nx; _mask=_mask)

    if slc_large === nothing
        return Float64[]
    end
    
    return (aper_weights .* data[slc_large])[pixel_mask]

end

end