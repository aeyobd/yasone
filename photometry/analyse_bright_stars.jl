using Optim
import FITSIO

import Printf: @printf
import StatsBase: median
import ImageMorphology: label_components
using ForwardDiff
using DataFrames


function moffat(r; alpha=1, gamma=1, amplitude=1)
    return @.  amplitude * (1 + r^2/gamma^2)^(-alpha)
end


"""
    moffat_2d(x, y; params...)

A generalized Moffat 2d light profile. 

Parameters
----------
- `x`, `y`: input positions. Either real or arrays
- `x0`, `y0` : center
- `amplitude`: Overall scale. 
- `alpha`: Moffat power index.
- `gamma`: Moffal scale radius
- `ellipticity`: e = 1 - b/a  
- `theta`: Position angle in radians
"""
function moffat_2d(
    x, y;
    x0::Real = 0,
    y0::Real = 0,
    ellipticity::Real = 0,
    theta::Real = 0,
    alpha::Real = 1,
    gamma::Real = 1,
    amplitude::Real = 1,
)
    q = 1 - ellipticity
    q < 0 && throw(ArgumentError("ellipticity must be < 1"))

    # shift
    dx = @. x - x0
    dy = @. y - y0

    c = cos(theta)
    s = sin(theta)

    # rotate
    xp = @.  dx * c + dy * s
    yp = @. -dx * s + dy * c

    r = @. sqrt(xp^2 * q + (yp^2) / q)

    return moffat(r; alpha=alpha, gamma=gamma, amplitude=amplitude)
end


function moffat_2d_model(x, y;
        x_c=0, y_c=0, params::AbstractVector{<:Real}, 
        alpha::Real=1.1, sat_level::Real=Inf
    )
    scale, dx, dy, gamma = params

    prof =  moffat_2d(x .- x_c .- dx, y .- y_c .- dy; amplitude=scale,gamma=gamma, alpha=alpha)
    return min.(prof, sat_level)
end



"""
    read_image(filename)

Reads a FITS image, return an array containing the image values.
"""
function read_image(filename)
    img = FITSIO.FITS(filename, "r") do f
        read(f[1])
    end

    return img
end



"""
    get_group_centroid_areas(img)

"""
function get_group_centroids_areas(img)
    labels = label_components(img)
    n_labels = maximum(labels)

    areas = zeros(Int, n_labels)
    sumrow = zeros(Float64, n_labels)
    sumcol = zeros(Float64, n_labels)

    rowmin = fill(size(img, 1) + 1, n_labels)
    rowmax = fill(0, n_labels)
    colmin = fill(size(img, 2) + 1, n_labels)
    colmax = fill(0, n_labels)

    @inbounds for col in axes(labels, 2), row in axes(labels, 1)
        ℓ = labels[row, col]
        ℓ == 0 && continue

        areas[ℓ] += 1
        sumrow[ℓ] += row
        sumcol[ℓ] += col

        rowmin[ℓ] = min(rowmin[ℓ], row)
        rowmax[ℓ] = max(rowmax[ℓ], row)
        colmin[ℓ] = min(colmin[ℓ], col)
        colmax[ℓ] = max(colmax[ℓ], col)
    end

    centroids = Vector{Tuple{Int,Int}}(undef, n_labels)
    rel_bboxes = Vector{Tuple{Int,Int,Int,Int}}(undef, n_labels)

    # now calculate  centres and bboxes
    for k in 1:n_labels
        cx = round(Int, sumrow[k] / areas[k])
        cy = round(Int, sumcol[k] / areas[k])

        centroids[k] = (cx, cy)


        rel_bboxes[k] = (
            rowmin[k] - cx,
            rowmax[k] - cx,
            colmin[k] - cy,
            colmax[k] - cy,
        )
    end

    return centroids, areas, rel_bboxes
end


function get_cutout_bounds(img_size, centroid, bbox, radius; pad=5)
    x_c, y_c = centroid
    Nx, Ny = img_size

    x_min = min(x_c - radius, x_c + bbox[1] - pad)
    x_max = max(x_c + radius, x_c + bbox[2] + pad)
    y_min = min(y_c - radius, y_c + bbox[3] - pad)
    y_max = max(y_c + radius, y_c + bbox[4] + pad)

    # apply boundary conditions
    x_min = max(x_min, 1)
    y_min = max(y_min, 1)
    x_max = min(x_max, Nx)
    y_max = min(y_max, Ny)

    return x_min:x_max, y_min:y_max
end



function image_grid(Nx::Int, Ny::Int)
    # --- Coordinate grids ---
    x = collect(1:Nx)
    y = collect(1:Ny)
    x_img = repeat(x, 1, Ny)
    y_img = repeat(y', Nx, 1)
    return x_img, y_img
end


function clean_image_mask(image, mask)
    finite_mask = isfinite.(image)
    combined_mask = finite_mask
    if !isnothing(mask)
        combined_mask .&= .!mask
    end

    return image, combined_mask
end


"""
Fit a Moffat PSF to an image cutout by minimizing chi-squared.

Optimized parameters:
(scale, dx, dy, alpha, gamma)

Returns a dict containing the optimized parameters.
"""
function fit_moff_simple_psf_chi2(
    image::AbstractMatrix;
    x_c,
    y_c,
    uncertainty = nothing,
    mask = nothing,
    initial_scale::Real = 1.0,
    bounds = nothing,
    sat_level = 62_000,
    alpha = 1.1,
    gamma_0 = 1.0,
)

    image_clean, combined_mask = clean_image_mask(image, mask)
    if isnothing(uncertainty)
        uncertainty_clean = 1
    else
        uncertainty_clean = uncertainty[combined_mask]
    end

    Nx, Ny = size(image_clean)
    x_img, y_img = image_grid(Nx, Ny)

    if bounds === nothing
        bounds = (-x_c, Nx-x_c+1, -y_c, Ny-y_c+1)
    end

    # scale, dx, dy, width
    lower = [0,  bounds[1][1], bounds[2][1], 0]
    upper = [Inf, bounds[1][end], bounds[2][end], 50]
    p0 = [initial_scale, 0.0, 0.0, gamma_0]

    function chi2_func(params)
        x = x_img .- x_c
        y = y_img .- y_c
        model = moffat_2d_model(x, y, params=params, sat_level=sat_level, alpha=alpha)

        resid = image_clean .- model
        return sum((resid[combined_mask] ./ uncertainty_clean) .^ 2)
    end

    opt = optimize(
        chi2_func,
        lower, upper, p0,
        Fminbox(LBFGS()),
        Optim.Options(show_trace = false),
        autodiff = :forward
    )

    params_best = Optim.minimizer(opt)
    min_chi2 = Optim.minimum(opt)

    best_scale, best_dx, best_dy, best_gamma = params_best

    return Dict(
        :scale => best_scale,
        :x_shift => best_dx,
        :y_shift => best_dy,
        :alpha => alpha,
        :gamma => best_gamma,
        :chi2 => min_chi2,
        :success => Optim.converged(opt),
    )
end


function fit_all_cutouts(img, img_err, cutout_cens, cutout_bbox; radius=20, pad=5, sat_level=62_000)
    # outputs
    img_residual = copy(img)
    img_model = zeros(size(img))
    fits = []

    x_img, y_img = image_grid(size(img)...)

    for i in eachindex(cutout_cens)

        println("processing $i / $(length(cutout_cens))")
        cen = cutout_cens[i]
        cutout_bounds = get_cutout_bounds(size(img), cen, cutout_bbox[i], 
                                          radius, pad=pad)

        cutout = img_residual[cutout_bounds...]
        cutout_err = img_err[cutout_bounds...]

        x_c = cen[1] - cutout_bounds[1][1]
        y_c = cen[2] - cutout_bounds[2][1]
        bounds = (cutout_bounds[1][[1, end]] .- cen[1], 
                  cutout_bounds[2][[1,end]] .- cen[2])

        fit = fit_moff_simple_psf_chi2(cutout, x_c=x_c, y_c=y_c, 
            sat_level=sat_level, bounds=bounds, uncertainty=cutout_err)

        params_best = [fit[:scale], fit[:x_shift], fit[:y_shift], fit[:gamma]]

        model = moffat_2d_model(x_img, y_img, params=params_best, 
            x_c=cen[1], y_c=cen[2], sat_level=sat_level)

        # update results
        img_residual .-= model
        img_model .+= model
        push!(fits, fit)
    end

    return img_model, img_residual, fits
end

function filter_and_sort_regions(centroids, areas, bboxes; size_min=20, n_max=nothing)
    idxs = eachindex(areas)[areas .> size_min]

    if isnothing(n_max)
        n_max = lastindex(idxs)
    end

    idxs_sorted = idxs[sortperm(areas[idxs], rev=true)][1:n_max]

    centroids = centroids[idxs_sorted]
    areas = areas[idxs_sorted]
    bboxes = bboxes[idxs_sorted]
    @assert issorted(areas, rev=true)
    @assert all(areas .> size_min)

    return centroids, areas, bboxes
end


function run_all(img_dir; n_max = 10)
    img = read_image(img_dir * "/nobkg.fits")
    flags = read_image(img_dir * "/flag.fits")
    img_err = read_image(img_dir * "/flat_fielded.weight.fits")
    flag_sat_extended = 1

    mask = (flags .& flag_sat_extended) .> 0
    centroids, areas, bboxes = get_group_centroids_areas(mask)
    img_masked = copy(img)
    img_masked[flags .> 0] .= NaN

    centroids, areas, bboxes = filter_and_sort_regions(centroids, areas, bboxes, n_max=n_max)

    img_model, img_residual, fits = fit_all_cutouts(img_masked, img_err, centroids, bboxes)


    return img_masked, img_model, img_residual, DataFrame(fits)
end

