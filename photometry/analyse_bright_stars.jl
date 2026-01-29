using Optim
import FITSIO

import Printf: @printf
import StatsBase: median
import ImageMorphology: label_components
using ADTypes: AutoForwardDiff
using DataFrames
using PyFITS


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
        params::AbstractVector{<:Real}, 
        alpha=1.1
        )
    scale, dx, dy, gamma = params

    prof =  moffat_2d(x .- dx, y .- dy; amplitude=scale, gamma=gamma, alpha=alpha)
    return prof
end

function moffat_2d_gradient(x, y; params::AbstractVector{<:Real}, alpha=1.1)
    scale, dx, dy, gamma = params
    prof = moffat_2d_model(x, y, params=params)

    r2 = @. (x-dx)^2 + (y-dy)^2
    d_scale = @. prof / scale
    d_dx = @. scale * (1 + r2/gamma^2)^(-alpha-1) * alpha * 2*(x-dx)/gamma
    d_dy = @. scale * (1 + r2/gamma^2)^(-alpha-1) * alpha * 2*(y-dy)/gamma
    d_gamma = @. scale * (1 + r2/gamma^2)^(-alpha-1) * alpha * r2/gamma^3 * 2

    return [d_scale, d_dx, d_dy, d_gamma]
end



function moffat_2d_model_alpha(x, y;
        params::AbstractVector{<:Real}, 
        )
    scale, dx, dy, gamma, alpha = params

    prof =  moffat_2d(x .- dx, y .- dy; amplitude=scale, gamma=gamma, alpha=alpha)
    return prof
end

function moffat_2d_alpha_gradient(x, y; params::AbstractVector{<:Real},)
    scale, dx, dy, gamma, alpha = params
    prof = moffat_2d_model(x, y, params=params)

    r2 = @. (x-dx)^2 + (y-dy)^2
    d_scale = @. prof / scale
    d_dx = @. scale * (1 + r2/gamma^2)^(-alpha-1) * alpha * 2*(x-dx)/gamma
    d_dy = @. scale * (1 + r2/gamma^2)^(-alpha-1) * alpha * 2*(y-dy)/gamma
    d_gamma = @. scale * (1 + r2/gamma^2)^(-alpha-1) * alpha * r2/gamma^3 * 2
    d_alpha = @. -prof * log(r2 / gamma^2 + 1)

    return [d_scale, d_dx, d_dy, d_gamma, d_alpha]
end

function moffat_2d_hessian(x, y; params::AbstractVector{<:Real}, alpha::Real=1.1)
    scale, dx, dy, gamma = params
    prof = moffat_2d_model(x, y, params=params, alpha=alpha)

    r2 = @. (x-dx)^2 + (y-dy)^2
    d_scale_scale = @. -prof / scale^2
    d_dx_scale = @.  (1 + r2/gamma^2)^(-alpha-1) * alpha * 2*(x-dx)/gamma
    d_dy_scale = @. (1 + r2/gamma^2)^(-alpha-1) * alpha * 2*(y-dy)/gamma
    d_gamma_scale = @. (1 + r2/gamma^2)^(-alpha-1) * alpha * r2/gamma^3 * 2

    return
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
    write_image(filename)

Reads a FITS image, return an array containing the image values.
"""
function write_image(filename, img)
    if isfile(filename)
        rm(filename)
    end
    
    img = FITSIO.FITS(filename, "w") do f
        write(f, img)
    end
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
    bboxes = Vector{Tuple{Int,Int,Int,Int}}(undef, n_labels)

    # now calculate  centres and bboxes
    for k in 1:n_labels
        cx = round(Int, sumrow[k] / areas[k])
        cy = round(Int, sumcol[k] / areas[k])

        centroids[k] = (cx, cy)


        bboxes[k] = (
            rowmin[k],
            rowmax[k],
            colmin[k],
            colmax[k],
        )
    end

    return centroids, areas, bboxes
end


function get_cutout_bounds(img_size, centroid, bbox, radius; pad=20)
    x_c, y_c = centroid
    Nx, Ny = img_size

    x_min = min(x_c - radius, bbox[1] - pad)
    x_max = max(x_c + radius, bbox[2] + pad)
    y_min = min(y_c - radius, bbox[3] - pad)
    y_max = max(y_c + radius, bbox[4] + pad)

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
    fit_model_psf_chi2(model, image; params)

Fit a function model to the psf of a given cutout image. 
The model should take x and y positions and a parameter vector for the first three
arguments. Additional keywords are passed to the model.
Return a dict containing the best-fit parameters and basic convergence information

# Parameters
## Required
- `lower` The lower parameter bounds
- `upper` The upper parameter bounds
- `p0` An initial guess for the parameters,
## Optional
- `x0`: The index corresponding to the first pixel row
- `y0`: The index corresponding to the first pixel column
- `mask`: A mask image
- `uncertainty`: An uncertainty image
- `chi2_f`: A functional form to transform z-scores before summing. For chi2, this is a quadratic
- `chi2_scale`: Scale the z-scores before passing to chi2_f?

"""
function fit_model_grad_psf_chi2(
    model, model_gradient,
    image::AbstractMatrix;
    uncertainty = nothing,
    mask = nothing,
    lower,
    upper,
    p0,
    x0 = 1,
    y0 = 1,
    chi2_f = x->x^2, 
    chi2_f_der = x->2x, 
    chi2_scale::Real = 1,
    iterations::Int = 100,
    x_reltol = 1e-10,
    f_reltol = 1e-10,
    kwargs...
)

    image_clean, combined_mask = clean_image_mask(image, mask)
    if isnothing(uncertainty)
        uncertainty_clean = 1
    else
        uncertainty_clean = uncertainty[combined_mask]
    end

    Nx, Ny = size(image_clean)
    x_img, y_img = image_grid(Nx, Ny)
    x_img = x_img[combined_mask] .+ x0 .- 1
    y_img = y_img[combined_mask] .+ y0 .- 1
    image_clean = image_clean[combined_mask]

    function chi2_func(params)
        img_model = model(x_img, y_img, params=params; kwargs...)

        resid = image_clean .- img_model
        return sum(@. chi2_f(resid / uncertainty_clean / chi2_scale))
    end

    function gradient_chi2!(G, params)
        img_grad = model_gradient(x_img, y_img, params=params; kwargs...)

        img_model = model(x_img, y_img, params=params; kwargs...)
        resid = image_clean .- img_model
        
        for i in eachindex(params)
            G[i] = sum(@. -chi2_f_der(resid / uncertainty_clean / chi2_scale) * img_grad[i] / uncertainty_clean / chi2_scale) 
        end
        return G
    end

    options = Optim.Options(show_trace = false,
        x_reltol = x_reltol,
        f_reltol = f_reltol,
        iterations = iterations
    )
    opt = optimize(
        chi2_func,
        gradient_chi2!,
        lower, upper, p0, 
        Fminbox(LBFGS()),
        options,
        autodiff = AutoForwardDiff()
    )

    params_best = Optim.minimizer(opt)
    min_chi2 = Optim.minimum(opt)

    return Dict(
        :params => params_best,
        :chi2 => min_chi2,
        :success => Optim.converged(opt),
    )
end

"""
    fit_model_psf_chi2(model, image; params)

Fit a function model to the psf of a given cutout image. 
The model should take x and y positions and a parameter vector for the first three
arguments. Additional keywords are passed to the model.
Return a dict containing the best-fit parameters and basic convergence information

# Parameters
## Required
- `lower` The lower parameter bounds
- `upper` The upper parameter bounds
- `p0` An initial guess for the parameters,
## Optional
- `x0`: The index corresponding to the first pixel row
- `y0`: The index corresponding to the first pixel column
- `mask`: A mask image
- `uncertainty`: An uncertainty image
- `chi2_f`: A functional form to transform z-scores before summing. For chi2, this is a quadratic
- `chi2_scale`: Scale the z-scores before passing to chi2_f?

"""
function fit_model_psf_chi2(
    model,
    image::AbstractMatrix;
    uncertainty = nothing,
    mask = nothing,
    lower,
    upper,
    p0,
    x0 = 1,
    y0 = 1,
    chi2_f = x->x^2, 
    chi2_f_der = x->2x, 
    chi2_scale::Real = 1,
    iterations::Int = 100,
    x_reltol = 1e-10,
    f_reltol = 1e-10,
    kwargs...
)

    image_clean, combined_mask = clean_image_mask(image, mask)
    if isnothing(uncertainty)
        uncertainty_clean = 1
    else
        uncertainty_clean = uncertainty[combined_mask]
    end

    Nx, Ny = size(image_clean)
    x_img, y_img = image_grid(Nx, Ny)
    x_img = x_img[combined_mask] .+ x0 .- 1
    y_img = y_img[combined_mask] .+ y0 .- 1
    image_clean = image_clean[combined_mask]

    function chi2_func(params)
        img_model = model(x_img, y_img, params=params; kwargs...)

        resid = image_clean .- img_model
        return sum(@. chi2_f(resid / uncertainty_clean / chi2_scale))
    end

    options = Optim.Options(show_trace = false,
        x_reltol = x_reltol,
        f_reltol = f_reltol,
        iterations = iterations
    )
    opt = optimize(
        chi2_func,
        lower, upper, p0, 
        Fminbox(LBFGS()),
        options,
        autodiff = AutoForwardDiff()
    )

    params_best = Optim.minimizer(opt)
    min_chi2 = Optim.minimum(opt)

    return Dict(
        :params => params_best,
        :chi2 => min_chi2,
        :success => Optim.converged(opt),
    )
end

function fit_all_cutouts(img, img_err, cutout_cens, cutout_bbox; 
        radius=30, pad=30, xi=10, alpha=1.1)
    # outputs
    img_residual = copy(img)
    img_model = zeros(size(img))
    fits = []

    x_img, y_img = image_grid(size(img)...)
    bkg_rms = median(img_err)
    @info "median err: $bkg_rms"

    for i in eachindex(cutout_cens)

        println("processing $i / $(length(cutout_cens))")
        cen = cutout_cens[i]
        bbox = cutout_bbox[i]
        cutout_bounds = get_cutout_bounds(size(img), cen, bbox, radius, pad=pad)

        cutout = img_residual[cutout_bounds...]
        cutout_err = img_err[cutout_bounds...]

        param_names = [:scale, :dx, :dy, :alpha, :gamma]
        x_low = cutout_bounds[1][1]
        y_low = cutout_bounds[2][1]
        lower = [0., bbox[1], bbox[3], 0.01]
        upper = [Inf, bbox[2], bbox[4], 30.0]
        p0 = [1e5, cen[1], cen[2], 0.065]

        #fit = fit_model_grad_psf_chi2(moffat_2d_model, moffat_2d_gradient, cutout,
        #                         uncertainty=cutout_err,
        #                         x0=x_low, y0=y_low,
        #                         lower=lower, upper=upper, p0=p0,
        #                         chi2_f = semilog2,
        #                         chi2_f_der = semilog2_der
        #                        #chi2_scale = 5,
        #                        )
        fit = fit_model_psf_chi2(moffat_2d_model, cutout,
                                 uncertainty=cutout_err,
                                 x0=x_low, y0=y_low,
                                 lower=lower, upper=upper, p0=p0,
                                 chi2_f = semilog2,
                                 chi2_scale = xi,
                                 alpha = alpha
                                )

        params_best = fit[:params]

        model = moffat_2d_model(x_img, y_img, params=params_best)


        if !fit[:success]
            @info "fit failed to converge"
        end

        img_residual .-= model
        img_model .+= model
        push!(fits, fit)

    end

    return img_model, img_residual, fits
end

function filter_and_sort_regions(centroids, areas, bboxes; size_min=50, n_max=nothing)
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


function run_all(imgname, weightname, flagsname; n_max = nothing)
    img = read_image(imgname)
    flags = read_image(flagsname)
    img_err = 1 ./ sqrt.(read_image(weightname))
    flag_sat_extended = 1

    mask = (flags .& flag_sat_extended) .> 0
    centroids, areas, bboxes = get_group_centroids_areas(mask)
    img_masked = copy(img)
    img_masked[flags .> 0] .= NaN

    centroids, areas, bboxes = filter_and_sort_regions(centroids, areas, bboxes, n_max=n_max)

    img_model, img_residual, fits = fit_all_cutouts(img_masked, img_err, centroids, bboxes)


    return img_model, img_residual, DataFrame(fits)
end

function semilog2(x)
    return @. log(1 + abs(x))^2
end

function semilog2_der(x)
    return @. 2*log(1 + abs(x)) /(1+abs(x)) * sign(x)
end

function (@main)(ARGS)
    imgname, weightname, flagsname = ARGS

    img_model, img_residual, df = run_all(imgname, weightname, flagsname)

    img_dir = dirname(imgname)
    if img_dir == ""
        img_dir = "."
    end
    write_image(img_dir * "/bright_model.fits", img_model)
    write_image(img_dir * "/bright_residual.fits", img_residual)
    write_fits(img_dir * "/bright_catalogue.fits", df, overwrite=true)
end
