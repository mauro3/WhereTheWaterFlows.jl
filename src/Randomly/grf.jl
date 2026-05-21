# Sampling Gaussian Random fields
#
# Follows Raess et al 2019, https://doi.org/10.1016/j.cageo.2019.06.007

using FFTW, Primes

export make_grf_sampler

"""
    find_next_2357(n::Integer) -> Integer

Finds the next integer following `n` that has only 2, 3, 5, and 7 as its prime factors.
This is useful for optimization in FFT operations, which are more efficient on sizes factorable into these primes.
"""
function find_next_2357(n::Integer)
    while true
        for f in reverse(collect(factor(n)))
            f[1] in (2, 3, 5, 7) || break
            f[1]>7 && break
            return n
        end
        n += 1
    end
end

"""
    pad(nx, ny, len) -> Tuple

Calculates the padding needed on both sides to facilitate efficient and non-periodic FFT processing.

# Arguments
- `nx, ny`::Int`: Original x and y dimension.
- `len::Number`: correlation length (in number of grid cells)

# Returns
- `Tuple`: A tuple containing the new x dimension, new y dimension, and the padding
           to be applied on each side to reach these new dimensions.
"""
function pad(nx, ny, len)
    # padding by a factor 2 on both sides is what Raess says
    padfac = 2
    padding = ((find_next_2357(2 * padfac * ceil(Int, len) + nx) - nx)÷2,
               (find_next_2357(2 * padfac * ceil(Int, len) + ny) - ny)÷2)
    nnx, nny = (nx, ny) .+ 2 .* padding
    return nnx, nny, padding
end

"""
    depad(out::Array, nx::Int, ny::Int, padding::Tuple) -> Array

Removes padding from an array.
"""
depad(out, nx, ny, padding) = out[1+padding[1]:nx+padding[1], 1+padding[2]:ny+padding[2]]

"""
    distance_squared(nx, ny) -> Array

Computes the square of the Euclidean distance from the center in a 2D grid.
"""
function distance_squared(nx, ny)
    xs = (1:nx) .- nx÷2
    ys = ((1:ny) .- ny÷2)'
    xs .^ 2 .+ ys .^ 2
end

"""
    gaussian_kernel(nx::Int, ny::Int, len::Number) -> Array

Generates a Gaussian kernel based on the input dimensions and correlation length.
"""
function gaussian_kernel(nx, ny, len)
    D = distance_squared(nx, ny)
    return exp.(-D / len^2)
end

"""
    exponential_kernel(nx::Int, ny::Int, len::Number) -> Array

Generates an exponential kernel based on the input dimensions and length scale.
"""
function exponential_kernel(nx, ny, len)
    D = sqrt.(distance_squared(nx, ny))
    return exp.(-D / len)
end

"""
    make_kernel(kernel_fn::Function, nx::Int, ny::Int, len::Number) -> Array

Precomputes and returns the Fourier transform of a kernel function over specified dimensions
and correlation length.
"""
function make_kernel(kernel_fn, nx, ny, len)
    nnx, nny = pad(nx, ny, len)
    sqrt.(rfft(kernel_fn(nnx, nny, len)))
end

"""
    eq18(nx::Int, ny::Int, len::Number, kernel::Array) -> Array

Draws a sample of the Gaussian random field using the precomputed kernel.
"""
function eq18(nx, ny, len, kernel)
    nnx, nny, padding = pad(nx, ny, len)
    Z = randn(typeof(len), nnx, nny)
    tmp = rfft(Z)
    tmp .*= kernel
    return depad(irfft(tmp, nnx), nx, ny, padding)
end

# Same as above but using planned FFTW for better performance
function eq18_plannedfftw(nx, ny, len, kernel, p_rfft, p_irfft)
    nnx, nny, padding = pad(nx, ny, len)
    Z = randn(typeof(len), nnx, nny)
    tmp = p_rfft * Z
    tmp .*= kernel
    return depad(p_irfft * tmp, nx, ny, padding)
end

"""
    make_grf_sampler(nx::Int, ny::Int, kernel_fn::Function, len::Number) -> Function

Creates a sampler function for generating Gaussian random fields (GRF) using pre-specified parameters.

# KW-arguments:
- `fftw_plan`: if not `nothing` then FFTW is planned for improved performance.
               Use `:estimate` or `:measure` (recommended) to use a plan.
"""
function make_grf_sampler(nx, ny, kernel_fn, len; T=Float64, fftw_plan=nothing)
    if min(nx÷3, ny÷3)<len
        @warn "Correlation length `len` should be at least 3x smaller than the domain length"
    end
    if fftw_plan==nothing
        let
            l = max(T(len), eps(T)) # l==0 leads to NaNs
            ker = make_kernel(kernel_fn, nx, ny, l) # pre-compute kernel
            return () -> eq18(nx, ny, l, ker)
        end
    elseif fftw_plan in [:estimate, :measure]
        let
            l = max(T(len), eps(T)) # l==0 leads to NaNs
            ker = make_kernel(kernel_fn, nx, ny, l) # pre-compute kernel
            flags = fftw_plan==:estimate ? FFTW.ESTIMATE : FFTW.MEASURE
            nnx, nny, _ = pad(nx, ny, l)
            # plan it
            tmp = randn(nnx, nny)
            p_rfft = plan_rfft(tmp; flags)
            p_irfft = plan_irfft(p_rfft*tmp, nnx; flags)
            return () -> eq18_plannedfftw(nx, ny, l, ker, p_rfft, p_irfft)
        end
    else
        error("Unrecognized fftw_plan value: $fftw_plan")
    end
end
