using Distributions
import Distributions: @check_args, ContinuousUnivariateDistribution,
@distr_support, rand, params, GLOBAL_RNG, pdf
struct Trapezoid{T<:Real} <: ContinuousUnivariateDistribution
    a::T
    b::T
    c::T
    d::T

    Trapezoid{T}(a::T, b::T, c::T, d::T) where {T} = (@check_args(Trapezoid, a < d); new{T}(a, b, c, d))
end
Trapezoid(a::T, b::T, c::T, d::T) where {T<:Real} = Trapezoid{T}(a, b, c, d)
Trapezoid(a::Real, b::Real, c::Real, d::Real) = Uniform(promote(a, b, c, d)...)
Trapezoid(a::Integer, b::Integer, c::Integer, d::Integer) = Trapezoid(Float64(a), Float64(b),
    Float64(c), Float64(d))
Trapezoid() = Trapezoid(0.0, 0.0, 1.0, 1.0)

@distr_support Trapezoid d.a d.b d.c d.d
params(d::Trapezoid) = (d.a, d.b, d.c, d.d)
rand(d::Trapezoid) = rand(GLOBAL_RNG, d)
function rand(rng::AbstractRNG, T::Trapezoid)
    (a, b, c, d) = params(T)
    b_m_a = b - a
    c_m_b = c - b
    d_m_c = d - c
    Cϕ = 4 / ((b_m_a * 2) + (c_m_b * 4) + (d_m_c * 2))
    pi1 = Cϕ * b_m_a/2
    pi2 = Cϕ * c_m_b
    pi3 = Cϕ * d_m_c/2
    u = rand(rng)
    if (u >= 0 && u <= pi1)
        return a + (u/pi1)^(1/2) * b_m_a
    elseif (u > pi1 && u <= (1 - pi3))
        return b + ((u - pi1)/pi2) * c_m_b
    else
        return d - ((1 - u)/pi3)^(1/2) * d_m_c
    end
end

function pdf(T::Trapezoid, x::Real)
    (a, b, c, d) = params(T)
    d_p_c = d + c
    a_p_b = a + b
    b_m_a = b - a
    d_m_c = d - c
    if (x >= a && x < b)
        return (2 / (d_p_c- a_p_b)) * ((x-a)/b_m_a)
    elseif (x >= b && x < c)
        return (2 / (d_p_c- a_p_b))
    else
        return (2 / (d_p_c- a_p_b)) * ((d-x)/d_m_c)
    end
end

using JLD
using AxisArrays
tavg = JLD.load("data/Bins/Binned_data_tavg.jld", "Binned_data_tavg")
tmax = JLD.load("data/Bins/Binned_data_tmax.jld", "Binned_data_tmax")
tmin = JLD.load("data/Bins/Binned_data_tmin.jld", "Binned_data_tmin")
distributions = AxisArray(zeros(Int64, size(tavg, 1), 4),
                    Axis{:Species}(tavg.axes[1].val), Axis{:Param}(["a","b","c","d"]))
for i in 1:length(tavg.axes[1].val)
    row = tavg[i,:]
    rowmin = tmin[i,:]
    rowmax = tmax[i,:]
    if any(row .> 0)
        start = minimum(find(row .> 0))[1] - 1
        ending = length(row) - minimum(find(reverse(row) .> 0))[1] + 1
        startmin = minimum(find(rowmin .> 0))[1] - 1
        endmax =  length(row) - minimum(find(reverse(rowmax) .> 0))[1] + 1
        a = tmin.axes[2].val[startmin]
        b = tavg.axes[2].val[start]
        c = tavg.axes[2].val[ending]
        d = tmax.axes[2].val[endmax]
        distributions[i, :] = [a, b, c, d]
    end
end
save("data/Temperature.jld", "Temperature", distributions)

using JLD
using AxisArrays
prec = JLD.load("data/Bins/Binned_data_prec.jld", "Binned_data_prec")
distributions = AxisArray(zeros(Int64, size(prec, 1), 2),
                    Axis{:Species}(prec.axes[1].val), Axis{:Param}(["a","b"]))
for i in 1:length(prec.axes[1].val)
    row = prec[i,:]
    step = prec.axes[2].val[2] - prec.axes[2].val[1]
    if any(row .> 0)
        start = minimum(find(row .> 0))[1]
        ending = length(row) - minimum(find(reverse(row) .> 0))[1] + 1
        a = prec.axes[2].val[start] - step
        b = prec.axes[2].val[ending]
        distributions[i, :] = [a, b]
    end
end
save("data/Rainfall.jld", "Rainfall", distributions)

tavg = JLD.load("data/Bins/Binned_data_monthly_tavg.jld", "Binned_data_monthly_tavg")
tmax = JLD.load("data/Bins/Binned_data_monthly_tmax.jld", "Binned_data_monthly_tmax")
tmin = JLD.load("data/Bins/Binned_data_monthly_tmin.jld", "Binned_data_monthly_tmin")
distributions = AxisArray(zeros(Int64, size(tavg, 1), 4),
                    Axis{:Species}(tavg.axes[1].val), Axis{:Param}(["a","b","c","d"]))
for i in 1:length(tavg.axes[1].val)
    row = tavg[i,:]
    rowmin = tmin[i,:]
    rowmax = tmax[i,:]
    if any(row .> 0)
        start = minimum(find(row .> 0))[1] - 1
        ending = length(row) - minimum(find(reverse(row) .> 0))[1] + 1
        startmin = minimum(find(rowmin .> 0))[1] - 1
        endmax =  length(row) - minimum(find(reverse(rowmax) .> 0))[1] + 1
        a = tmin.axes[2].val[startmin]
        b = tavg.axes[2].val[start]
        c = tavg.axes[2].val[ending]
        d = tmax.axes[2].val[endmax]
        distributions[i, :] = [a, b, c, d]
    end
end
save("data/Temperature.jld", "Temperature", distributions)
