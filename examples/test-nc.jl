using NetCDF
using Unitful
using Unitful.DefaultSymbols
using AxisArrays
using myunitful

ncinfo("/Users/claireh/Downloads/test.nc")

x = ncread("/Users/claireh/Downloads/test.nc", "t2m")

y = x * 1.0

y[y .≈ -32767] = NaN
temps = y .* 0.0009362609303530759K .+ 270.9726153656286K

ncinfo("data/era_interim_moda_1990")
twomtemp = ncread("data/era_interim_moda_1990", "t2m")
twomtemp = twomtemp * 1.0

twomtemp[twomtemp .≈ -32767] = NaN
temps = twomtemp .* 0.0017312391138308897K .+ 270.9726153656286K
tempsC = uconvert.(°C, temps)
tempsC=ustrip.(tempsC)
step1 = 0.75°
step2 = 180°/241
tempax = AxisArray(uconvert.(°C, temps),
                       Axis{:latitude}(-180.0°:step1:(180.0° - step1 / 2)),
                       Axis{:longitude}(-90.0°:step2:(90.0°-step2/2)),
                       Axis{:time}(1month:1month:10year))
using RCall
@rput tempsC
step1u = ustrip(step1)
step2u = ustrip(step2)
@rput step1u; @rput step2u
R"pdf(file='plots/testERAint.pdf', paper = 'a4r', height= 8.27, width=11.69 )
image.plot(seq(-180, 180, by=step1u),
seq(-90, 90, by=step2u), tempsC[,,1], xlab='', ylab='')
dev.off()"
