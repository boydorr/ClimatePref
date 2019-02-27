using MyUnitful
using ClimatePref
using Simulation
using Distributions
using Unitful
using Unitful.DefaultSymbols
using AxisArrays
using Profile
function create_eco(numSpecies::Int64, grid::Tuple{Int64, Int64}, area::Unitful.Area{Float64}, totalK::Unitful.Quantity{Float64}, req::Unitful.Quantity{Float64}, individuals::Int64)
    # Set up initial parameters for ecosystem

    # Set up how much energy each species consumes
    energy_vec = SolarRequirement(fill(req, numSpecies))
    #energy_vec = SimpleRequirement([2.0])
    # Set probabilities
    birth = 0.6/year
    death = 0.6/year
    l = 1.0
    s = 0.0
    boost = 1.0

    # Collect model parameters together (in this order!!)
    param = EqualPop(birth, death, l, s, boost)

    #dir = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/wc"
    #dir = "/home/claireh/Documents/RawWorldclim"
    #srad = readworldclim(joinpath(dir, "wc2.0_5m_srad"))
    #srad.array = srad.array[-10° .. 60°, 35° .. 80°]
    #meansrad = mean(srad.array[.!isnan.(srad.array)])
    #totalK = uconvert(kJ, meansrad * month * (area/(grid[1]*grid[2])))/5

    # Create ecosystem
    kernel = GaussianKernel(0.1km, numSpecies, 10e-10)
    movement = BirthOnlyMovement(kernel, Torus())

    opts = fill(274.0K, numSpecies)
    vars = fill(0.5K, numSpecies)
    traits = GaussTrait(opts, vars)
    native = fill(true, numSpecies)
    abun = rand(Multinomial(individuals, numSpecies))
    sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
        movement, param, native)
    abenv = simplehabitatAE(274.0K, grid, totalK, area)
    rel = Gauss{typeof(1.0K)}()
    eco = Ecosystem(sppl,abenv,rel)
end
function runsim(eco::Ecosystem, times::Unitful.Time, reps::Int64)
    burnin = 1month; interval = 1month; timestep = 1month
    lensim = length(0month:interval:times)
    abun = generate_storage(eco, lensim, reps)
    totalK = sum(eco.abenv.budget.matrix)
    grid = size(eco.abenv.habitat.matrix)
    for j in 1:reps
        repopulate!(eco)
        reenergise!(eco, totalK, grid)
        thisstore = view(abun, :, :, :, j)
        simulate_record!(thisstore, eco, times, interval, timestep)
    end
    return abun
end

function plot_abun(abun::Array{Int64, 4}, grd::Tuple{Int64, Int64}, spp::UnitRange{Int64} = 1:size(abun,1), thin::Int64=size(abun, 4))
    pick = rand(1:size(abun, 4), thin)
    if size(abun, 2) >16
        abun = mapslices(sum, abun, dims = 2)
        grd = (1,1)
    end
    for k in 1:size(abun, 2)
        for n in spp
            if k ==1 && n == spp[1]
                display(plot(abun[1,k,:,pick], ylabel = "Abundance", xlabel = "Months", label="", grid = false, color = :($k), linealpha = 0.1, layout = grd, subplot = k,
                ylim = (0, maximum(abun))))
            else
                display(plot!(abun[n,k,:,pick], ylabel = "Abundance", xlabel = "Months", label="", color = :($k), linealpha = 0.1, grid = false, subplot=k, ylim = (0, maximum(abun))))
            end
        end
    end
end


function plot_abun(abun::Array{Int64, 4}, spp::Int64, thin::Int64=size(abun, 4))
    pick = rand(1:size(abun, 4), thin)
    for k in 1:size(abun, 2)
        if k ==1
            display(plot(abun[spp,k,:,pick], ylabel = "Abundance", xlabel = "Months", label="", grid = false, color = :($k), linealpha = 0.1, layout = grd, subplot = k,
            ylim = (0, maximum(abun))))
        else
            display(plot!(abun[spp,k,:,pick], ylabel = "Abundance", xlabel = "Months", label="", color = :($k), linealpha = 0.1, grid = false, subplot=k, ylim = (0, maximum(abun))))
        end
    end
end
