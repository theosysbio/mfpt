# This script produces the data given in sub-section "3.2 The telegraph
# model" of the main text.



# SECTION 1: PREAMBLE

# 1a. Import the necessary packages
using ModelingToolkit
using Catalyst
using DifferentialEquations
using SparseArrays
using FiniteStateProjection
using FiniteStateProjection: singleindices, pairedindices, netstoichmat

# 1b. Import "mfpt_functions.jl" script as follows if they are saved in the 
# same directory, otherwise amend the path accordingly
include("mfpt_functions.jl")



# SECTION 2: THE TELEGRAPH MODEL

# parameters & variables are defined in Equation 13 of the main text
@parameters λ, μ, K, δ
@variables t, G, M
model = @reaction_network begin
    λ * (1 - G), 0 --> G
    μ, G --> 0
    K, G --> G + M
    δ, M --> 0
end λ μ K δ

# assign parameter values to the model
# by changing the value of λ one can reproduce all three
# plots in Figure 3(c) e.g. λ = (10.0, 1.0, 0.1)
λ, μ, K, δ = 1.0, 20.0, 1000.0, 1.0
params = (λ, μ, K, δ)

# define what fraction of the model's steady-state mean you are 
# interested in 
frac_ss_mean = 0.8



# SECTION 3: COMPUTING THE DATA IN FIGURE 3

# compute the average time taken to reach some percentage (frac_ss_mean)
# of the steady-state mean in both the stochastic and deterministic regimes
# for increasing values of the parameter μ

μs = collect(0.0:0.1:100.0)

det_wts = []
sto_wts = []

# set the steady-state mean of the system
ss_mean = 100

for μ in μs

    params = (λ, μ, K, δ)
    # adjust K such that the steady-state mean stays constant for different values of μ
    K = ss_mean * ((params[1] + params[2]) / params[1])

    params = (λ, μ, K, δ)

    # compute 80% of the steady-state mean
    target = ss_mean * frac_ss_mean

    # compute the deterministic FPT
    firsttimeprob_det, firsttimecb_det = let
        # Initial conditions
        setdefaults!(model, [:G => 0, :M => 0, :λ => params[1], :μ => params[2], :K => params[3], :δ => params[4]])

        # Convert to an ODESystem
        odesys = convert(ODESystem, model)

        # Find its integer index within the solution, u
        Sidx = 2

        # Condition function to stop the simulation
        function MFPT_S_det(u, t, integrator, Sidx = Sidx)
            u[Sidx] >= target
        end

        # An affect! to stop the simulation when the condition is true
        function stop_sim(integrator)
            savevalues!(integrator, true)
            terminate!(integrator)
            nothing
        end

        # Create the callback
        cb_det = DiscreteCallback(MFPT_S_det, stop_sim)

        # Define the ODEProblem
        oprob = ODEProblem(model, [], (0.0, 500.0), params)

        oprob, cb_det
    end

    # solve the ODEProblem and terminate when the callback is reached i.e.
    # when 80% of the steady-state mean is achieved
    sol = solve(firsttimeprob_det, Rosenbrock23(), adaptive=false, dt=0.001, save_end=false, callback=firsttimecb_det)
    det_wt = sol.t[end]

    # compute the (stochcastic) MFPT using the functions in "mfpt_functions.jl"
    u0 = zeros(size(model.states)[1], Int((round(ss_mean * frac_ss_mean, digits=0))))
    sto_wt = fsp_time(model, params, Int((round(ss_mean * frac_ss_mean, digits=0))), 0, u0)[1]
   
    push!(det_wts, det_wt)
    push!(sto_wts, sto_wt)
end



# SECTION 4: PLOTTING THE DATA

using LaTeXStrings
using Plots
using Plots.PlotMeasures

msw = 0
ms = 5
p_det_wts = plot(log10.(μs), det_wts; seriestype=:scatter, label="Deterministic",
    xlabel=L"log(μ)", ylabel="Waiting time", color=:darkolivegreen, 
    markerstrokewidth=msw, markersize=ms)
# NOTE: depending on the value of λ = (0.1, 1.0, 10.0) you can change the argument 
# "ylims" to best represent the differences in the data
p_stoc_wts = plot!(log10.(μs), sto_wts; seriestype=:scatter, label="Stochastic",
    framestyle=:box, color=:coral1, markerstrokewidth=msw, markersize=ms, ylims=(1.5,3.0))