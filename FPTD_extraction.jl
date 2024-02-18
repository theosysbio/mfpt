# This script shows how to extract the full time-dependent first-passage time distribution
# (FPTD) for the telegraph model given in reaction scheme 13 of the main text. 

# SECTION 1: PREAMBLE
# 1a. Import the necessary packages
using ModelingToolkit
using Catalyst
using DifferentialEquations
using SparseArrays
using LinearAlgebra
using Plots
using FiniteStateProjection
using FiniteStateProjection: singleindices, pairedindices, netstoichmat

# SECTION 2: defines the functions necessary to compute the FPTD
# 2a. This functions computes the sparse matrix used to extract the FPTD 
function create_sparsematrix_FPTD(sys, dims::NTuple, ps, t; idx_filter=x -> true)
    Ntot = prod(dims)
    lind = LinearIndices(sys.ih, dims)

    I = Int[]
    J = Int[]
    V = Float64[]

    predsize = Ntot * (length(Catalyst.get_eqs(sys.rs)) + 1)

    sizehint!(I, predsize)
    sizehint!(J, predsize)
    sizehint!(V, predsize)

    for idx_cart in singleindices(sys.ih, dims)
        idx_lin = lind[idx_cart]
        push!(I, idx_lin)
        push!(J, idx_lin)

        if !idx_filter(idx_cart)
            push!(V, -1.0)
            continue
        end

        rate = 0.0
        for rf in sys.rfs
            rate -= rf(idx_cart, t, ps...)
        end

        push!(V, -rate)
    end

    S::Matrix{Int64} = netstoichmat(sys.rs)
    for (i, rf) in enumerate(sys.rfs)
        for (idx_cin, idx_cout) in pairedindices(sys.ih, dims, CartesianIndex(S[:, i]...))

            if !idx_filter(idx_cin) || !idx_filter(idx_cout)
                continue
            end

            idx_lin = lind[idx_cin]
            idx_lout = lind[idx_cout]
            push!(I, idx_lout)
            push!(J, idx_lin)
            rate = rf(idx_cin, t, ps...)
            push!(V, -rate)
        end
    end

    sparse(I, J, V)
end;

# 2b. Override Base.convert from FiniteStateProjection.jl
function Base.convert(::Type{SparseMatrixCSC}, sys::FSPSystem, dims::NTuple, ps, t::Real)
    create_sparsematrix_FPTD(sys, dims, ps, t)
end;

# SECTION 3: THE TELEGRAPH MODEL
# parameters & variables are defined in Equation 13 of the main text
@parameters λ, μ, K, δ
@variables t, G, M
model = @reaction_network begin
    λ * (1 - G), 0 --> G
    μ, G --> 0
    K, G --> G + M
    δ, M --> 0
end

# assign parameter values to the model
# by changing the value of λ one can reproduce all three
# ridgeline plots in the Supplementary material e.g. λ = 
# (10.0, 1.0, 0.1)
# We fix the steady-state mean of the system and vary K
# accordingly
λ, μ, δ = 10.0, 0.5, 1.0
ss_mean = 100
K = ss_mean * ((λ + μ) / λ)
params = (λ, μ, K, δ)

# define what fraction of the model's steady-state mean you are 
# interested in 
frac_ss_mean = 0.8

# compute 80% of the steady-state mean
target = ss_mean * frac_ss_mean

# SECTION 4: THE FPTD OF THE TELEGRAPH MODEL
# compute the (sparse) matrix A
# initial distribution (over two species)
u0 = ones(2, Int(target))
# here we start with 0 copies of G_on and M
u0[1, 1] = 1
sys = FSPSystem(model);
A = convert(SparseMatrixCSC, sys, size(u0), params, 0)

# define the function f for use with ODE Problem ito the matrix (A__Y)^T 
dt = 0.001
f = (du, u, p, t) -> mul!(vec(du), -transpose(A), vec(u))
prob = ODEProblem(f, u0, 100, params)
sol = solve(prob, Rosenbrock23(), adaptive=false, dt=dt)

# extract the FPTD from the solution to the ODE Problem
deriv = sol(sol.t, Val{1})
pdf = .-deriv[1, :]
AUC = sum(pdf) * dt

# plot the FPTD from the solution
plot(sol.t, pdf, label="FPTD")

# END