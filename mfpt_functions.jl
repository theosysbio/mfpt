# This script contains the functions necessary to reproduce the data given in "3. Stochstic vs. 
# deterministic calculations" of the main text. We include an example of how to use these 
# functions in the script called "telegraph_example.jl" which produces the data given in
# sub-section "3.2 The telegraph model" of the main text.



# Section 1 defines the functions necessary to compute the moments of the FPT distribution given
# in "2. Mathematical framework" of the main text

# 1a. This functions computes the (sparse) matrix (A__Y)^T given in Equation 6 of the main text 
function create_sparsematrix(sys, dims::NTuple, ps, t; idx_filter=x -> true)
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

    I, J, V
end;


# 1b. This function solves the matrix recurrence relation given in Equation 7 of the main text
function fsp_time(model, p, target, m, u0)
    sys = FSPSystem(model)

    charfunc_Y(idx) = idx[1] - 1 == target

    A_I, A_J, A_V = create_sparsematrix(sys, size(u0), p, 0.0; idx_filter=!charfunc_Y)
    AYT = SparseArrays.sparse(A_J, A_I, A_V)
    sol_up = AYT \ ones(length(u0))
    # sol_up = AYT \ ones(size(u0))
end;