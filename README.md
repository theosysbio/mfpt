# Mean first passage times for stochastic biochemical reaction networks
A framework to compute the MFPT of inherently stochastic biochemical reaction networks. Traditional methods to estimate the MFPT rely on numerical integration of deterministic models, however, these ignore intrinsic noise and hence their predictions may be inaccurate. Here we provide an accurate and efficient computational framework to compute the MFPT for reaction networks in the presence of (intrinsic) noise.

# Description
To compute MFPTs numerically, we use an adaptation of [FiniteStateProjection.jl](https://github.com/kaandocal/FiniteStateProjection.jl) to construct the transition matrix of a stochastic reaction network and solve a modified system of equations using the standard sparse solvers in Julia. For more details see the supplementary material for [The timing of cellular events: a stochastic vs deterministic perspective](https://www.biorxiv.org/content/10.1101/2023.07.20.549956v1), specifically the section "Finite State Projection for the modified CME". 

# Examples
The script [mfpt_functions.jl](https://github.com/theosysbio/mfpt/blob/main/mfpt_functions.jl) contains the core functions to compute the MFPT for any chemical reaction network and the script [telegraph_example.jl](https://github.com/theosysbio/mfpt/blob/main/telegraph_example.jl) implements a working example of computing MFPTs for a telegraph process. The script [FPTD_extraction.jl](https://github.com/theosysbio/mfpt/blob/main/FPTD_extractionn.jl) shows how to obtain the full time-dependant FPTD using the telegraph model as an example.

# Literature
For more details on the method and additional examples please refer to the associated paper: [The timing of cellular events: a stochastic vs deterministic perspective](https://www.biorxiv.org/content/10.1101/2023.07.20.549956v1)
