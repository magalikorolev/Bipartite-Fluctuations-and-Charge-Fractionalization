using ITensors, ITensorMPS
#using Plots

#BLAS.set_num_threads(3)

arg1 = parse(Float64, ARGS[1])

# Parameters of the system
t1 = 1
t2 = 1
V = 0 
U = 0.2 
tperp = 1.5 
N = 400
mu = arg1

# Get Hamiltonian and sites
function hamiltonian(L)
    sites = siteinds("S=1/2", L)
    h = OpSum()
    for i in 1:2:L-3
        # Intra-chain
        add!(h, -2*t1, "Sx", i, "Sx", i+2)
        add!(h, -2*t1, "Sy", i, "Sy", i+2)
        add!(h, V, "Sz", i, "Sz", i+2)
        add!(h, -2*t2, "Sx", i+1, "Sx", i+3)
        add!(h, -2*t2, "Sy", i+1, "Sy", i+3)
        add!(h, V, "Sz", i+1, "Sz", i+3)
        # Inter-chains
        add!(h, U, "Sz", i, "Sz", i+1)
        add!(h, -2*tperp, "Sx", i, "Sx", i+1)
        add!(h, -2*tperp, "Sy", i, "Sy", i+1)
        # Add the chemical potential
        add!(h, mu, "Sz", i)
        add!(h, mu, "Sz", i+1)
    end
    # PBC
    add!(h, -2*t1, "Sx", L-1, "Sx", 1) 
    add!(h, -2*t1, "Sy", L-1, "Sy", 1) 
    add!(h, V, "Sz", L-1, "Sz", 1) 
    add!(h, -2*t2, "Sx", L, "Sx", 2)
    add!(h, -2*t2, "Sy", L, "Sy", 2) 
    add!(h, V, "Sz", L, "Sz", 2)
    add!(h, U, "Sz", L-1, "Sz", L)
    add!(h, -2*tperp, "Sx", L-1, "Sx", L)
    add!(h, -2*tperp, "Sy", L-1, "Sy", L)
    add!(h, mu, "Sz", L-1)
    add!(h, mu, "Sz", L)
    return MPO(h, sites), sites
end

# Operators O
function operator_O_1(L, sites)
    opO_1 = OpSum()
    for i in 3:2:L-1
        add!(opO_1, 1, "Sz", i)
    end
    return MPO(opO_1, sites)
end

function operator_O_2(L, sites)
    opO_2 = OpSum()
    for i in 2:2:L
        add!(opO_2, 1, "Sz", i)
    end
    return MPO(opO_2, sites)
end

# Operators O2
function operator_O2_1(L, sites)
    opO2_1 = OpSum()
    for i in 3:2:L-1
        for j in 3:2:L-1
            add!(opO2_1, 1, "Sz", i, "Sz", j)
        end
    end
    return MPO(opO2_1, sites)
end

function operator_O2_2(L, sites)
    opO2_2 = OpSum()
    for i in 2:2:L
        for j in 2:2:L
            add!(opO2_2, 1, "Sz", i, "Sz", j)
        end
    end
    return MPO(opO2_2, sites)
end

function operator_O2_3(L, sites)
    opO2_3 = OpSum()
    for i in 3:2:L-1
        for j in 2:2:L
            add!(opO2_3, 1, "Sz", i, "Sz", j)
        end
    end
    return MPO(opO2_3, sites)
end



function operator_O_bis(L, sites)
    opO = OpSum()
    for i in 1:L
        add!(opO, 1/2, "Sz", i)
    end
    return MPO(opO,sites)
end

function operator_O2_bis(L, sites)
    opO2 = OpSum()
    for i in 1:L-1
        for j in 1:L-1
            add!(opO2, 1/4, "Sz", i, "Sz", j)
            add!(opO2, 1/4, "Sz", i, "Sz", j+1)
            add!(opO2, 1/4, "Sz", i+1, "Sz", j)
            add!(opO2, 1/4, "Sz", i+1, "Sz", j+1)
        end
    end
    return MPO(opO2, sites)
end

# Perform DMRG to find the ground state and calculate fluctuations
H, sites = hamiltonian(N)
psi = randomMPS(sites)
sweeps = Sweeps(25)
maxdim!(sweeps, 10, 20, 40, 80, 100, 200)
cutoff!(sweeps, 1e-10)
energy, psi0 = dmrg(H, psi, sweeps)

function calculate_fluctuations(L)
    # Calculate <O>
    #println("Test O1")
    #O_1 = operator_O_1(L, sites)
    #println("Test avg O1")
    #avgO_1 = inner(psi0', O_1, psi0)
    #println("Test O2")
    #O_2 = operator_O_2(L, sites)
    #avgO_2 = inner(psi0', O_2, psi0)

    # Calculate <O^2>
    #println("Test O1^2")
    #O2_1 = operator_O2_1(L, sites)
    #avgO2_1 = inner(psi0', O2_1, psi0)
    #println("Test O2^2")
    #O2_2 = operator_O2_2(L, sites)
    #avgO2_2 = inner(psi0', O2_2, psi0)
    #println("Test O3^2")
    #O2_3 = operator_O2_3(L, sites)
    #avgO2_3 = inner(psi0', O2_3, psi0)

    O_bis = operator_O_bis(L, sites)
    avgO_bis = inner(psi0', O_bis, psi0)
    O2_bis = operator_O2_bis(L, sites)
    avgO2_bis = inner(psi0', O2_bis, psi0)


    # Calculate fluctuations
    #f1 = avgO2_1 - avgO_1^2
    #f2 = avgO2_2 - avgO_2^2
    fluctuations = avgO2_bis - avgO_bis^2
    println("system size =", L)
    return fluctuations
end 

System_sizes = 4:2:N
fluct = [calculate_fluctuations(x) for x in System_sizes]
fluct_re = pi^2*real(fluct)


#fluct_phi_1 = []
#fluct_phi_2 = []
#fluct_phi_plus = []
#for x in System_sizes
    #c = calculate_fluctuations(x)
    #push!(fluct_phi_1, a)
    #push!(fluct_phi_2, b)
    #push!(fluct_phi_plus, c)
#end



#fluct_re_1 = pi^2*real(fluct_phi_1)
#fluct_re_2 = pi^2*real(fluct_phi_2)
#fluct_re_plus = pi^2*real(fluct_phi_plus)

# Calculate the average value of O for different system sizes L

#fluct = [calculate_fluctuations(x, sites, psi0)[3] for x in System_sizes]
#fluct_re = pi^2*real(fluct)/4

#plot(System_sizes, fluct_re_1, xlabel="x", label="F_ϕ1(x)", title="Fluctuations of ϕ+ N=$N, t1=$t1, t2=$t2, U=$U, V=$V, 20 sweeps")
#plot!(System_sizes, fluct_re_2, xlabel="x", label="F_ϕ2(x)")
#plot!(System_sizes, fluct_re_plus, xlabel="x", label="F_ϕ+(x)")


# Open a file in write mode
open("results_two_chains/mu/Ut02_mu_$arg1.txt", "w") do file
    for i in 1:length(System_sizes)
        println(file, "$(System_sizes[i]),$(fluct_re[i])")
   end
end
