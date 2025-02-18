using ITensors, ITensorMPS

arg1 = parse(Float64, ARGS[1])

# Parameters of the system
t = 1
U = arg1
N = 400
delta = 1.0

# Get Hamiltonian and sites
function tight_binding_chain(L)
    sites = siteinds("S=1/2", L)
    h = OpSum()
    for i in 1:L-1
        add!(h, -2*t, "Sx", i, "Sx", i+1)
        add!(h, -2*t, "Sy", i, "Sy", i+1)
        add!(h, U, "Sz", i, "Sz", i+1)
    end
    #add!(h, -2*t, "Sx", L, "Sx", 1) #PBC
    #add!(h, -2*t, "Sy", L, "Sy", 1) #PBC
    #add!(h, U, "Sz", L, "Sz", 1) #PBC
    for i in 1:L
        add!(h, delta*tanh(i-Int(L/2)), "Sz", i)
    end
    return MPO(h, sites), sites
end


# Get operator phi
function operator_phi(L, sites)
    opO = OpSum()
    for i in 2:L
        add!(opO, 1, "Sz", i)
    end
    return MPO(opO,sites)
end

# Get operator phi2
function operator_phi2(L, sites)
    opO2 = OpSum()
    for i in 2:L
        for j in 2:L
            add!(opO2, 1, "Sz", i, "Sz", j)
        end
    end
    return MPO(opO2,sites)
end


# Get operator theta
function operator_theta(L, sites)
    opO = OpSum()
    for i in 2:L
        add!(opO, 1, "Sx", i, "Sy", i-1)
        add!(opO, -1, "Sy", i, "Sx", i-1)
    end
    return MPO(opO,sites)
end

# Get operator theta2
function operator_theta2(L, sites)
    opO2 = OpSum()
    for i in 2:L
        for j in 2:L
            add!(opO2, 1, "Sx", i, "Sy", i-1, "Sx", j, "Sy", j-1)
            add!(opO2, -1, "Sx", i, "Sy", i-1, "Sy", j, "Sx", j-1)
            add!(opO2, -1, "Sy", i, "Sx", i-1, "Sx", j, "Sy", j-1)
            add!(opO2, 1, "Sy", i, "Sx", i-1, "Sy", j, "Sx", j-1)
        end
    end
    return MPO(opO2,sites)
end

function dens(L, sites)
    op = OpSum()
    add!(op, 1, "Sz",L)
    return MPO(op, sites)
end 

# Perform DMRG to find the ground state and calculate fluctuations

# Create the Hamiltonian and the sites and get the GS
H, sites = tight_binding_chain(N)
psi = randomMPS(sites)
sweeps = Sweeps(30)
maxdim!(sweeps, 10, 20, 80, 100, 150, 250)
cutoff!(sweeps, 1e-10)
energy, psi0 = dmrg(H, psi, sweeps)

System_sizes = 2:2:N  # Change this range to explore more values of L

norm_squared = [real(norm(psi0[i])).^2 ./N for i in 1:N]

open("res_JR_400/norm_squared_delta_1_U_$arg1.txt", "w") do file
    for i in 1:N
        println(file, "$(i),$(norm_squared[i])")
    end
end


function calculate_fluctuations(L)

    th1 = operator_theta(L, sites)
    avgth1 = inner(psi0', th1, psi0)
    ph1 = operator_phi(L, sites)
    avgph1 = inner(psi0', ph1, psi0)

    th2 = operator_theta2(L, sites)
    avgth2 = inner(psi0', th2, psi0)
    ph2 = operator_phi2(L, sites)
    avgph2 = inner(psi0', ph2, psi0)

    # Calculate fluctuations
    flc_th = avgth2 - avgth1^2
    flc_ph = avgph2 - avgph1^2

    pbdens = inner(psi0', dens(L,sites), psi0)

    println("system size =", L)

    return flc_th, flc_ph, pbdens
end 

fluct_th = []
fluct_ph = []
proba_dens = []
for i in System_sizes
    f2, f3, pbd = calculate_fluctuations(i)
    push!(fluct_th, pi^2 * real(f2))
    push!(fluct_ph, pi^2 * real(f3))
    push!(proba_dens, pbd)
end

open("res_JR_400/flc_theta_delta_1_U_$arg1.txt", "w") do file
    for i in 1:length(System_sizes)
        println(file, "$(System_sizes[i]),$(fluct_th[i])")
    end
end

open("res_JR_400/flc_phi_delta_1_U_$arg1.txt", "w") do file
    for i in 1:length(System_sizes)
        println(file, "$(System_sizes[i]),$(fluct_ph[i])")
    end
end

open("res_JR_400/prob_dens_delta_1_U_$arg1.txt", "w") do file
    for i in 1:length(System_sizes)
        println(file, "$(System_sizes[i]),$(proba_dens[i])")
    end
end
