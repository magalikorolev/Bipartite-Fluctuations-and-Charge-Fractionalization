using ITensors, ITensorMPS

arg1 = parse(Float64, ARGS[1])

# Parameters of the system
t = 1
U = 0
N = 300
mu = arg1

# Get Hamiltonian and sites
function tight_binding_chain(L)
    sites = siteinds("S=1/2", L)
    h = OpSum()
    for i in 1:L-1
        add!(h, -2*t, "Sx", i, "Sx", i+1)
        add!(h, -2*t, "Sy", i, "Sy", i+1)
        add!(h, U, "Sz", i, "Sz", i+1)
    end
    add!(h, -2*t, "Sx", L, "Sx", 1) #PBC
    add!(h, -2*t, "Sy", L, "Sy", 1) #PBC
    add!(h, U, "Sz", L, "Sz", 1) #PBC
    for i in 1:L
        add!(h, mu, "Sz", i)
    end
    return MPO(h, sites), sites
end


# Get operator θ-ϕ and θ+ϕ
function operator_theta_R(L, sites)
    opOR = OpSum()
    for i in 2:L
        add!(opOR, 1, "Sx", i, "Sy", i-1)
        add!(opOR, -1, "Sy", i, "Sx", i-1)
        add!(opOR, -1, "Sz", i)
    end
    return MPO(opOR,sites)
end

function operator_theta_L(L, sites)
    opOL = OpSum()
    for i in 2:L
        add!(opOL, 1, "Sx", i, "Sy", i-1)
        add!(opOL, -1, "Sy", i, "Sx", i-1)
        add!(opOL, 1, "Sz", i)
    end
    return MPO(opOL,sites)
end

# Get operator (θ-ϕ)^2 and (θ+ϕ)^2
function operator_theta_R2(L, sites)
    opOR2 = OpSum()
    for i in 2:L
        for j in 2:L
            add!(opOR2, 1, "Sx", i, "Sy", i-1, "Sx", j, "Sy", j-1)
            add!(opOR2, -1, "Sx", i, "Sy", i-1, "Sy", j, "Sx", j-1)
            add!(opOR2, -1, "Sy", i, "Sx", i-1, "Sx", j, "Sy", j-1)
            add!(opOR2, 1, "Sy", i, "Sx", i-1, "Sy", j, "Sx", j-1)
            add!(opOR2, -1, "Sz", i, "Sx", j, "Sy", j-1)
            add!(opOR2, 1,"Sz", i, "Sy", j, "Sx", j-1)
            add!(opOR2, 1, "Sz", i, "Sz", j)
        end
    end
    return MPO(opOR2,sites)
end

function operator_theta_L2(L, sites)
    opOL2 = OpSum()
    for i in 2:L
        for j in 2:L
            add!(opOL2, 1, "Sx", i, "Sy", i-1, "Sx", j, "Sy", j-1)
            add!(opOL2, -1, "Sx", i, "Sy", i-1, "Sy", j, "Sx", j-1)
            add!(opOL2, -1, "Sy", i, "Sx", i-1, "Sx", j, "Sy", j-1)
            add!(opOL2, 1, "Sy", i, "Sx", i-1, "Sy", j, "Sx", j-1)
            add!(opOL2, 1,"Sz", i, "Sx", j, "Sy", j-1)
            add!(opOL2, -1, "Sz", i, "Sy", j, "Sx", j-1)
            add!(opOL2, -1, "Sz", i, "Sz", j)
        end
    end
    return MPO(opOL2,sites)
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


# Perform DMRG to find the ground state and calculate fluctuations

# Create the Hamiltonian and the sites and get the GS
H, sites = tight_binding_chain(N)
psi = randomMPS(sites)
sweeps = Sweeps(25)
maxdim!(sweeps, 10, 20, 80, 100, 200, 250)
cutoff!(sweeps, 1e-10)
energy, psi0 = dmrg(H, psi, sweeps)

function calculate_fluctuations(L)

    # Calculate <theta>
    theta1 = operator_theta_R(L, sites)
    avgtheta1 = inner(psi0', theta1, psi0)

    # Calculate <theta^2>
    theta2 = operator_theta_R2(L, sites)
    avgtheta2 = inner(psi0', theta2, psi0)

    O1th = operator_theta(L, sites)
    avgOth = inner(psi0', O1th, psi0)
    O2th = operator_theta2(L, sites)
    avgO2th = inner(psi0', O2th, psi0)

    O1ph = operator_phi(L, sites)
    avgOph = inner(psi0', O1ph, psi0)
    O2ph = operator_phi2(L, sites)
    avgO2ph = inner(psi0', O2ph, psi0)

    # Calculate fluctuations
    flc_chir = avgtheta2 - avgtheta1^2 
    flc_th = avgO2th - avgOth^2
    flc_ph = avgO2ph - avgOph^2

    println("system size =", L)

    return flc_chir,flc_th,flc_ph
end 

# Calculate the average value of O for different system sizes L
System_sizes = 2:2:300  # Change this range to explore more values of L
fluct = [calculate_fluctuations(x) for x in System_sizes]
#fluct_re = pi^2*real(fluct)

open("results_mu/flc_nint_mu_$arg1.txt", "w") do file
    for i in 1:length(System_sizes)
        println(file, "$(System_sizes[i]),$(fluct[i])")
    end
end 
