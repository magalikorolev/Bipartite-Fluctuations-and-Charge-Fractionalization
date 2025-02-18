using ITensors, ITensorMPS
#using Plots
#using LsqFit
#using Statistics
using Serialization

arg1 = parse(Float64, ARGS[1])

# Parameters of the system
t = 1
U = 0
N = 300
delta = arg1

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
        add!(h, delta*tanh(i-(L/2)), "Sz", i)
    end
    return MPO(h, sites), sites
end

function pos(x, sites)
    op = OpSum()
    add!(op, "Sz", x)
    return MPO(op,sites)
end

# Perform DMRG to find the ground state and calculate fluctuations

# Create the Hamiltonian and the sites and get the GS
H, sites = tight_binding_chain(N)
psi = randomMPS(sites)
sweeps = Sweeps(25)
maxdim!(sweeps, 10, 20, 80, 100, 200)
cutoff!(sweeps, 1e-10)
energy, psi0 = dmrg(H, psi, sweeps)

serialize("results_psi0/nint_tanh_delta_$delta.dat", psi0)
#psi1 = deserialize("psi0_delta_$delta.dat")

#prob_density = [real(norm(psi1[i])).^2 ./N for i in 1:N]

#p1 = [1.0,1.0,1.0]
#function expfit(x,k)
#    return k[1] .* exp.( .- abs.(x .- Int(N/2)) ./ k[2]) .+ k[3]
#end

#fit_exp = curve_fit(expfit,1:N,prob_density,p1)
#a,b,c = coef(fit_exp)
#a_round = round(a,digits=4)
#b_round = round(b,digits=4)
#c_round = round(c,digits=4)
#fit_y_exp = expfit(1:N,coef(fit_exp))

#plot(1:N, prob_density, xlabel="Site", ylabel="|ψ(x)|^2", legend=false)
#plot!(1:N,fit_y_exp,linestyle=:dash)


