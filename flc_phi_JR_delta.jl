using ITensors, ITensorMPS
#using Plots
#using LaTeXStrings
#using GLM
#using DataFrames

#BLAS.set_num_threads(3)

arg1 = parse(Float64, ARGS[1])

# Parameters of the system
t = 1
U = 0
N = 300
D = arg1


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
    for i in 1:Int(L/2)
        add!(h, -D, "Sz", i)
        add!(h, D, "Sz", i+Int(L/2))
    end
    return MPO(h, sites), sites
end

# Get operator O
function operator_O(L, sites)
    opO = OpSum()
    for i in 2:L
        add!(opO, 1, "Sz", i)
    end
    return MPO(opO,sites)
end

# Get operator O2
function operator_O2(L, sites)
    opO2 = OpSum()
    for i in 2:L
        for j in 2:L
            add!(opO2, 1, "Sz", i, "Sz", j)
        end
    end
    return MPO(opO2,sites)
end

# Perform DMRG to find the ground state and calculate fluctuations

# Create the Hamiltonian and the sites and get the GS
H, sites = tight_binding_chain(N)
psi = randomMPS(sites)
sweeps = Sweeps(25)
maxdim!(sweeps, 10, 20, 40, 80, 100, 200, 250)
cutoff!(sweeps, 1e-10)
energy, psi0 = dmrg(H, psi, sweeps)

function calculate_fluctuations(L)

    # Calculate <O>
    O1 = operator_O(L, sites)
    avgO = inner(psi0', O1, psi0)

    # Calculate <O^2>
    O2 = operator_O2(L, sites)
    avgO2 = inner(psi0', O2, psi0)

    # Calculate fluctuations
    fluctuations = avgO2 - avgO^2 

    println("system size =", L)

    return fluctuations
end 

# Calculate the average value of O for different system sizes L
System_sizes = 2:2:300  # Change this range to explore more values of L
fluct = [calculate_fluctuations(x) for x in System_sizes]
fluct_re = pi^2*real(fluct)

# Calculate log(L)
#log_L = log.((N/pi)*sin.(pi*System_sizes/N))

# Plot the average value of O as a function of L and log(L)
#plot(System_sizes, fluct_re, xlabel="L", label="F_θ(L)", title="Fluctuations of ϕ and log(L) as a function of L")
#plot!(System_sizes, log_L, label=L"\log(L)", linestyle=:dash)


# Create a DataFrame for the linear regression
#data = DataFrame(logL=log_L, F=fluct_re)
#model = lm(@formula(F ~ logL), data)

# Extract the fit parameters
#intercept = coef(model)[1]
#slope = coef(model)[2]

#println(intercept)
#println(slope)
#intercept_round = round(intercept,digits=3)
#slope_round = round(slope,digits=3)

# Generate the fitted values
#fit_F = intercept .+ slope .* log_L

# Plot the original data and the fitted curve
#plot(System_sizes, fluct_re, xlabel=L"L", label=L"\mathcal{F}_ϕ(L)",title="Fluctuations of ϕ (\$S^z\$) with U=$arg1",dpi=400)
#plot!(System_sizes, fit_F, label="Log Fit: $slope_round*log(L) + $intercept_round ", linestyle=:dash)
#savefig("plot_fluc_phi_spin.pdf")


# Open a file in write mode
open("results_JR/flc_phi_nint_delta_$arg1.txt", "w") do file
    for i in 1:length(System_sizes)
        println(file, "$(System_sizes[i]),$(fluct_re[i])")
   end
    #println(slope, intercept)
end


#=
# Initialize empty arrays to store the read data
read_x = []
read_y = []
# Open the file in read mode
open("flc_phi_U_05.txt", "r") do file
    println(eachline(file))
    for line in eachline(file)
        values = split(line, ",")
        push!(read_x, parse(Float64, values[1]))
        push!(read_y, parse(Float64, values[2]))
    end
end
# Print the read data to verify
#println(read_x)
#println(read_y)
# Plot y as a function of x
#plot(read_x, read_y, title="Plot of y as a function of x", xlabel="x", ylabel="y", label="y(x)")
  =#
