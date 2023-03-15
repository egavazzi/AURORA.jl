function v_of_E(E)
	mₑ = 9.10939e-31;
	qₑ = 1.6021773e-19;

	v = (2 * qₑ * abs(E) / mₑ) .^ (1/2) .* sign(E);
	return v
end

function CFL_criteria(t, h_atm, v, CFL_number=64)
    dt = t[2] - t[1]
    dz = h_atm[2] - h_atm[1]

    # The Courant-Freidrichs-Lewy (CFL) number normally hase to be small (<4) to ensure numerical
    # stability. However, as a Crank-Nicolson scheme is always stable, we can take a bigger CFL. We
    # should be careful about numerical accuracy though.
    # For Gaussian inputs (or similar), it seems that the CFL can be set to 64 without major effects
    # on the results, while reducing computational time tremendously
    CFL = v * dt / dz
    n_factors = 2 .^ collect(0:22)
    iFactor = 1
    # This while loop effectively reduces dt by a factor of 2 at each iteration and check if the new
    # CFL is < 64. If not, it continues reducing dt.
    t_finer = t
    while (CFL > CFL_number) && (iFactor < length(n_factors))
        t_finer = range(t[1], t[end], length(t) * n_factors[iFactor] + 1 - n_factors[iFactor])
        dt = t_finer[2] - t_finer[1]
        CFL = v * dt / dz
        iFactor += 1
    end
    CFL_factor = n_factors[max(1, iFactor - 1)]

    return t_finer, CFL_factor
end


using HCubature
"""
    mu_avg(θ_lims)

Calculate the cosinus of the center of each pitch-angle beams delimited by θ_lims. This is
for isotropically distributed fluxes within each beam, i.e the fluxes are weighted by sin(θ)

# Calling
`μ_center = mu_avg(θ_lims) `

# Inputs
- `θ_lims` : pitch-angle limits of all the beams, range or vector [n_beams + 1]

# Outputs
- `μ_center` : cosine of the center of all the pitch-angle beams, vector [n_beams + 1]
"""
function mu_avg(θ_lims)
    μ_center = zeros(length(θ_lims) - 1)
    for i in eachindex(μ_center)
        μ_center[i] = hcubature(x -> cosd.(x) .* sind.(x), [θ_lims[i]], [θ_lims[i + 1]])[1][1] ./
                      hcubature(x -> sind.(x), [θ_lims[i]], [θ_lims[i + 1]])[1][1]
    end
    return μ_center
end


## ====================================================================================== ##


function save_parameters(altitude_max, θ_lims, E_max, B_angle_to_zenith, t_sampling, t, n_loop, INPUT_OPTIONS, savedir)
	savefile = string(savedir, "/", "parameters.txt")
    open(savefile, "w") do f
        write(f, "altitude_max = $altitude_max \n")
        write(f, "θ_lims = $θ_lims \n")
        write(f, "E_max = $E_max \n")
        write(f, "B_angle_to_zenith = $B_angle_to_zenith \n")
        write(f, "\n")
        write(f, "t_sampling = $t_sampling \n")
        write(f, "t = $t \n")
        write(f, "n_loop = $n_loop \n")
        write(f, "\n")
        write(f, "input_options = $INPUT_OPTIONS")
    end
end


using MAT
function save_neutrals(h_atm, n_neutrals, ne, Te, savedir)
    savefile = string(savedir, "/", "neutral_atm.mat")
    file = matopen(savefile, "w")
        write(file, "h_atm", h_atm)
        write(file, "nN2", n_neutrals.nN2)
        write(file, "nO2", n_neutrals.nO2)
        write(file, "nO", n_neutrals.nO)
        write(file, "ne", ne)
        write(file, "Te", Te)
    close(file)
end


using MAT
using Printf
function save_results(Ie, E, t, μ_lims, h_atm, I0, μ_scatterings, i, CFL_factor, savedir)
    # Extract the time array for the current loop
	t_run = collect(t .+ t[end] * (i - 1))

    # Reduce t_run and Ie to match the t_sampling
    t_run = t_run[1:CFL_factor:end]
    Ie_save = Ie[:, 1:CFL_factor:end, :]

    savefile = joinpath(savedir, (@sprintf "IeFlickering-%02d.mat" i))
	file = matopen(savefile, "w")
		write(file, "Ie_ztE", Ie_save)
		write(file, "E", E)
		write(file, "t_run", t_run)
		write(file, "mu_lims", μ_lims)
		write(file, "h_atm", h_atm)
		write(file, "I0", I0)
		write(file, "mu_scatterings", μ_scatterings)
	close(file)
end


## ====================================================================================== ##


using Interpolations

function interp1(X, V, Xq)
    knots = (X,)
    itp = interpolate(knots, V, Gridded(Linear()))
    itp[Xq]
end

function interp2(X, Y, V, Xq, Yq)
    knots = (X,Y)
    itp = interpolate(knots, V, Gridded(Linear()))
    itp[Xq, Yq]
end


## ====================================================================================== ##


function square(x)
    ifelse(mod2pi(x) < π, 1.0, -1.0)
end


## ====================================================================================== ##

# function from Björn for smooth input onset
function f_smooth_transition(x, a = 0, b = 1)
    if (b - a) == 0
        y = 1
        return y
    end

    x = (x - a) / (b - a)
    psi0P = psi(x)
    psi1N = psi(1 - x)

    y = psi0P / (psi0P + psi1N)

    return y
end

function psi(x)
    PSI = exp(-1 / x)
    if x <= 0
        PSI = 0
    end
    return PSI
end
f_smooth_transition(0, 1, 0)