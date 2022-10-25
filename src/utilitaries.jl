function v_of_E(E)
	mₑ = 9.10939e-31;
	qₑ = 1.6021773e-19;

	v = (2 * qₑ * abs(E) / mₑ) .^ (1/2) .* sign(E);
	return v
end


## ====================================================================================== ##


function save_parameters(altitude_max, θ_lims, E_max, B_angle_to_zenith, t, n_loop, input_file, savedir)
	savefile = string(savedir, "/", "parameters.txt")
    open(savefile, "w") do f
        write(f, "altitude_max = $altitude_max \n")
        write(f, "θ_lims = $θ_lims \n")
        write(f, "E_max = $E_max \n")
        write(f, "B_angle_to_zenith = $B_angle_to_zenith \n")
        write(f, "\n")
        write(f, "t = $t \n")
        write(f, "n_loop = $n_loop \n")
        write(f, "\n")
        write(f, "input_file = $input_file")
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
function save_results(Ie, E, t, μ_lims, h_atm, I0, μ_scatterings, n_loop, savedir, i)
	t_run = collect(t .+ t[end] * (n_loop - 1))

	savefile = string(savedir, "/", (@sprintf "IeFlickering-%02d.mat" i))
	file = matopen(savefile, "w")
		write(file, "Ie_ztE", Ie)
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