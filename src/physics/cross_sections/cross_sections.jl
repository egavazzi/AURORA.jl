"""
    zero_below_energy_loss!(σ, E_centers, E_loss, label)

Set `σ` to zero in the energy bins whose center is below `E_loss`, where the collision
cannot happen, and warn (once per session for each `label`) when any of those values was
non-zero: it points to a cross-section fit that starts below its energy loss.
"""
function zero_below_energy_loss!(σ, E_centers, E_loss, label)
    below = (E_centers .< E_loss) .& (σ .> 0)
    if any(below)
        @warn "$(label) has a non-zero cross section below its energy loss of $(E_loss) eV, \
               at bin centers $(round.(E_centers[below]; digits = 3)) eV. These values are set \
               to zero." maxlog = 1 _id = (label, E_loss)
        σ[below] .= 0
    end
    return σ
end
