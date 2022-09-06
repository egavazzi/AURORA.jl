function v_of_E(E)
    mₑ = 9.10939e-31;
    qₑ = 1.6021773e-19;

    v = (2 * qₑ * abs(E) / mₑ) .^ (1/2) .* sign(E);
    return v
end