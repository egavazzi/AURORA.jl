function v_of_E(E)
  mₑ = 9.10939e-31;
  qₑ = 1.6021773e-19;

  v = (2 * qₑ * abs(E) / mₑ) .^ (1/2) .* sign(E);
  return v
end


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