# Input Flux

Types and functions for specifying the precipitating electron flux. See the
[Custom input flux](@ref "Custom Input Flux") tutorial for usage examples.

## Spectrum types

```@docs; canonical=false
AbstractSpectrum
FlatSpectrum
GaussianSpectrum
MaxwellianSpectrum
FileSpectrum
evaluate_spectrum
```

## Modulation types

```@docs; canonical=false
AbstractModulation
ConstantModulation
SinusoidalFlickering
SquareFlickering
SmoothOnset
apply_modulation
```

## InputFlux

```@docs; canonical=false
InputFlux
compute_flux
```

## File loading

```@docs; canonical=false
Ie_top_from_file
```
