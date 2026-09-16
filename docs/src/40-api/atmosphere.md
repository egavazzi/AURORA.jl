# Ionospheric Background Models

## MSIS (Neutral Atmosphere)
```@docs; canonical=false
find_msis_file
NeutralAtmosphere
run_msis
read_msis_file
read_ccmc_msis
```

## IRI (Ionosphere)
```@docs; canonical=false
find_iri_file
ElectronProfile
run_iri
read_iri_file
read_ccmc_iri
```

## Species density profiles and laws
```@docs; canonical=false
DensityProfile
@law
ExprLaw
```

## Neutral species and their collision channels
```@docs; canonical=false
NeutralSpecies
N2Species
O2Species
OSpecies
CollisionChannel
AURORA.channel
channel_names
ionizing_channels
default_channels
default_elastic_cross_section
default_secondary_law
default_cascading_spec
```
