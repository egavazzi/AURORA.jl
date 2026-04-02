# [Internals](@id Internals)

This page describes how AURORA implements the transport model in code: matrix assembly,
time stepping, the energy loop, and the main data flow through the solver.

The physical meaning of the transport equation and the degradation processes is summarized
in the [Physics](@ref "Transport Equation") pages. Here the focus is on the numerical
realization.

## The transport equation

AURORA solves the multi-stream electron transport equation for the differential electron
flux ``I_e(z, \mu, t, E)`` as a function of altitude ``z``, pitch-angle cosine ``\mu``,
time ``t`` and energy ``E``. In its semi-discretized form (continuous in ``z`` and ``t``,
discrete in ``\mu`` and ``E``):

```math
\frac{\partial I_e}{\partial t} =
    \underbrace{-\mu \, v \, \frac{\partial I_e}{\partial z}}_{\text{streaming}}
    + \underbrace{\frac{\partial}{\partial z}\left(D \frac{\partial I_e}{\partial z}\right)}_{\text{diffusion}}
    - \underbrace{A \cdot I_e}_{\text{losses}}
    + \underbrace{B \cdot I_e}_{\text{pitch-angle scattering}}
    + \underbrace{Q}_{\text{sources}}
```

where:
- ``v = v(E)`` is the electron speed at energy ``E``
- ``\mu`` is the cosine of the pitch angle
- ``A(z, E)`` contains the total loss frequency from collisions with neutrals and thermal
  electrons
- ``B(z, E, \mu \to \mu')`` is the pitch-angle scattering matrix
- ``D(z, E, \mu)`` is the spatial diffusion coefficient
- ``Q(z, \mu, t, E)`` is the source term from energy cascading (ionization secondaries and
  degraded primaries from higher energies)

The equation is solved on a staggered grid: the flux ``I_e`` is defined at cell centers in
altitude, while the diffusion operator uses a standard second-order finite difference
stencil.

## Matrix representation

The state vector at a given energy ``E`` is the flux flattened over altitude and pitch
angle: ``\vec{I_e} \in \mathbb{R}^{n_z \times n_\mu}``. The spatial operator is assembled
as a block matrix of size ``(n_z n_\mu) \times (n_z n_\mu)``:

- **Diagonal blocks** (one per beam ``\mu_i``): tridiagonal in the altitude dimension,
  containing the streaming, diffusion, and loss terms.
- **Off-diagonal blocks**: sparse, connecting different pitch-angle beams via the
  scattering matrix ``B``.

This gives the system:

```math
\frac{d\vec{I_e}}{dt} = -\mathcal{L} \, \vec{I_e} + \vec{Q}
```

where ``\mathcal{L}`` is the combined spatial/loss/scattering operator.

## Crank-Nicolson discretization

AURORA uses the Crank-Nicolson scheme — a second-order, unconditionally stable
semi-implicit time integrator:

```math
\left(I + \frac{\Delta t}{2} \mathcal{L}\right) \vec{I_e}^{n+1}
= \left(I - \frac{\Delta t}{2} \mathcal{L}\right) \vec{I_e}^{n}
+ \Delta t \, \vec{Q}
```

This is rewritten as:

```math
M_{\text{lhs}} \, \vec{I_e}^{n+1} = M_{\text{rhs}} \, \vec{I_e}^{n} + \Delta t \, \vec{Q}
```

At each energy level, the matrices ``M_{\text{lhs}}`` and ``M_{\text{rhs}}`` are assembled
once and the system is solved using a sparse LU factorization (KLU from SuiteSparse).

For the **steady-state** case, the time derivative vanishes and we solve:

```math
\mathcal{L} \, \vec{I_e} = \vec{Q}
```

directly via the same KLU factorization.

## Matrix construction

### Loss matrix A

`update_A!(A, iE, cross_sections, ionosphere)` computes the total loss frequency at each
altitude for energy bin `iE`:

```math
A(z) = \sum_{\text{species}} n_s(z) \sum_j \sigma_j^s(E) \cdot v(E) + \nu_{ee}(z, E)
```

where the sum over ``j`` includes elastic, inelastic, and ionization cross sections for
each neutral species (N₂, O₂, O), and ``\nu_{ee}`` is the electron-electron collision
frequency from the Swartz & Nisbet (1971) formula.

### Scattering matrix B

`update_B!(B, iE, ...)` constructs the pitch-angle redistribution matrix. For each
collision type, the differential cross section (from Itikawa data) is converted from
scattering-angle space to beam-to-beam transfer probabilities:

```math
B(\mu_i \to \mu_j) = \sum_s n_s(z) \cdot v(E) \sum_k \sigma_k^s(E) \cdot P_k(\mu_i \to \mu_j)
```

where ``P_k`` is the scattering probability for collision process ``k``.

### Diffusion coefficient D

The diffusion coefficient accounts for pitch-angle mixing due to small-angle scattering.
It is computed from the elastic scattering cross sections.

### Source term Q

`update_Q!(Q, iE, ...)` accumulates contributions from higher energies into the source
vector for the current energy ``iE``:

1. **Electron-electron losses**: Energy transferred from bin ``iE+1`` appears as a source
   in bin ``iE``.
2. **Inelastic scattering**: Electrons losing discrete amounts of energy to excitation of
   neutral states.
3. **Ionization fragments**:
   - *Secondary electrons*: produced when a neutral is ionized, with a Lorentzian-like
     energy distribution (width ~11.4 eV per Itikawa 1986). The pre-computed *cascading
     transfer matrices* (one per species) map the ionization rate at energy ``E'`` to the
     secondary electron production at energy ``E < E'``.
   - *Degraded primaries*: the ionizing electron loses its ionization energy and continues
     at lower energy.

These cascading transfer matrices are computed once (by `cascading_N2()`, `cascading_O2()`,
`cascading_O()`) and cached in `internal_data/e_cascading/` for reuse.

## The energy-descending loop

The energy loop in `energy_loop!` iterates from the highest energy bin to the lowest. This
ordering is essential because of the cascading mechanism: the source term ``Q(E)`` at
energy ``E`` includes contributions from *all* higher energies. By solving from high to
low, these contributions are fully computed before they are needed.

At each energy step:

1. **Update matrices** — recompute A, B, D for the current energy (cross sections are
   energy-dependent).
2. **Solve** — advance the Crank-Nicolson scheme (time-dependent) or invert the system
   (steady-state) to obtain ``I_e`` at this energy.
3. **Update Q** — compute the cascading/degradation contributions from this energy to all
   lower energies.

## Implementation notes

- **Sparse factorizations**: The block-structured transport matrices are solved with KLU
  factorizations from SuiteSparse.
- **Loop partitioning**: Long time-dependent simulations may be split into `n_loop` chunks
  so that the working arrays fit within the requested memory budget.
- **Cached data**: Scattering and cascading data are saved to disk and reused when the 
  relevant grids match, avoiding repeated setup work.

## SimulationCache data flow

The `SimulationCache` is the central workspace allocated by `initialize!`. Here is how
data flows through it during a time-dependent simulation:

```
initialize!(sim)
│
├── Compute Ie_top[n_μ, n_t, n_E]     # Full input flux boundary condition
├── Allocate Ie[n_z·n_μ, n_t_chunk, n_E]   # Solution array for one time chunk
├── Allocate I0[n_z·n_μ, n_E]         # Initial condition (carried between loops)
├── Build TransportMatrices (A, B, D, Q, Ddiffusion)
└── Build SolverCache (KLU factorizations)

solve! loop:
│
├── Loop i_loop = 1:n_loop
│   ├── Extract Ie_top_chunk from Ie_top for this time window
│   ├── Reset Ie to zero (or carry I0 from previous loop)
│   │
│   ├── energy_loop!(descending in E)
│   │   └── For each energy:
│   │       ├── Update A, B, D, Ddiffusion ← cross_sections[iE]
│   │       ├── Assemble Mlhs, Mrhs using TransportMatrices
│   │       ├── Solve → writes into Ie[:, :, iE]
│   │       └── Update Q[:, :, iE-1:1] ← cascading from Ie[:, :, iE]
│   │
│   ├── Downsample Ie → Ie_save (keep every CFL_factor-th time step)
│   ├── Save Ie_save to disk
│   └── Store I0 ← Ie[:, end, :] for next loop
│
└── Done
```

## Boundary conditions

- **Top boundary** (``z = z_{\max}``): The incoming (downward) flux is prescribed by the
  `Ie_top` array from the [`InputFlux`](@ref). Upward-going beams at the top boundary are
  free (outflow).

- **Bottom boundary** (``z = z_{\min}``): Upward-going flux is determined by reflection
  of the downward flux (backscatter from the dense atmosphere below).
