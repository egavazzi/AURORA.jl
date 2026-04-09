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
- ``D(z, E, \mu)`` is the spatial diffusion coefficient (due to finite pitch-angle widths)
- ``Q(z, \mu, t, E)`` is the source term from energy cascading (ionization secondaries and
  degraded primaries from higher energies)

The spatial domain is discretized on a non-uniform altitude grid with ``I_e`` defined at
cell centers. Different finite difference schemes are used depending on the nature of each
term: the streaming term uses directional (upwind) differences, while the diffusion term
uses central differences.

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

## Spatial discretization

The altitude grid in AURORA is non-uniform, with finer spacing at lower altitudes where
collision rates are highest.

The **streaming term** ``-\mu v \, \partial I_e / \partial z`` uses directional
(upwind) finite differences: a backward difference for upward-going beams (``\mu > 0``)
and a forward difference for downward-going beams (``\mu < 0``). This choice is dictated
by the physics — the derivative should be evaluated using information from the direction
the beam is coming from.

The **diffusion term** ``\partial_z(D \, \partial_z I_e)`` uses central differences,
as diffusion spreads information symmetrically in both directions.

## Crank-Nicolson discretization

AURORA uses the Crank-Nicolson scheme — a second-order, unconditionally stable
semi-implicit time integrator. It descretizes the time derivative as:

```math
\frac{1}{v}\frac{\vec{I_e}^{n+1} - \vec{I_e}^{n}}{\Delta t}
= \frac{1}{2}\left(-\mathcal{L}\,\vec{I_e}^{n+1} + \vec{Q}^{n+1}\right)
+ \frac{1}{2}\left(-\mathcal{L}\,\vec{I_e}^{n} + \vec{Q}^{n}\right)
```

i.e., the right-hand side is the average at the current and next time steps.
Rearranging leads to the linear system solved at each time step:

```math
M_{\text{lhs}} \, \vec{I_e}^{n+1} = M_{\text{rhs}} \, \vec{I_e}^{n} + \frac{\vec{Q}^{n+1} + \vec{Q}^{n}}{2}
```

where:
- ``M_{\text{lhs}} = D_{\Delta t} + \mathcal{L}/2``
- ``M_{\text{rhs}} = D_{\Delta t} - \mathcal{L}/2``
- ``D_{\Delta t} = \mathrm{diag}(1 / (v\,\Delta t))`` is the time-step diagonal (with units 1/m)

At each energy level, these matrices are assembled once and the system is solved using a
sparse LU factorization (KLU from SuiteSparse).

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

`update_D!(D, model)` computes the spatial diffusion coefficient for each energy and
pitch-angle bin. Within a discrete bin, electrons span a range of parallel velocities
``v\cos\theta``. This spread causes electrons in the same beam to arrive at different
altitudes at different times. The diffusion coefficient D models this spread of travel times:

```math
D[iE, i\mu] = \frac{(\Delta t_{\text{arrival}} / 4)^2}{\bar{t}_{\text{arrival}}}
```

where ``\Delta t_{\text{arrival}}`` is the range of arrival times across the bin and
``\bar{t}_{\text{arrival}}`` is the mean. ``D`` vanishes for infinitely narrow bins.

### Source term Q

`update_Q!(Q, iE, ...)` accumulates contributions from higher energies into the source
vector for the current energy ``iE``:

1. **Electron-electron losses**: Energy transferred from bin ``iE+1`` appears as a source
   in bin ``iE``.
2. **Inelastic scattering**: Electrons losing discrete amounts of energy to excitation of
   neutral states.
3. **Ionization**:
   - *Secondary electrons*: produced when a neutral is ionized and emitted isotropically
     over all pitch angles. The energy distribution of secondaries is species-dependent.
     Pre-computed *cascading transfer matrices* (one per species) map the ionization rate
     at energy ``E'`` to the secondary electron production at energy ``E < E'``.
   - *Degraded primaries*: the ionizing electron continues in the same beam at reduced
     energy (ionization potential + energy given to the secondary).

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

## Boundary conditions

- **Top boundary** (``z = z_{\max}``):
  - *Downward-going beams* (``\mu < 0``): Dirichlet condition — the flux is prescribed by
    the `Ie_top` array from the [`InputFlux`](@ref).
  - *Upward-going beams* (``\mu > 0``): Zero-gradient condition
    (``\partial I_e / \partial z = 0``), allowing electrons to exit freely at the top of
    the simulation domain.

- **Bottom boundary** (``z = z_{\min}``): Dirichlet condition — the flux is held at its
  initial value.

## Implementation notes

- **Sparse factorizations**: The block-structured transport matrices are solved with KLU
  factorizations from SuiteSparse.
- **Loop partitioning**: Long time-dependent simulations may be split into `n_loop` chunks
  so that the working arrays fit within the requested memory budget.
- **Cached data**: Scattering and cascading data are saved to disk and reused when the 
  relevant grids match, avoiding repeated setup work.