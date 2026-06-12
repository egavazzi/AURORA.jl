import JLD2

"""
    ExprLaw

Serialization-friendly wrapper around a user-supplied physical law. This can be a density profile,
a cascading secondary-distribution law, a phase-function generator, etc. It stores the law's
source as a string so the model can round-trip through JLD2: only the string is written to disk,
and the callable is rebuilt by parsing and evaluating that string on load.

Build one with the [`@law`](@ref) macro:

```julia
density         = @law z -> 1e18 .* exp.(-(z .- 100e3) ./ 30e3)
secondary_e_law = @law (E_s, E_p) -> 1 / (11.4^2 + E_s^2)
```

To ensure reproducibility, laws must be self-contained: closed-form expressions referencing
only their arguments and names available in `AURORA`/`Base`. A law that needs to carry e.g.
parameters should be constructed as a functor struct instead. To ensure this, trying to build
a `@law` that captures local variables will throw an error.

!!! warning "Loading evaluates code"
    Reloading a saved model evaluates the stored law source, so a `physics_state.jld2` file
    from an untrusted source can execute arbitrary code. Only load files you trust.
"""
struct ExprLaw
    src::String
    f
    function ExprLaw(src::AbstractString, f)
        if f isa Function && !isempty(fieldnames(typeof(f)))
            captured = join(fieldnames(typeof(f)), ", ")
            throw(ArgumentError(
                    "@law captures local variable(s) ($captured), which won't be possible to reconstruct " *
                    "from the saved source. Hold parameters in a functor struct instead. See the documentation " *
                    "about modifying species and grids."))
        end
        return new(String(src), f)
    end
end

ExprLaw(src::AbstractString) = ExprLaw(src, build_law(src))

# Rebuild the callable from its source. Evaluated in the AURORA module so that package
# helpers (e.g. phase_fcn_N2) resolve when a saved model is reloaded in a fresh session.
build_law(src::AbstractString) = Core.eval(@__MODULE__, Meta.parse(src))

# invokelatest keeps calls valid after a reload, where `f` is eval'd in a newer world age
# than the calling code. Laws run only at initialize! (densities/phase once, the cascading
# law inside a disk-cached integration), so the dynamic-dispatch cost is negligible.
(l::ExprLaw)(args...) = Base.invokelatest(l.f, args...)

# Render a captured law expression as clean source: strip line-number metadata and unwrap the
# single-statement block that `->` wraps its body in, so `@law z -> exp(-z)` is stored as
# "z->exp(-z)" rather than "z->begin … end".
law_source(expr) = expr isa Expr ? string(collapse_blocks(Base.remove_linenums!(copy(expr)))) :
                                   string(expr)

function collapse_blocks(x)
    x isa Expr || return x
    x.head === :block && length(x.args) == 1 && return collapse_blocks(x.args[1])
    return Expr(x.head, Any[collapse_blocks(a) for a in x.args]...)
end

"""
    @law expr

Wrap a law `expr` (e.g. `z -> exp(-z)`) into an [`ExprLaw`](@ref), capturing its source text
so it can be saved and reloaded verbatim. Use this instead of a bare anonymous function for
any density profile, cascading law, or phase-function generator.
"""
macro law(expr)
    return :( ExprLaw($(law_source(expr)), $(esc(expr))) )
end

# JLD2 persists only the source string; `f` is rebuilt by ExprLaw(src) on load.
struct ExprLawSerialization
    src::String
end
JLD2.writeas(::Type{ExprLaw}) = ExprLawSerialization
Base.convert(::Type{ExprLawSerialization}, l::ExprLaw) = ExprLawSerialization(l.src)
Base.convert(::Type{ExprLaw}, s::ExprLawSerialization) = ExprLaw(s.src)

Base.show(io::IO, l::ExprLaw) = print(io, "@law ", l.src)

# A law is reproducible unless it is a bare anonymous function/closure. Functors (callable
# structs, including ExprLaw), named functions, MSISDensity and VectorDensity all pass.
is_anonymous(f) = f isa Function && startswith(string(nameof(f)), "#")

function require_reproducible(law, role::AbstractString)
    is_anonymous(law) && throw(ArgumentError(
        "$role is a bare anonymous function and cannot be saved for reproducibility. \
         Wrap it with @law (e.g. `@law $role`), or pass a functor or named function."))
    return law
end
