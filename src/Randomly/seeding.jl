"""
    SeededSample(fields, seeds; rng_type=Xoshiro)

Attach per-field seeds and the RNG type to sampled model inputs.

Return this from a sampler called as `sample(; seeds)` by seeded [`map_mc`](@ref).

- `fields`: tuple of model arguments; replay with `model(draw.fields...)`.
- `seeds`: named tuple of the seeds used, one per field.
- `rng_type`: RNG type seeded with them.
"""
struct SeededSample{F,S<:NamedTuple,R<:AbstractRNG}
    fields::F
    seeds::S
    rng_type::Type{R}
end

SeededSample(fields, seeds; rng_type=Xoshiro) = SeededSample(fields, seeds, rng_type)

"""
    make_sample_seeds(n, fields; seed)

Create `n` per-draw seed records for the named input `fields`.

Return a vector of named tuples with one `UInt64` seed per field. Each field's seeds
come from its own stream, seeded by the master `seed` and the field name. Thus changing
`n` or adding, removing or reordering other fields leaves a field's seeds unchanged.

# Examples

```jldoctest
julia> using WhereTheWaterFlows

julia> const WWF = WhereTheWaterFlows; const WWFR = WWF.Randomly;

julia> a = WWFR.make_sample_seeds(3, (:bed, :pressure); seed=42);

julia> b = WWFR.make_sample_seeds(5, (:melt, :pressure, :bed); seed=42);

julia> all(a[i].bed == b[i].bed && a[i].pressure == b[i].pressure for i in 1:3)
true
```
"""
function make_sample_seeds(n::Integer, fields::Tuple{Vararg{Symbol}}; seed::Integer)
    # one RNG per field, seeded from master seed and the name's bytes (not `hash`)
    rngs = NamedTuple{fields}(map(f -> Xoshiro([UInt64(seed); codeunits(String(f))]), fields))
    return [map(rng -> rand(rng, UInt64), rngs) for _ in 1:n]
end
