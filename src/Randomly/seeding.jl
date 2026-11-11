"""
    SeededSample(fields, seeds; rng_type=Xoshiro)

Attach per-field seeds and the RNG type to sampled model inputs.

Return this record from custom samplers opting into seeded [`map_mc`](@ref) calls.
The record carries metadata; constructing it does not create or advance an RNG.

- `fields`: tuple of positional model arguments; replay with `model(draw.fields...)`.
- `seeds`: named tuple of the integer seeds actually used, one per input field.
- `rng_type`: RNG type used for these fields, defaulting to `Random.Xoshiro`.

Keep the seed schema and RNG type fixed within an ensemble. Built-in seeded samplers
use fresh Xoshiro generators, but custom samplers may report another `AbstractRNG` type.
See also [`make_sample_seeds`](@ref).
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

Pass a tuple of field-name symbols and a nonnegative master `seed` representable as
`UInt64`. Return a vector of named tuples containing one `UInt64` seed per field.
`n` may be zero; [`map_mc`](@ref) requires a positive sample count.

Each field has its own Xoshiro stream, seeded by the master seed and field-name bytes,
not Julia's `hash`. Extending `n`, adding/removing other fields or reordering names
preserves existing columns; renaming a field changes its stream. The helper creates
only seed records, not GRFs, and does not consume the default RNG.

Store the actual table and RNG/environment provenance for replay. Equal seeds do not
imply equal fields across different grids or covariance settings, and cross-version
random-stream stability is not promised. See also [`SeededSample`](@ref).

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
    n >= 0 || throw(ArgumentError("Sample count must be nonnegative"))
    generators = NamedTuple{fields}(map(fields) do field
        Xoshiro([UInt64(seed); UInt64.(codeunits(String(field)))])
    end)
    return [map(rng -> rand(rng, UInt64), generators) for _ in 1:n]
end
