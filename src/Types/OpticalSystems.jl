"""
    OpticalSystem(elements)

An ordered sequence of optical elements representing a complete
optical system. When called on a `ScalarField`, applies each
element in sequence.

# Examples
```jldoctest
julia> system = OpticalSystem([
           FreeSpace(0.05),
           ThinLens(0.1),
           FreeSpace(0.1),
       ]);

julia> length(system.elements)
3
```
"""
struct OpticalSystem
    elements :: Vector{AbstractOpticalElement}
end

function (system::OpticalSystem)(field::ScalarField) :: ScalarField
  for element in system.elements
    field = apply(element, field)
  end
  return field
end

apply(e::AbstractOpticalElement, f::ScalarField) =
  error("$(typeof(e)) must implement apply()")
