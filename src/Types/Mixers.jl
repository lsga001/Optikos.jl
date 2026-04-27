"""
    SFGCrystal()

Sum Frequency Generation crystal. Models the nonlinear interaction where
two input fields E₁ and E₂ are mixed to produce an output field proportional
to their product:

```math
E_{out}(\\mathbf{r}) = E_1(\\mathbf{r}) \\cdot E_2(\\mathbf{r})
```

In the context of coherence measurements, placing E₂(r) = E₁*(-r) gives
access to the cross-spectral density W(r, -r).

# Notes
This is a simplified model of SFG that assumes:
- Perfect phase matching
- Undepleted pump approximation
- Thin crystal (no propagation effects within the crystal)

# References
Boyd, *Nonlinear Optics*, 3rd ed., §2.4
"""
struct SFGCrystal <: AbstractMixer end
