const STATIC_STRUCTURE_FACTOR_DEFINITION = (
    fields=("structureFactorBySpeciesPair", "structureFactorByMomentum"),
    optional_fields=(), complex_fields=("structureFactorBySpeciesPair",),
    metadata=("momentum_indices", "momenta"), reduction="mean",
)

"""Validate saved component-pair and total histories against their Fourier grid."""
function validate_static_structure_factor(histories, metadata)
    pairs = histories["structureFactorBySpeciesPair"]
    total = histories["structureFactorByMomentum"]
    pairs isa AbstractArray{<:Complex} || throw(ArgumentError(
        "StaticStructureFactor/structureFactorBySpeciesPair must be a complex history"))
    ndims(pairs) == 4 && ndims(total) == 2 || throw(DimensionMismatch(
        "StaticStructureFactor requires component × component × momentum × block and momentum × block histories"))
    size(pairs, 1) == size(pairs, 2) && size(pairs, 1) > 0 && size(pairs, 3) > 0 ||
        throw(DimensionMismatch("StaticStructureFactor requires nonempty, matching component axes and a nonempty momentum axis"))
    size(total) == (size(pairs, 3), size(pairs, 4)) || throw(DimensionMismatch(
        "StaticStructureFactor pair and total histories must have matching momentum and block axes"))
    indices, momenta = metadata["momentum_indices"], metadata["momenta"]
    indices isa AbstractMatrix{<:Integer} && !(eltype(indices) <: Bool) ||
        throw(ArgumentError("StaticStructureFactor momentum_indices must be an integer matrix"))
    momenta isa AbstractMatrix{<:Real} && all(isfinite, momenta) ||
        throw(ArgumentError("StaticStructureFactor momenta must be a finite real matrix"))
    size(indices) == size(momenta) && size(indices, 1) == size(pairs, 3) && size(indices, 2) > 0 ||
        throw(DimensionMismatch("StaticStructureFactor momentum metadata must have matching momentum × dimension shapes"))
    return nothing
end

"""
    static_structure_factor_summary(histories, counts, metadata; min_blocks,
        blocking_diagnostics=false)

Summarize validated raw, site-normalized block means with sample-count weights.
Analyze each pair entry's real and imaginary parts separately and the saved total
independently using maximum repeated-blocking standard errors. Covariances are
not retained; no normalization, subtraction, or symmetrization is applied.
"""
function static_structure_factor_summary(histories, counts, metadata;
                                         min_blocks, blocking_diagnostics::Bool=false)
    pairs = histories["structureFactorBySpeciesPair"]
    summary = mean_observable_summary(Dict(
        "structureFactorBySpeciesPairReal" => real.(pairs),
        "structureFactorBySpeciesPairImaginary" => imag.(pairs),
        "structureFactorByMomentum" => histories["structureFactorByMomentum"],
    ), counts, metadata; min_blocks, blocking_diagnostics)
    summary["structureFactorBySpeciesPair"] = complex.(
        pop!(summary, "structureFactorBySpeciesPairReal"),
        pop!(summary, "structureFactorBySpeciesPairImaginary"))
    return summary
end
