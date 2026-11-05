"""
    describe_density(summary; io=stdout) -> nothing

Print an observable inventory, all Density estimates, production counts, and
sample-count-weighted mean/repeated-blocking conventions in aligned sections.
Pass a Density dictionary such as `table.Density[1]`,
or `missing` to report that no measurements are available. The input is not modified.

Species IDs label the rows when available; otherwise use component numbers and
state that species IDs are unavailable. Nonfinite uncertainties are displayed as
"uncertainty unavailable". Values must be Measurements, with matching species
vectors for density and particle number. No estimates are recomputed.
"""
function describe_density(summary::Union{Missing,AbstractDict}; io::IO=stdout)
    if ismissing(summary)
        println(io, "No Density measurements available.")
        return nothing
    end
    density, particles = summary["densityBySpecies"], summary["particleNumberBySpecies"]
    density isa AbstractVector && particles isa AbstractVector || throw(ArgumentError(
        "Density and particle number must be species vectors"))
    length(density) == length(particles) || throw(DimensionMismatch(
        "Density and particle number must have matching species counts"))
    species = get(summary, "species", nothing)
    has_species = species !== nothing && !ismissing(species)
    if has_species
        species isa AbstractVector{<:Integer} && !(eltype(species) <: Bool) &&
            all(>(0), species) && allunique(species) || throw(ArgumentError(
            "species must contain positive unique integer IDs"))
        length(species) == length(density) || throw(DimensionMismatch(
            "species IDs must match the Density species count"))
    end
    headers = ["Species", "Density (particles/site)", "Particle number (particles)"]
    rows = [[has_species ? string(species[i]) : "Component $i",
             format_observable_estimate(density[i]), format_observable_estimate(particles[i])]
            for i in eachindex(density)]
    total = format_observable_estimate(summary["totalParticleNumber"])

    inventory = observable_table_lines(["Quantity", "Shape", "Axes", "Units", "Description"], [
        ["densityBySpecies", string(length(density)), "species", "particles/site", "Species density"],
        ["particleNumberBySpecies", string(length(particles)), "species", "particles", "Species particle number"],
        ["totalParticleNumber", "scalar", "—", "particles", "Total particle number"],
    ])
    measurements = observable_table_lines(headers, rows)
    has_species || push!(measurements, "Species IDs unavailable; rows are labeled by component.")
    append!(measurements, ["", "Total particle number: $total particles"])
    statistics = observable_production_lines(summary, "Density")
    conventions = mean_observable_conventions(summary,
        "Density is particles/site; particle numbers are in particles.")
    print_observable_report(io, "Density", (
        ("OBSERVABLES", inventory), ("MEASUREMENTS", measurements),
        ("PRODUCTION STATISTICS", statistics), ("ESTIMATOR CONVENTIONS", conventions),
    ); width=maximum(textwidth, vcat(inventory, measurements, statistics)))
    return nothing
end

"""
    describe_energy(summary; io=stdout) -> nothing

Print an observable inventory and all saved diagonal, kinetic, and total Energy
estimates, both extensive and per site, in model energy units. Include production
counts and sample-count-weighted mean/repeated-blocking conventions in aligned
sections. Accept `missing` for unavailable measurements; show nonfinite errors
as "uncertainty unavailable". Estimates must
be scalar Measurements. No values are recomputed and the input is not modified.
"""
function describe_energy(summary::Union{Missing,AbstractDict}; io::IO=stdout)
    if ismissing(summary)
        println(io, "No Energy measurements available.")
        return nothing
    end
    components = (("Diagonal", "diagonalEnergy"), ("Kinetic", "kineticEnergy"), ("Total", "totalEnergy"))
    rows = [[label, format_observable_estimate(summary[field]),
             format_observable_estimate(summary[field * "PerSite"])]
            for (label, field) in components]
    inventory = observable_table_lines(["Quantity", "Shape", "Axes", "Units", "Description"],
        [[field * suffix, "scalar", "—", isempty(suffix) ? "model energy" : "model energy/site",
          "$label energy" * (isempty(suffix) ? "" : " per site")]
         for (label, field) in components for suffix in ("", "PerSite")])
    measurements = observable_table_lines(["Component", "Energy (model energy)",
                                           "Energy per site (model energy/site)"], rows)
    statistics = observable_production_lines(summary, "Energy")
    conventions = mean_observable_conventions(summary,
        "Extensive and per-site estimates are retained as stored, in model energy and model energy/site.")
    print_observable_report(io, "Energy", (
        ("OBSERVABLES", inventory), ("MEASUREMENTS", measurements),
        ("PRODUCTION STATISTICS", statistics), ("ESTIMATOR CONVENTIONS", conventions),
    ); width=maximum(textwidth, vcat(inventory, measurements, statistics)))
    return nothing
end

"""
    describe_windings(summary; io=stdout) -> nothing

Print one simulation's dimensionless winding means and squares by component and
spatial axis, saved totals, and available pair correlations/coflow/counterflow.
Rows use numeric subscripts (W₁, W₂, …) for component indices, not species IDs.
Headings use Unicode expectation brackets, squares, and sums; no LaTeX rendering
is required. Boundary conditions and pair mappings
are identified when available; otherwise explicitly marked unknown/unavailable.
An inventory uses mathematical labels rather than dictionary field names. Print
all estimates, component/axis information, production counts, and weighted-mean
and repeated-blocking conventions in aligned sections. No stiffness factors or
derived estimates are calculated. Single-component summaries have no component pairs.

Accept `missing` for unavailable measurements. Estimates must be Measurements in
component × axis or pair × axis matrices, with axis vectors for totals. Nonfinite
errors are shown as "uncertainty unavailable". The input is not modified.
"""
function describe_windings(summary::Union{Missing,AbstractDict}; io::IO=stdout)
    if ismissing(summary)
        println(io, "No Windings measurements available.")
        return nothing
    end
    means, squares = summary["windingBySpecies"], summary["squaredWindingBySpecies"]
    means isa AbstractMatrix && squares isa AbstractMatrix && size(means) == size(squares) ||
        throw(DimensionMismatch("Windings component means and squares must be matching matrices"))
    ncomponents, naxes = size(means)
    ncomponents > 0 && naxes > 0 || throw(DimensionMismatch("Windings component and axis dimensions must not be empty"))
    for field in ("totalWinding", "squaredTotalWinding")
        summary[field] isa AbstractVector && length(summary[field]) == naxes ||
            throw(DimensionMismatch("$field must have one entry per spatial axis"))
    end
    periodic = get(summary, "periodic_axes", nothing)
    has_periodic = periodic !== nothing && !ismissing(periodic)
    if has_periodic
        periodic isa AbstractVector{Bool} && length(periodic) == naxes ||
            throw(DimensionMismatch("periodic_axes must contain one Boolean per spatial axis"))
    end
    pairs = get(summary, "species_pairs", nothing)
    has_pairs = pairs !== nothing && !ismissing(pairs)
    has_pairs && validate_winding_pairs(pairs, ncomponents)
    fields = ("windingCorrelationBySpecies", "squaredCoflowBySpeciesPair", "squaredCounterflowBySpeciesPair")
    channels = [get(summary, field, nothing) for field in fields]
    available = findall(value -> value !== nothing && !ismissing(value), channels)
    for i in available
        channels[i] isa AbstractMatrix && size(channels[i], 2) == naxes ||
            throw(DimensionMismatch("$(fields[i]) must be a pair × spatial axis matrix"))
    end
    npairs = has_pairs ? size(pairs, 1) : isempty(available) ? 0 : size(channels[first(available)], 1)
    npairs <= ncomponents * (ncomponents - 1) ÷ 2 || throw(DimensionMismatch("Too many Windings component pairs"))
    all(i -> size(channels[i], 1) == npairs, available) ||
        throw(DimensionMismatch("Windings pair matrices and mapping must have matching row counts"))

    axes = ["Axis $a" for a in 1:naxes]
    component_labels = ["W" * join('₀' + (digit - '0') for digit in string(s)) for s in 1:ncomponents]
    component_rows = [[component_labels[s], axes[a], format_observable_estimate(means[s, a]),
                       format_observable_estimate(squares[s, a])]
                      for s in 1:ncomponents for a in 1:naxes]
    total_rows = [[axes[a], format_observable_estimate(summary["totalWinding"][a]),
                   format_observable_estimate(summary["squaredTotalWinding"][a])] for a in 1:naxes]
    pair_rows = Vector{String}[]
    for p in 1:npairs, a in 1:naxes
        label = has_pairs ? "($(component_labels[pairs[p, 1]]), $(component_labels[pairs[p, 2]]))" :
            "Pair row $p — component mapping unavailable"
        values = [i in available ? format_observable_estimate(channels[i][p, a]) : "Unavailable" for i in 1:3]
        push!(pair_rows, [label, axes[a], values...])
    end

    inventory_rows = [
        ["⟨Wₐ⟩", join(size(means), "×"), "component × axis", "dimensionless", "Component mean"],
        ["⟨Wₐ²⟩", join(size(squares), "×"), "component × axis", "dimensionless", "Component square"],
        ["⟨Σₐ Wₐ⟩", string(naxes), "axis", "dimensionless", "Total mean"],
        ["⟨(Σₐ Wₐ)²⟩", string(naxes), "axis", "dimensionless", "Total square"],
    ]
    for (i, (symbol, meaning)) in enumerate((("⟨Wₐ Wᵦ⟩", "Pair correlation"),
            ("⟨(Wₐ+Wᵦ)²⟩", "Pair coflow"), ("⟨(Wₐ−Wᵦ)²⟩", "Pair counterflow")))
        shape = i in available ? join(size(channels[i]), "×") : "Not present in this summary"
        push!(inventory_rows, [symbol, shape, "pair × axis", "dimensionless", meaning])
    end
    inventory = observable_table_lines(["Quantity", "Shape", "Axes", "Units", "Description"], inventory_rows)
    component_info = ["W₁, W₂, … label component indices; species IDs are not stored."]
    for a in 1:naxes
        boundary = has_periodic ? (periodic[a] ? "periodic" : "open") : "Boundary condition unknown"
        push!(component_info, "$(axes[a]): $boundary")
    end
    measurements = ["Per component"]
    append!(measurements, observable_table_lines(
        ["Component (a)", "Axis", "⟨Wₐ⟩", "⟨Wₐ²⟩"], component_rows))
    append!(measurements, ["", "Total over all components"])
    append!(measurements, observable_table_lines(["Axis", "⟨Σₐ Wₐ⟩", "⟨(Σₐ Wₐ)²⟩"], total_rows))
    append!(measurements, ["", "Per component pair"])
    if ncomponents == 1
        push!(measurements, "No component pairs")
    elseif npairs == 0
        push!(measurements, "Pair quantities: Unavailable")
    else
        append!(measurements, observable_table_lines(
            ["Pair (a,b)", "Axis", "⟨Wₐ Wᵦ⟩", "Coflow: ⟨(Wₐ+Wᵦ)²⟩",
             "Counterflow: ⟨(Wₐ−Wᵦ)²⟩"], pair_rows))
    end
    statistics = observable_production_lines(summary, "Windings")
    conventions = mean_observable_conventions(summary,
        "No stiffness normalization applied.")
    append!(conventions, [
        "  Total moments  Moments of Σₐ Wₐ on each spatial axis.",
        "  Pair moments   ⟨(Wₐ±Wᵦ)²⟩; no division by the number of components.",
    ])
    print_observable_report(io, "Windings (dimensionless)", (
        ("OBSERVABLES", inventory), ("COMPONENTS AND AXES", component_info),
        ("MEASUREMENTS", measurements), ("PRODUCTION STATISTICS", statistics),
        ("ESTIMATOR CONVENTIONS", conventions),
    ); width=maximum(textwidth, vcat(inventory, measurements, statistics)))
    return nothing
end

"""
    describe_particles_histogram(summary; io=stdout) -> nothing

Print one aggregated ParticlesHistogram inventory, bin metadata, all saved first
and second particle-number moments, a compact distribution overview, production
counts, and estimator conventions. Accept a summary dictionary or `missing`.
Rows retain species order; missing IDs use component labels. Missing metadata
is shown as "Unavailable"; nonfinite errors as "uncertainty unavailable".

The overview counts nonzero bins and shows their particle-number range, the
smallest particle number with the largest raw count (flagging tied peaks), and
its saved probability and standard error. It does not print full distributions.
Counts are raw sums; probabilities are counts divided by configuration count.
Moment units are particles and particles². Saved estimates and repeated-blocking
errors are displayed without recomputation. No files are read or input modified.
"""
function describe_particles_histogram(summary::Union{Missing,AbstractDict}; io::IO=stdout)
    if ismissing(summary)
        println(io, "No ParticlesHistogram measurements available.")
        return nothing
    end
    counts = summary["particleNumberBySpeciesCounts"]
    total_counts = summary["totalParticleNumberCounts"]
    probabilities = summary["particleNumberBySpeciesProbability"]
    total_probabilities = summary["totalParticleNumberProbability"]
    counts isa AbstractMatrix{<:Integer} && total_counts isa AbstractVector{<:Integer} &&
        !(eltype(counts) <: Bool) && !(eltype(total_counts) <: Bool) &&
        all(>=(0), counts) && all(>=(0), total_counts) ||
        throw(ArgumentError("Histogram counts must be nonnegative integer matrices/vectors"))
    nspecies, nbins = size(counts)
    nspecies > 0 && nbins > 0 && !isempty(total_counts) ||
        throw(DimensionMismatch("Histogram species and bin axes must not be empty"))
    probabilities isa AbstractMatrix && size(probabilities) == size(counts) &&
        total_probabilities isa AbstractVector && size(total_probabilities) == size(total_counts) ||
        throw(DimensionMismatch("Histogram probabilities must match their count arrays"))
    means = summary["meanParticleNumberBySpecies"]
    squares = summary["meanSquaredParticleNumberBySpecies"]
    means isa AbstractVector && squares isa AbstractVector &&
        length(means) == length(squares) == nspecies ||
        throw(DimensionMismatch("Histogram moments must have one entry per species"))

    metadata = [get(summary, key, nothing) for key in
        ("species", "species_particle_number_bins", "total_particle_number_bins", "max_particle_number_by_species")]
    for (i, expected) in enumerate((nspecies, nbins, length(total_counts), nspecies))
        value = metadata[i]
        if value === nothing || ismissing(value)
            metadata[i] = nothing
            continue
        end
        value isa AbstractVector{<:Integer} && !(eltype(value) <: Bool) ||
            throw(ArgumentError("Histogram labels and capacities must be integer vectors"))
        length(value) == expected || throw(DimensionMismatch("Histogram metadata does not match its axis"))
    end
    species, bins, total_bins, capacities = metadata
    species === nothing || (all(>(0), species) && allunique(species)) ||
        throw(ArgumentError("Histogram species IDs must be positive and unique"))
    for grid in (bins, total_bins)
        grid === nothing || grid == collect(0:(length(grid) - 1)) ||
            throw(ArgumentError("Histogram bins must enumerate consecutive particle numbers starting at zero"))
    end
    if capacities !== nothing
        all(capacity -> 0 <= capacity < nbins, capacities) ||
            throw(ArgumentError("Histogram capacities must fit the species bin axis"))
        all(s -> all(iszero, @view counts[s, (capacities[s] + 2):nbins]), 1:nspecies) ||
            throw(ArgumentError("Histogram counts exceed a species capacity"))
    end
    labels = [species === nothing ? "Component $s" : "Species $(species[s])" for s in 1:nspecies]
    inventory_rows = Vector{String}[]
    for (field, axes, units, meaning) in (
        ("particleNumberBySpeciesCounts", "species × bin", "counts", "Raw species counts"),
        ("totalParticleNumberCounts", "bin", "counts", "Raw total counts"),
        ("particleNumberBySpeciesProbability", "species × bin", "dimensionless", "Species probabilities"),
        ("totalParticleNumberProbability", "bin", "dimensionless", "Total probabilities"),
        ("meanParticleNumberBySpecies", "species", "particles", "⟨Nₛ⟩"),
        ("meanSquaredParticleNumberBySpecies", "species", "particles²", "⟨Nₛ²⟩"),
        ("meanTotalParticleNumber", "—", "particles", "⟨Ntotal⟩"),
        ("meanSquaredTotalParticleNumber", "—", "particles²", "⟨Ntotal²⟩"),
        ("zSectorSampleCount", "—", "configurations", "Configuration count"),
    )
        value = summary[field]
        shape = value isa AbstractArray ? join(size(value), "×") : "scalar"
        push!(inventory_rows, [field, shape, axes, units, meaning])
    end
    inventory = observable_table_lines(["Quantity", "Shape", "Axes", "Units", "Description"], inventory_rows)
    bin_range(grid) = grid === nothing ? "Unavailable" : "$(first(grid))–$(last(grid))"
    bin_rows = [[labels[s], format_observable_count(nbins), bin_range(bins),
                 capacities === nothing ? "Unavailable" : string(capacities[s])] for s in 1:nspecies]
    push!(bin_rows, ["Total", format_observable_count(length(total_counts)), bin_range(total_bins), "—"])
    bin_info = observable_table_lines(["Distribution", "Stored bins", "Stored N range", "Species capacity"], bin_rows)
    species === nothing && push!(bin_info, "Species IDs unavailable; rows use component indices.")
    push!(bin_info, "Bins above a species capacity are padding; zero-count allowed bins were not observed.")

    moment_rows = [[labels[s], format_observable_estimate(means[s]), format_observable_estimate(squares[s])]
                   for s in 1:nspecies]
    push!(moment_rows, ["Total", format_observable_estimate(summary["meanTotalParticleNumber"]),
                       format_observable_estimate(summary["meanSquaredTotalParticleNumber"])])
    moments = observable_table_lines(["Distribution", "⟨N⟩ (particles)", "⟨N²⟩ (particles²)"], moment_rows)
    overview_rows = [particle_histogram_overview_row(labels[s], @view(counts[s, :]),
                    @view(probabilities[s, :]), bins) for s in 1:nspecies]
    push!(overview_rows, particle_histogram_overview_row("Total", total_counts, total_probabilities, total_bins))
    overview = observable_table_lines(
        ["Distribution", "Occupied bins", "Observed N range", "Peak N", "P(peak)"], overview_rows)
    overview_width = maximum(textwidth, overview)
    append!(overview, [
        "Occupied bins have nonzero raw counts; the observed range may contain empty bins.",
        "Tied peaks share the largest raw count; the smallest N is shown.",
        "Full distributions remain in the summary count and probability arrays.",
    ])
    statistics = observable_table_lines(["", ""], [
        ["Saved blocks", format_observable_count(get(summary, "nblocks", nothing))],
        ["ParticlesHistogram samples", format_observable_count(get(summary, "nsamples", nothing))],
        ["Z-sector configurations", format_observable_count(summary["zSectorSampleCount"])],
    ]; right_columns=(2,))[2:end]
    saved(key) = get(summary, key, nothing) === nothing || ismissing(summary[key]) ? "Unavailable" : summary[key]
    conventions = [
        "  Raw aggregation  " * observable_method(saved("aggregation_method"), "sum", "Sum counts across saved blocks; no sample-count multiplication"),
        "  Probability      Raw bin count / nsamples (configuration count)",
        "  Moments          Σ Nᵖ counts / nsamples, p=1,2; total moments use the saved total histogram",
        "  Uncertainty      " * observable_method(saved("error_method"),
            "repeated blocking: maximum standard error", "Max repeated-blocking SE; sample-count-weighted normalized block histories"),
        "  Raw-count errors None stored",
        "  Second moments   ⟨N²⟩, not ⟨N⟩²; estimates and uncertainties are retained as stored",
    ]
    print_observable_report(io, "ParticlesHistogram", (
        ("OBSERVABLES", inventory), ("SPECIES AND BINS", bin_info), ("MOMENTS", moments),
        ("DISTRIBUTION OVERVIEW", overview), ("PRODUCTION STATISTICS", statistics),
        ("ESTIMATOR CONVENTIONS", conventions),
    ); width=max(overview_width, maximum(textwidth, vcat(inventory, bin_info, moments, statistics))))
    return nothing
end

"""Summarize occupied bins and choose the first maximum on the ordered bin axis."""
function particle_histogram_overview_row(label, counts, probabilities, bins)
    occupied = count(>(0), counts)
    occupied == 0 && return [label, "0", "Unavailable", "Unavailable", "Unavailable"]
    peak = argmax(counts)
    ties = count(==(counts[peak]), counts)
    observed = bins === nothing ? "Unavailable" : "$(bins[findfirst(>(0), counts)])–$(bins[findlast(>(0), counts)])"
    peak_label = bins === nothing ? "Unavailable" : string(bins[peak])
    ties > 1 && (peak_label *= " (tie: $ties bins)")
    return [label, format_observable_count(occupied), observed, peak_label,
            format_observable_estimate(probabilities[peak])]
end

"""
    describe_single_particle_green(summary; io=stdout) -> nothing

Describe one simulation's retained `gKTau`, `grEqualTime`, and `gREqualTimeHits`
in an aligned report containing an observable inventory, coordinate grids, first
stored examples, production statistics, and estimator conventions. Complex
`gKTau` errors are reported separately for real and imaginary parts.

Absent quantities say "Not present in this summary", since older aggregates may
omit measured data. Missing metadata is "Unavailable". Accept `missing` for an
unavailable summary. Only the supplied dictionary is inspected; no files are read,
values normalized, or input data modified. Nonfinite errors are displayed as
"uncertainty unavailable". Momentum/displacement grids show first/last rows;
consecutive `tau_edges` bound each τ bin. Return `nothing`.
"""
function describe_single_particle_green(summary::Union{Missing,AbstractDict}; io::IO=stdout)
    if ismissing(summary)
        println(io, "No SingleParticleGreenFunction measurements available.")
        return nothing
    end
    saved(key) = haskey(summary, key) && summary[key] !== nothing && !ismissing(summary[key]) ?
        summary[key] : "Unavailable"
    quantities = (("gKTau", 3, "species × k × τ", "Raw weighted mean"),
                  ("grEqualTime", 2, "species × r", "Normalized estimate"),
                  ("gREqualTimeHits", 2, "species × r", "Raw hits"))
    present = [haskey(summary, key) && summary[key] !== nothing && !ismissing(summary[key])
               for (key, _, _, _) in quantities]
    rows = Vector{String}[]
    for (i, (key, rank, axes, meaning)) in enumerate(quantities)
        if present[i]
            data = summary[key]
            data isa AbstractArray && ndims(data) == rank ||
                throw(DimensionMismatch("$key must have $rank dimensions ($axes)"))
            push!(rows, [key, join(size(data), "×"), axes, meaning])
        else
            push!(rows, [key, "Not present in this summary", axes, "—"])
        end
    end
    species = get(summary, "species", nothing)
    species = ismissing(species) ? nothing : species
    if species !== nothing
        species isa AbstractVector || throw(ArgumentError("species must be a vector"))
        all(i -> !present[i] || size(summary[quantities[i][1]], 1) == length(species), 1:3) ||
            throw(DimensionMismatch("Species metadata does not match the measurement arrays"))
    end
    nk = present[1] ? size(summary["gKTau"], 2) : nothing
    nr = present[2] ? size(summary["grEqualTime"], 2) : present[3] ? size(summary["gREqualTimeHits"], 2) : nothing
    if present[2] && present[3]
        size(summary["grEqualTime"]) == size(summary["gREqualTimeHits"]) ||
            throw(DimensionMismatch("Equal-time estimates and hits must have matching axes"))
    end
    edges = get(summary, "tau_edges", nothing)
    edges = ismissing(edges) ? nothing : edges
    if edges !== nothing
        edges isa AbstractVector{<:Real} && length(edges) >= 2 && all(isfinite, edges) &&
            all(>(0), diff(edges)) || throw(ArgumentError("tau_edges must be finite and strictly increasing"))
        !present[1] || length(edges) == size(summary["gKTau"], 3) + 1 ||
            throw(DimensionMismatch("tau_edges does not match gKTau time bins"))
    end

    momenta = green_coordinate_grid(summary, "momenta", nk)
    displacements = green_coordinate_grid(summary, "displacements", nr)
    displacement_vectors = green_coordinate_grid(summary, "displacement_vectors", nr)

    observable_lines = observable_table_lines(
        ["Quantity", "Shape", "Axes", "Description"], rows)
    push!(observable_lines, "Axes: k = momentum, r = displacement, τ = imaginary time.")

    coordinate_rows = Vector{String}[]
    if species === nothing
        push!(coordinate_rows, ["Species IDs", "Unavailable", "Unavailable", "Unavailable"])
    elseif isempty(species)
        push!(coordinate_rows, ["Species IDs", "0", "Unavailable", "Unavailable"])
    else
        push!(coordinate_rows, ["Species IDs", string(length(species)), string(first(species)), string(last(species))])
    end
    for (label, grid) in (("Momentum", momenta), ("Displacement", displacements),
                          ("Displacement vectors", displacement_vectors))
        if grid === nothing
            push!(coordinate_rows, [label, "Unavailable", "Unavailable", "Unavailable"])
        elseif isempty(grid)
            push!(coordinate_rows, [label, "0", "Unavailable", "Unavailable"])
        else
            push!(coordinate_rows, [label, string(size(grid, 1)),
                                    format_green_coordinate_row(grid, 1),
                                    format_green_coordinate_row(grid, size(grid, 1))])
        end
    end
    if edges === nothing
        push!(coordinate_rows, ["τ bins", "Unavailable", "Unavailable", "Unavailable"])
    else
        push!(coordinate_rows, ["τ bins", string(length(edges) - 1),
                                format_green_coordinate(first(edges)),
                                format_green_coordinate(last(edges))])
    end
    coordinate_lines = observable_table_lines(
        ["Grid", "Count", "First", "Last"], coordinate_rows; right_columns=(2,))
    append!(coordinate_lines, [
        "Note: Coordinate grids do not imply measurement availability.",
        "Full grids: summary[\"species\"], summary[\"momenta\"], summary[\"displacements\"], " *
        "summary[\"displacement_vectors\"], and summary[\"tau_edges\"].",
    ])

    example_lines = String[]
    for (i, (key, _, _, _)) in enumerate(quantities)
        present[i] || continue
        data = summary[key]
        if isempty(data)
            push!(example_lines, "$key: Unavailable (empty array)")
            continue
        end
        indices = join(Tuple(first(CartesianIndices(data))), ",")
        species_label = species === nothing ? "Unavailable" : string(first(species))
        if key == "gKTau"
            momentum = momenta === nothing || isempty(momenta) ? "Unavailable" :
                format_green_coordinate_row(momenta, 1)
            tau = edges === nothing ? "Unavailable" :
                "[" * join(format_green_coordinate.(edges[1:2]), ", ") * "]"
            push!(example_lines, "$key[$indices]    species=$species_label  k=$momentum  τ=$tau")
            value = first(data)
            value isa Complex || throw(ArgumentError("gKTau estimates must be complex Measurements"))
            append!(example_lines, indent_green_example_rows([
                green_estimate_row("Real", real(value)),
                green_estimate_row("Imaginary", imag(value)),
            ]))
        else
            displacement = displacements === nothing || isempty(displacements) ? "Unavailable" :
                format_green_coordinate_row(displacements, 1)
            push!(example_lines, "$key[$indices]    species=$species_label  r=$displacement")
            if key == "grEqualTime"
                value_row = green_estimate_row("Value", first(data))
                value_row[4] = "error status: " * string(saved("error_status"))
                append!(example_lines, indent_green_example_rows([
                    value_row,
                ]))
            else
                first(data) isa Real && !(first(data) isa Measurements.Measurement) ||
                    throw(ArgumentError("gREqualTimeHits must contain raw numeric counts"))
                append!(example_lines, indent_green_example_rows([
                    ["Hits", format_observable_count(first(data)), "", "raw"],
                ]))
            end
        end
        push!(example_lines, "")
    end
    isempty(example_lines) && push!(example_lines, "No retained measurement arrays.")
    push!(example_lines, "Note: Examples show first stored entries, not typical values.")

    statistics_lines = observable_table_lines(["", ""], [
        ["Saved blocks", format_observable_count(saved("nblocks"))],
        ["Green-function samples", format_observable_count(saved("sampleCount"))],
        ["Z-sector visits", format_observable_count(saved("zSectorVisitCount"))],
    ]; right_columns=(2,))[2:end]

    convention_lines = String[]
    if present[1]
        append!(convention_lines, ["gKTau",
            "  Aggregation    " * observable_method(saved("gKTauAggregationMethod"),
                "sample-count-weighted mean", "sampleCount-weighted block means"),
            "  Uncertainty    " * observable_method(saved("gKTauErrorMethod"),
                "repeated blocking: maximum standard error", "Max repeated-blocking SE (Re/Im separately)"),
            "  Normalization  None; no Z-visit division or physical normalization", ""])
    end
    if present[2]
        append!(convention_lines, ["grEqualTime",
            "  Uncertainty    " * observable_method(saved("error_method"),
                "repeated blocking: maximum paired jackknife standard error",
                "Max repeated-blocking paired jackknife SE"),
            "  Window         " * string(saved("equal_time_window")),
            "  Convention     " * string(saved("equal_time_convention")),
            "  Normalization  " * string(saved("normalization")), ""])
    end
    if present[3]
        append!(convention_lines, ["gREqualTimeHits", "  Uncertainty    None (raw hit counts)"])
    end
    isempty(convention_lines) && push!(convention_lines, "No estimator conventions available.")
    !isempty(convention_lines) && isempty(last(convention_lines)) && pop!(convention_lines)

    sections = (("OBSERVABLES", observable_lines), ("COORDINATE GRIDS", coordinate_lines),
                ("MEASUREMENT EXAMPLES", example_lines), ("PRODUCTION STATISTICS", statistics_lines),
                ("ESTIMATOR CONVENTIONS", convention_lines))
    width_candidates = vcat(observable_lines[1:length(rows) + 1],
                            coordinate_lines[1:length(coordinate_rows) + 1], statistics_lines)
    print_observable_report(io, "SingleParticleGreenFunction", sections;
                            width=maximum(textwidth, width_candidates))
    return nothing
end

"""Render report sections with a shared width and wrapped explanatory text."""
function print_observable_report(io, title, sections; width)
    report_width = max(68, width)
    println(io, title)
    println(io, repeat("═", report_width))
    for (section_index, (heading, lines)) in enumerate(sections)
        section_index == 1 || println(io)
        println(io, heading)
        println(io, repeat("─", report_width))
        for line in lines
            print_observable_wrapped_line(io, line, report_width)
        end
    end
    return nothing
end

"""Validate and return one saved Green-function coordinate grid, when available."""
function green_coordinate_grid(summary, key, expected_count)
    grid = get(summary, key, nothing)
    if grid === nothing || ismissing(grid)
        return nothing
    end
    grid isa AbstractMatrix{<:Real} && size(grid, 2) > 0 && all(isfinite, grid) ||
        throw(ArgumentError("$key must be a finite numeric coordinate matrix"))
    expected_count === nothing || size(grid, 1) == expected_count ||
        throw(DimensionMismatch("$key coordinates do not match the measurement axis"))
    return grid
end

function format_green_coordinate(value::Real)
    value isa Integer && return string(value)
    numeric = Float64(value)
    if numeric != 0 && (abs(numeric) < 0.00005 || abs(numeric) >= 1.0e6)
        return Printf.@sprintf("%.4e", numeric)
    end
    return Printf.@sprintf("%.4f", numeric)
end

format_green_coordinate_row(grid, row) =
    "[" * join(format_green_coordinate.(grid[row, :]), ", ") * "]"

function format_observable_count(value)
    (value === nothing || ismissing(value)) && return "Unavailable"
    if value isa Integer
        text = string(value)
    elseif value isa Real && isfinite(value) && isinteger(value)
        text = Printf.@sprintf("%.0f", Float64(value))
    else
        return string(value)
    end
    return replace(text, r"(?<=\d)(?=(?:\d{3})+$)" => ",")
end

observable_method(value, stored, display) =
    value isa AbstractString && value == stored ? display : string(value)

function green_estimate_row(label, value)
    parts = split(format_observable_estimate(value), " ± "; limit=2)
    return [label, first(parts), length(parts) == 2 ? "± " * last(parts) : "", ""]
end

function indent_green_example_rows(rows)
    lines = observable_table_lines(["", "", "", ""], rows; right_columns=(2, 3))[2:end]
    return ["  " * line for line in lines]
end

"""Print one report line, wrapping without truncation at spaces or formula operators."""
function print_observable_wrapped_line(io, line, width)
    if textwidth(line) <= width
        println(io, line)
        return nothing
    end
    prefix_match = match(r"^(\s*\S+(?:\s+\S+)*?\s{2,})", line)
    continuation = repeat(" ", prefix_match === nothing ? 2 : textwidth(prefix_match.captures[1]))
    remaining, prefix = line, ""
    while !isempty(remaining)
        available = width - textwidth(prefix)
        characters = collect(remaining)
        cumulative = 0
        cut = 0
        for (i, character) in enumerate(characters)
            next_width = cumulative + textwidth(string(character))
            next_width > available && break
            cumulative = next_width
            cut = i
        end
        cut > 0 || throw(ArgumentError("report width is too narrow"))
        if cut < length(characters)
            space = findlast(isspace, characters[1:cut])
            if space !== nothing && space > 1
                cut = space - 1
            else
                operator = findlast(character -> character in ('*', '/'), characters[1:cut])
                operator !== nothing && (cut = operator)
            end
        end
        println(io, prefix, rstrip(String(characters[1:cut])))
        remaining = lstrip(String(characters[cut + 1:end]))
        prefix = continuation
    end
    return nothing
end

"""Format a stored estimate without interpreting nonfinite errors as usable uncertainties."""
function format_observable_estimate(value)
    value isa Measurements.Measurement || throw(ArgumentError("Observable estimates must be Measurements"))
    return isfinite(Measurements.uncertainty(value)) ? string(value) :
        "$(Measurements.value(value)) (uncertainty unavailable)"
end

"""Return aligned table lines, optionally right-aligning selected columns."""
function observable_table_lines(headers, rows; right_columns=())
    widths = [maximum(textwidth, [headers[j]; [row[j] for row in rows]]) for j in eachindex(headers)]
    return [rstrip(join([j in right_columns ? lpad(row[j], widths[j]) : rpad(row[j], widths[j])
                         for j in eachindex(headers)], "  ")) for row in [headers, rows...]]
end

function observable_production_lines(summary, observable)
    return observable_table_lines(["", ""], [
        ["Saved blocks", format_observable_count(get(summary, "nblocks", nothing))],
        ["$observable samples", format_observable_count(get(summary, "nsamples", nothing))],
    ]; right_columns=(2,))[2:end]
end

"""Describe the common block-mean estimator using the saved uncertainty method."""
function mean_observable_conventions(summary, normalization)
    method = get(summary, "error_method", nothing)
    uncertainty = method === nothing || ismissing(method) ? "Unavailable" :
        observable_method(method, "repeated blocking: maximum standard error",
                          "Max repeated-blocking SE")
    return [
        "  Aggregation    sampleCount-weighted block means",
        "  Uncertainty    " * uncertainty,
        "  Normalization  " * normalization,
    ]
end
