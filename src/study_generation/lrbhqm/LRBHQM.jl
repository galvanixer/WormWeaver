module LRBHQM

import DataFrames
import JLD2
import YAML

import ..Study
import ..WormWeaver: inspect_g2_entry_modes

export generate_study

include("unit_density.jl")

const UNIT_DENSITY_DATA_DIR = normpath(joinpath(
    @__DIR__, "..", "..", "..", "data", "lrbhqm", "unit_density"))
const MODEL_COLUMNS = (
    :geometry, :model_family, :ensemble, :nspecies, :hopping_kind,
    :hopping_alpha_1, :hopping_alpha_2, :hopping_distance_cutoff_1,
    :hopping_distance_cutoff_2, :U_1_1, :U_1_2, :U_2_1, :U_2_2,
)

function _data_path(name)
    name isa AbstractString && !isempty(strip(name)) ||
        throw(ArgumentError("unit_density_data must be a nonempty bundled filename"))
    basename(name) == name && !occursin('\0', name) || throw(ArgumentError(
        "unit_density_data must be a filename without directory components"))
    endswith(name, ".jld2") ||
        throw(ArgumentError("unit_density_data must name a .jld2 file"))
    path = joinpath(UNIT_DENSITY_DATA_DIR, name)
    isfile(path) || throw(ArgumentError("Unknown unit-density dataset: $name"))
    return path
end

function _load_crossings(path)
    data = JLD2.load(path, "curated_crossings")
    data isa DataFrames.AbstractDataFrame ||
        throw(ArgumentError("curated_crossings must be a DataFrame"))
    required = [:model_id, MODEL_COLUMNS..., :g2_mode, :nsites, :beta,
        :t_over_U, :mu_at_unit_density, :mu_lower, :mu_upper, :include]
    absent = setdiff(required, propertynames(data))
    isempty(absent) || throw(ArgumentError(
        "Curated dataset is missing columns: $(join(absent, ", "))"))
    isempty(data) && throw(ArgumentError("Curated dataset must not be empty"))
    all(value -> value isa Bool, data.include) ||
        throw(ArgumentError("Curated include values must be Boolean"))
    all(value -> value isa AbstractString && !isempty(value), data.model_id) ||
        throw(ArgumentError("Curated model IDs must be nonempty strings"))
    all(value -> value isa AbstractString && !isempty(value), data.g2_mode) ||
        throw(ArgumentError("Curated lookup modes must be nonempty strings"))
    all(value -> value isa Integer && !(value isa Bool) && value > 0, data.nsites) ||
        throw(ArgumentError("Curated lattice sizes must be positive integers"))
    all(value -> value isa Real && isfinite(value), data.beta) ||
        throw(ArgumentError("Curated beta values must be finite"))
    all(value -> value isa Real && isfinite(value), data.t_over_U) ||
        throw(ArgumentError("Curated t/U values must be finite"))
    length(unique(data.model_id)) == 1 ||
        throw(ArgumentError("Curated dataset must contain exactly one model_id"))
    reference = first(eachrow(data))
    for row in eachrow(data), column in MODEL_COLUMNS
        _matches(row[column], reference[column]) || throw(ArgumentError(
            "Curated dataset must contain exactly one physical model"))
    end
    for row in eachrow(DataFrames.filter(:include => identity, data))
        all(column -> row[column] isa Real && isfinite(row[column]),
            (:mu_at_unit_density, :mu_lower, :mu_upper)) || throw(ArgumentError(
            "Included crossings must have finite chemical-potential estimates"))
    end
    keys = [:model_id, :g2_mode, :nsites, :beta, :t_over_U]
    any(DataFrames.nonunique(data, keys)) &&
        throw(ArgumentError("Curated dataset contains duplicate crossing keys"))
    return data
end

function _study_fields(study)
    fields = String[]
    for section in ("overrides", "sweep")
        values = get(study, section, Dict())
        values isa AbstractDict || continue
        append!(fields, String.(keys(values)))
    end
    groups = get(study, "correlated_sweeps", [])
    groups isa AbstractVector || return fields
    for group in groups
        group isa AbstractDict && append!(fields, String.(keys(group)))
    end
    return fields
end

function _reject_mu_settings(study)
    target = "model.onsite.mu"
    for field in _study_fields(study)
        if field == target || startswith(target, field * ".") ||
           startswith(field, target * ".")
            throw(ArgumentError(
                "$field conflicts with unit-density chemical-potential resolution"))
        end
    end
end

function _required(config, path)
    value = Study.get_field(config, path)
    ismissing(value) && throw(ArgumentError("Configuration must define $path"))
    return value
end

function _finite_vector(value, length_, label)
    value isa AbstractVector && length(value) == length_ &&
        all(item -> item isa Real && isfinite(item), value) ||
        throw(ArgumentError("$label must contain $length_ finite numbers"))
    return value
end

function _model_parameters(config)
    geometry = _required(config, "lattice.geometry")
    geometry == "chain" || throw(ArgumentError("LRBHQM requires lattice.geometry: chain"))
    extent = _required(config, "lattice.extent")
    extent isa AbstractVector && length(extent) == 1 && only(extent) isa Integer &&
        !(only(extent) isa Bool) && only(extent) > 0 || throw(ArgumentError(
        "LRBHQM lattice.extent must contain one positive integer"))
    nspecies = _required(config, "model.nspecies")
    nspecies == 2 || throw(ArgumentError("Curated LRBHQM data require two species"))
    hopping = _finite_vector(_required(config, "model.hopping.J"), nspecies,
        "model.hopping.J")
    _matches(hopping[1], hopping[2]) ||
        throw(ArgumentError("LRBHQM unit-density lookup requires equal hoppings"))
    alpha = _finite_vector(_required(config, "model.hopping.alpha"), nspecies,
        "model.hopping.alpha")
    cutoffs = _required(config, "model.hopping.distance_cutoff")
    cutoffs isa AbstractVector && length(cutoffs) == nspecies || throw(ArgumentError(
        "model.hopping.distance_cutoff must contain one value per species"))
    onsite = _required(config, "model.onsite.U")
    onsite isa AbstractVector && length(onsite) == nspecies &&
        all(row -> row isa AbstractVector && length(row) == nspecies &&
            all(value -> value isa Real && isfinite(value), row), onsite) ||
        throw(ArgumentError("model.onsite.U must be a finite 2×2 matrix"))
    U11 = onsite[1][1]
    U11 > 0 || throw(ArgumentError("model.onsite.U[1][1] must be positive"))
    beta = _required(config, "system.beta")
    beta isa Real && isfinite(beta) ||
        throw(ArgumentError("system.beta must be finite"))
    model = (
        geometry=geometry,
        model_family=_required(config, "model.family"),
        ensemble=_required(config, "system.ensemble"),
        nspecies=nspecies,
        hopping_kind=_required(config, "model.hopping.kind"),
        hopping_alpha_1=alpha[1],
        hopping_alpha_2=alpha[2],
        hopping_distance_cutoff_1=cutoffs[1],
        hopping_distance_cutoff_2=cutoffs[2],
        U_1_1=onsite[1][1],
        U_1_2=onsite[1][2],
        U_2_1=onsite[2][1],
        U_2_2=onsite[2][2],
    )
    return (; model, nsites=only(extent), beta, t_over_U=hopping[1] / U11,
        U11, nspecies)
end

function _lookup_mode(config, match_sampling_mode)
    match_sampling_mode || return "none"
    modes = inspect_g2_entry_modes(config)
    length(modes) <= 1 || throw(ArgumentError(
        "Sampling-mode matching requires at most one scheduled G2 entry mode"))
    return isempty(modes) ? "none" : string(only(modes))
end

function _resolver(data, data_name, match_sampling_mode)
    return function (config)
        parameters = _model_parameters(config)
        mode = _lookup_mode(config, match_sampling_mode)
        estimate = estimate_mu_over_U(data; parameters.model, g2_mode=mode,
            parameters.nsites, parameters.beta, parameters.t_over_U)
        mu = parameters.U11 * estimate.mu_over_U
        isfinite(mu) || throw(ArgumentError("Resolved chemical potential must be finite"))
        Study.set_field!(config, "model.onsite.mu", fill(mu, parameters.nspecies))
        manifest = Dict{String,Any}(
            "lrbhqm.unit_density_data" => data_name,
            "lrbhqm.model_id" => estimate.model_id,
            "lrbhqm.lookup_mode" => mode,
            "lrbhqm.match_sampling_mode" => match_sampling_mode,
            "lrbhqm.t_over_U" => estimate.t_over_U,
            "lrbhqm.requested_beta" => estimate.requested_beta,
            "lrbhqm.selected_beta" => estimate.beta,
            "lrbhqm.mu_over_U" => estimate.mu_over_U,
            "lrbhqm.estimation_method" => string(estimate.method),
            "lrbhqm.supporting_t_over_U" => collect(estimate.supporting_t_over_U),
        )
        return (; config, manifest)
    end
end

"""
    generate_study(path; overwrite=false, validate_config=nothing,
                   match_sampling_mode=false) -> String

Generate an LRBHQM study at unit density using the bundled JLD2 file named by
`unit_density_data`. Physical hoppings remain ordinary configuration values;
the generator derives `t/U₁₁` and resolves equal species chemical potentials.
By default, use crossings with lookup mode `none`; enable `match_sampling_mode`
to select crossings from the configured G2 entry mode.
"""
function generate_study(path::AbstractString; overwrite::Bool=false,
                        validate_config=nothing, match_sampling_mode::Bool=false)
    source = abspath(path)
    study = Study.mapping(YAML.load(read(source, String)), "Study")
    data_name = get(study, "unit_density_data", nothing)
    data_path = _data_path(data_name)
    _reject_mu_settings(study)
    data = _load_crossings(data_path)
    resolver = _resolver(data, String(data_name), match_sampling_mode)
    return Study._generate_study(source; overwrite, validate_config, resolve_config=resolver,
        extra_allowed=("unit_density_data",))
end

end
