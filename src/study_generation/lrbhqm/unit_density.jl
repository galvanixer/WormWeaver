_absent(value) = value === nothing || ismissing(value) || value == "null"

function _matches(a, b)
    _absent(a) && _absent(b) && return true
    a isa Real && b isa Real && return isapprox(a, b; rtol=1e-8, atol=1e-12)
    return isequal(a, b) || (a isa AbstractString && b isa AbstractString && a == b)
end

"""
    estimate_mu_over_U(data; model, g2_mode="none", nsites, beta, t_over_U)

Return a sampled, interpolated, or lower-side extrapolated unit-density chemical
potential from a curated crossings table. Only rows with `include == true` are
anchors. An unavailable β and lower-side extrapolation use the largest available
β. Upper-side extrapolation is unsupported.
"""
function estimate_mu_over_U(data; model, g2_mode="none", nsites, beta, t_over_U)
    model isa NamedTuple || throw(ArgumentError("model must be a named tuple"))
    required = [:model_id, :g2_mode, :nsites, :beta, :t_over_U,
        :mu_at_unit_density, :mu_lower, :mu_upper, :include]
    absent = setdiff([required; collect(keys(model))], propertynames(data))
    isempty(absent) || throw(ArgumentError("Missing columns: $(join(absent, ", "))"))
    beta isa Real && isfinite(beta) || throw(ArgumentError("beta must be finite"))
    t_over_U isa Real && isfinite(t_over_U) ||
        throw(ArgumentError("t_over_U must be finite"))

    candidates = DataFrames.filter(row -> row.g2_mode == g2_mode &&
        row.nsites == nsites &&
        all(column -> _matches(row[column], model[column]), keys(model)), data)
    isempty(candidates) &&
        throw(ArgumentError("No curated data match the requested model, mode, and L"))
    length(unique(candidates.model_id)) == 1 ||
        throw(ArgumentError("Model parameters match multiple model IDs"))
    available_betas = unique(candidates.beta)
    all(value -> value isa Real && isfinite(value), available_betas) ||
        throw(ArgumentError("Available beta values must be finite"))
    matching_betas = filter(value -> _matches(value, beta), available_betas)
    length(matching_betas) <= 1 ||
        throw(ArgumentError("Requested beta matches multiple stored values"))
    selected_beta = isempty(matching_betas) ? maximum(available_betas) : only(matching_betas)
    selected = DataFrames.filter(row -> isequal(row.beta, selected_beta), candidates)
    length(unique(selected.t_over_U)) == DataFrames.nrow(selected) ||
        throw(ArgumentError("Duplicate t/U values in the selected curve"))
    points = sort(DataFrames.filter(:include => identity, selected), :t_over_U)
    largest_beta = maximum(available_betas)
    if (isempty(points) || t_over_U < first(points.t_over_U)) &&
       !_matches(selected_beta, largest_beta)
        selected_beta = largest_beta
        selected = DataFrames.filter(row -> isequal(row.beta, selected_beta), candidates)
        length(unique(selected.t_over_U)) == DataFrames.nrow(selected) ||
            throw(ArgumentError("Duplicate t/U values in the selected curve"))
        points = sort(DataFrames.filter(:include => identity, selected), :t_over_U)
    end
    isempty(points) && throw(ArgumentError("No included crossings remain"))
    xs = points.t_over_U
    sampled = findall(x -> _matches(x, t_over_U), xs)

    if length(sampled) == 1
        left = right = only(sampled)
        weight, method = 0.0, :sampled
    elseif first(xs) < t_over_U < last(xs)
        right = searchsortedfirst(xs, t_over_U)
        left = right - 1
        weight, method = (t_over_U - xs[left]) / (xs[right] - xs[left]), :interpolated
    elseif t_over_U < first(xs)
        DataFrames.nrow(points) >= 2 ||
            throw(ArgumentError("Extrapolation requires at least two included crossings"))
        left, right = 1, 2
        weight, method = (t_over_U - xs[left]) / (xs[right] - xs[left]), :extrapolated
    else
        throw(ArgumentError("Upper-side extrapolation is not supported"))
    end

    interpolate(column) = left == right ? points[left, column] :
        (1 - weight) * points[left, column] + weight * points[right, column]
    mu = interpolate(:mu_at_unit_density)
    if method == :extrapolated
        mu_lower = (1 - weight) * points.mu_lower[left] + weight * points.mu_upper[right]
        mu_upper = (1 - weight) * points.mu_upper[left] + weight * points.mu_lower[right]
    else
        mu_lower, mu_upper = interpolate(:mu_lower), interpolate(:mu_upper)
    end

    return (; model_id=only(unique(selected.model_id)), nsites,
        requested_beta=beta, beta=selected_beta, t_over_U, mu_over_U=mu,
        mu_lower, mu_upper, method, supporting_t_over_U=(xs[left], xs[right]))
end
