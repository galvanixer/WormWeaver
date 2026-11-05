# Aggregate saved raw static structure factors

Aggregate the saved `StaticStructureFactor` block histories directly, preserving
their site-normalized, imaginary-time-averaged equal-time convention and Fourier
grid. Reconstructing from density-correlation summaries would discard the block
associations needed for statistical errors; connected subtraction also requires
Fourier density means that these histories do not provide.

Analyze the saved total-density history independently so cancellations between
component channels contribute correctly to its uncertainty. Each stored complex
pair entry receives separate real and imaginary repeated-blocking errors, matching
the existing Green-function representation. This retains marginal uncertainties
without a covariance matrix; derived channels must use the original joint block
histories for covariance-aware analysis. Do not symmetrize or infer physical
species IDs. Preserve native field names and append the observable to defaults;
callers requiring a fixed output schema can select observables explicitly.
