# WormWeaver

Terms for combining and interpreting quantum Monte Carlo measurements.

## Language

**Scheduled G2 entry**:
A positive-weight opening from a one-worm sector reachable from Z into a two-worm sector. It describes access permitted by the move schedule, rather than an observed or accepted transition.

**G2 entry mode**:
The proposal type permitted for scheduled G2 entry: regular, colocated aligned, or colocated anti-aligned. A schedule may permit multiple modes; an unrestricted colocated opening permits both alignments.

**Simulation aggregate**:
A summary of production measurements from one simulation run, retaining observable estimates and their within-run statistical uncertainties. Separate runs remain separate summaries until combined in a subsequent analysis.

**Density report**:
A readable account of one simulation's aggregated species densities, species particle numbers, and total particle number, with units and uncertainties. It includes the saved production block count, Density sample count, and uncertainty method.

**Raw density correlation**:
The imaginary-time average of an equal-time occupation product, averaged over site pairs with the same signed displacement. It retains disconnected contributions and the same-site occupation square; the total-density correlation includes all ordered species pairs.
_Avoid_: Connected density correlation, time-separated susceptibility

**Raw static structure factor**:
The site-normalized imaginary-time average of an equal-time component-density Fourier product. It retains disconnected contributions; the total-density channel includes all ordered component pairs.
_Avoid_: Connected structure factor, density susceptibility

**Component-pair static structure factor**:
The ordered-component matrix of raw static structure factors. Off-diagonal entries may be complex, with reversing the component order giving the complex conjugate.

**Total-density static structure factor**:
The real static structure factor of the total component density, equal to the sum over all ordered component-pair channels.

**Energy report**:
A readable account of one simulation's diagonal, kinetic, and total energy estimates, both extensive and per site, with uncertainties in the model's energy scale. It includes the saved production block count, Energy sample count, and uncertainty method.

**Particle-number histogram report**:
A readable account of one simulation's species and total particle-number distributions, including raw counts, normalized probabilities, and particle-number moments. It identifies species and bin labels, available uncertainties, production counts, and estimator conventions.

**Occupied particle-number bin**:
A particle-number bin with a nonzero observed configuration count. The observed particle-number range spans the smallest and largest occupied bin labels and may contain unoccupied bins.

**Peak particle number**:
A particle number whose histogram bin has the largest observed configuration count. A tied peak occurs when multiple bins share that largest count; the reported representative is the smallest such particle number, with the tie explicitly identified.

**Windings report**:
A readable account of one simulation's dimensionless winding means and squares, their totals, and available pair correlations, coflow, and counterflow moments, separated by spatial axis and labeled by component indices. It includes uncertainties, available boundary information, saved production block count, Windings sample count, and uncertainty method.

**Winding component**:
An indexed contribution Wₐ to the winding on a given spatial axis. Its component index identifies its position, not a physical species ID.

**Pair coflow moment**:
The dimensionless expectation ⟨(Wₐ+Wᵦ)²⟩ for two distinct components on the same spatial axis, without a stiffness prefactor or division by the number of components.

**Pair counterflow moment**:
The dimensionless expectation ⟨(Wₐ−Wᵦ)²⟩ for two distinct components on the same spatial axis, without a stiffness prefactor or division by the number of components.

**Raw G(k,τ) aggregate**:
The sample-count-weighted mean of the sector-reweighted block estimates of the momentum–imaginary-time single-particle Green function. Z-visit division and physical normalization have not been applied.
_Avoid_: Physically normalized Green function

**Counterflow equal-time hit histogram**:
Unweighted counts of colocated counterflow endpoints by signed displacement within a one-sided finite imaginary-time window, directed by the first species in an ordered pair. These counts retain sampling-sector and window dependence and are not a physically normalized counterflow correlator.

**Single-particle Green-function report**:
A readable inventory of one simulation's retained single-particle Green-function measurements, explaining their shapes, axes, coordinates, and raw or normalized interpretation. It includes the saved block/sample counts, Z-visit totals, and uncertainty methods for the retained estimates.

**Equal-time Green-function error status**:
A description of the availability and validity of paired jackknife error estimates across the evaluated blocking levels. An "ok" status means those levels were valid; it does not establish convergence or quantify finite-window bias.
