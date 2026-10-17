# Data dictionary

## `data/extracted/kinetic_parameters.csv`

Rate constants, half-lives and concentrations used by `scripts/12`, one row per parameter per source. 67 rows from 15 primary sources.

| Column | Type | Description |
|---|---|---|
| `param_id` | string | Identifier (K001…). `scripts/12` loads parameters by this id, never by position |
| `kind` | enum | `numeric` (a number printed in the source), `qualitative_constraint` (a statement with no usable number, e.g. "undetectable"), `declared_gap` (searched for and not found), `derived` (computed by us from a primary table) |
| `axis` | enum | `miR-29/BACE1/Abeta`, `miR-7/SNCA/alpha-synuclein`, `both` |
| `parameter` / `symbol` | string | What the value is, and the symbol the model uses |
| `value_as_written` | string | The value exactly as printed in the source. For `numeric` and `qualitative_constraint` it must appear inside `verbatim_quote` |
| `value_si` | float | The same value converted to the model's units. Blank for qualitative constraints and gaps |
| `unit_si` | enum | `1/hour`, `M`, `M^-1 s^-1`, `molecules per cell`, `uM`, `um^3`, `fold`, `fractional change`, `fractional increase`, `fraction of patients`, `dimensionless` |
| `unit_as_written` | string | Unit as printed in the source |
| `population` / `condition` | free text | Who or what was measured, and under which experimental condition. A rate constant without its condition is rejected by `scripts/08` |
| `n` | string | Sample size as reported |
| `method` | free text | Measurement technique |
| `species` | string | Organism or system (human CSF, cell line, rat or mouse neurons, recombinant protein in vitro…) |
| `pmid` / `doi` / `first_author` / `year` / `journal` | | Citation of the primary source |
| `verbatim_quote` | free text | **The sentence the value was read from.** Empty for `derived` and `declared_gap` rows, on purpose |
| `source` | string | Where the text was read (PMC full text, supplementary PDF or table supplied by the author of this repository, Scite full text) |
| `note` | free text | Conversions, caveats, and for `derived` and `declared_gap` rows, how the value was computed or how the search was done |

**Reading rules.**
- A `declared_gap` row carries no value, no quote and no citation. It records that a number was looked for and not found, and where. The ODE code refuses to load it.
- A `derived` row is **not** a number printed by its source. Its derivation is in `note`, and `scripts/08` recomputes it from the archived table in `data/raw/kinetics_2026/`.
- K038 and K039 record values from a model that the source paper itself rejects. They are kept as warnings and not used.
- A value printed in a figure or in a figure's table is recorded with that row or label as its quote, and `note` says which figure it was read from. A value that exists only as a position on a plotted curve is not recorded at all.
- K052 and K054 are the reason the `derived` kind exists. Wilhelm et al. report α- and β-synuclein together (K052) and, in a footnote, the ratio between them (K053); neither is an α-synuclein concentration on its own. K054 combines the two and says so, and `scripts/08` recomputes it. Quoting K052 as α-synuclein would overstate it twofold.

## `data/raw/kinetics_2026/tushev_2018_table_S1.xls`

Supplementary Table S1 of Tushev et al. 2018 (Neuron, DOI 10.1016/j.neuron.2018.03.030): 24,435 3'UTR isoforms with gene symbol, cell-type enrichment, compartment localisation and a `half.life[hours]` column, the processed data the paper's own availability statement points to. `scripts/08` looks up the SNCA, BACE1 and APP rows again from this file and recomputes the pooled BACE1 decay constant. Read it with two cautions, both recorded in the rows that use it: the decay was observed over 16 h, so a quarter of the isoforms carry half-lives beyond the window and the largest reach 18,858 h; and the paper's neuron-enriched summary reproduces only when values above about 25 h are excluded, while the glia summary does not reproduce at all.

## `data/raw/kinetics_2026/wilhelm_2014_table_S1.xlsx`

Additional Data Table S1 of Wilhelm et al. 2014 (Science, DOI 10.1126/science.1252884): quantitative immunoblot measurements of 64 presynaptic proteins, with percentage of total protein, copy number per synapse (mean ± SEM of four preparations), molar concentration and footnotes. The legend is embedded in the sheet as an image. Concentrations were computed by the authors over the synaptic volume minus the mitochondrial volume; `scripts/08` confirms this by recomputing the implied volume from each copy number and concentration pair and comparing it against K050 minus K057, read from Figure 1C of the same paper.

## `data/raw/kinetics_2026/schwanhausser_2011_supplementary_table.xls`

Supplementary table of Schwanhäusser et al. 2011 (Nature, DOI 10.1038/nature10098): genome-wide mRNA and protein half-lives, copy numbers and synthesis rates in mouse NIH3T3 fibroblasts, 5028 rows. Archived unchanged. The file metadata show it was last modified in November 2012; the paper's corrigendum (DOI 10.1038/nature11848) came out in February 2013, and it has not been confirmed whether this is the corrected version. Used only for the genome-wide medians (K040, K041) and to confirm that BACE1 and SNCA are absent (K030).

## `data/raw/structures_2026/`

Deposited coordinates, unchanged: `6N4O.pdb` (human Argonaute2 with miR-122 and target, X-ray 2.9 Å), `4D8C.cif` (BACE1 with a cyclic sulfone inhibitor, X-ray 2.07 Å), `6CU7.cif` (full-length α-synuclein fibril, rod polymorph, cryo-EM), `5OQV.cif` (Aβ(1-42) fibril, cryo-EM 4.0 Å).

## `results/tables/ode_calibrated_results.json`

Output of `scripts/12`: the free parameters, their ranges and the central value of each (measured where a comparable measurement exists); the illustrative aggregation constants and why no source gives them; the comparison between the required mimic dose and the measured BACE1 knockdown; the AD clearance-versus-production experiment; time for a miRNA mimic to wash out at each measured half-life; the α-synuclein pH-gate experiment; human Aβ42 aggregate loads expressed as multiples of M\*; and a summary of the free-parameter scan.

## `results/tables/ode_free_parameter_sensitivity.csv`

One row per point of the free-parameter scan: `free_parameter`, `value`, `abeta_AD_over_control`, `verdict` (`AD_above_control`, `AD_at_or_below_control`, `not_evaluable`). A run that fails numerically is reported as `not_evaluable`, never counted as a reversal. The runs integrate the monomer cascade only (the illustrative aggregation block is left out, because it cannot change the monomer pool and diverged at two corners of the scan); all 18 are evaluable.

## `results/tables/structure_figure_provenance.json`

One entry per structure, read from the coordinate file by `scripts/13`: deposited title, method, resolution and resolution criterion, primary-citation DOI and PMID, the figure's role, and checks computed from the coordinates: for 4D8C the catalytic aspartates found both by sequence motif (DTGS, DSGT) and by distance to the inhibitor; for the fibrils, the chains in each protofilament, the layer spacing and the distance between protofilaments. Residue numbers are those of the deposition.

## `results/tables/mimic_dosing_feasibility.csv`

One row per species and dosing interval: `species`, `param_id` (the kinetic-table row the decay constant came from), `half_life_hours`, `dosing_interval_hours`, `peak_over_average` and `fraction_of_interval_above_half_peak`. No free parameter and no assumed dose enter: the dose cancels out of the peak-to-average ratio.

## `results/tables/mimic_dosing_feasibility.json`

The same calculation with the comparisons that use it: each species' decay constant and its `param_id`, the daily-dosing comparison against the median-stability reference (penalty, equivalent interval, fold stabilisation required), and the stated assumption that the delivered mimic is cleared at the endogenous rate — which is why the output is a required fold stabilisation rather than a verdict on feasibility.

## `results/tables/mimic_target_response.csv`

One row per scenario of `scripts/36`: `axis` (AD or PD), `rho_case` (`scan low`, `scan midpoint`, `scan high`, and for PD `K065 lower bound`, the weakest repression that K065 allows for the pair), `k_repress`, `rho`, the delivered `species` with its `param_id` and `half_life_h`, `interval_h`, `protein_case` and `protein_half_life_h`, the `target` and `target_knockdown`, and whether it is `reachable` at a dose up to 10^4 times the baseline repressor pool. For a reachable scenario: the bolus `dose` and the resulting `peak_excess`, `peak_pool_over_baseline`, `trough_pool_over_baseline` and `mean_excess` (all in units of the baseline pool), the `linear_excess` and `linear_dose` that the steady-state closed form asks for, `efficiency` (`linear_dose` over `dose`, at most 1), `knockdown_reached` (the cycle-averaged protein knockdown, equal to the target), `protein_ripple` (peak-to-trough swing of the protein over its mean), `abeta_ripple` (AD only), `elongation_reduction_m49_8` and `elongation_reduction_m46_0` (PD only: the mean slowdown of fibril elongation on the saturation curve with K009 and K010, at the bouton concentration K054) and `residual` (how closely the periodic state closed). For an unreachable scenario the dose columns are empty and `knockdown_reached` is what the largest dose reaches. The dose does not depend on `protein_case`: over a periodic cycle the mean protein equals the mean mRNA.

## `results/tables/manuscript2_numbers.json`

Written by `scripts/37_manuscript2_numbers.py`: every number that the Article 2 manuscript quotes, as a flat key-to-string map already rounded and formatted for the text, with the source file of each value under `_sources`. The manuscript itself is kept outside this repository and refers to these keys; a build that finds a key missing stops. It selects, rounds and counts from the outputs of `scripts/12`, `17`, `36`, `39` and `40`, from `scripts/_provenance.py` and from the kinetic parameter table, and computes plainly arithmetic quantities (a ratio, a percentage, a time to 90%) that it labels as such. An assert stands behind each sentence of the manuscript that relies on a pattern in the data.

## `results/tables/mimic_target_rho_sweep.csv`

The same calculation at a dosing interval of 24 h over the repression strength ρ from 0.01 to 10⁴ (31 points, five per decade), for every axis, delivered species and target: `rho`, whether the target is `reachable`, the `dose`, `peak_pool_over_baseline`, `linear_dose`, `efficiency`, `protein_ripple` and `knockdown_reached`. It exists so that a statement about the dose does not have to rest on the declared scan of `k_repress`: weak repression raises the dose roughly as 1/ρ, strong repression approaches a limit, and the smallest peak over ρ is a bound that holds whatever ρ is.

## `results/tables/mimic_target_response.json`

The rest of `scripts/36`: the software versions, the integration settings, the repression strength ρ at the ends and midpoint of the declared scan for each axis and, for PD, at the K065 floor (with the mimic factor the closed form asks for in each case), the summary of the ρ sweep (smallest daily peak over ρ, and the peak at the weakest and strongest repression tried), the nucleation-rate scaling from the measured Aβ42 reaction orders, the six verification results with their errors and tolerances, the assumptions in English and Portuguese, and the counts of scenarios and of unreachable scenarios.

## `results/tables/mimic_design_space.csv` and `mimic_design_space.json`

Written by `scripts/40_mimic_design_space.py`. The CSV has one row per grid point, 8,379 in all (nine series, 49 values of ρ, 19 dosing intervals): `axis` (AD or PD), `species` with its `param_id` and `half_life_h`, the `target` and `target_knockdown`, `rho`, `interval_h`, whether the target is `reachable` at a dose up to the cap, the `peak_pool_over_baseline` that reaches it and, where it is not reachable, the `knockdown_reached` at the cap. The JSON holds the grid (`rho_grid`, `interval_grid_h`), the three reference ceilings (`thresholds`: 3, 10 and 30 times the baseline pool), the longest interval within each ceiling for each series at the low, middle and high ends of the declared scan of ρ and, for PD, at the K065 floor (`interval_within_threshold`, with a flag: `interp`, `beyond` when the longest interval of the grid is still within the ceiling, `none` when even the shortest is not), the count of `unreached_points` and the `verification` block (agreement with `scripts/36` at 24 h, monotonicity in the interval, and the closed-form check of the shortest unreached interval).

## `results/tables/structure_renders_article2.json`

Written by `scripts/39_structure_renders_article2.py`: for each PyMOL render drawn for the figures of the mechanistic study (under `results/figures/structures/a2_*.png`), the PDB entry, the pixel size of the crop and what was selected (chains, seed nucleotides, catalytic aspartates). For the side views of the two fibrils, `layer_transformation` lists, per protofilament, the chain that was moved, the chain it was fitted to, the number of C-alpha pairs, the rms of the fit, the rise per layer (angstrom), the rotation per layer (degrees) and the chains along the axis. The translucent layer past the last chain is constructed from that transformation; it is illustration, not a result.

## `results/tables/parameter_provenance_classes.csv`

Written by `scripts/_provenance.py` (through `scripts/37` and `scripts/38`): one row per row of the kinetic parameter table, with `param_id`, `kind`, `axis` (AD, PD or both), `source` (first author and year, empty for the declared gap), `system` (where it was measured: `people`, `human_pm`, `animal_brain`, `neurons`, `cells`, `cell_free` or `none`) and `block` (what it feeds: `mirna`, `mrna`, `repression`, `protein`, `flux` or `aggregation`). The classes are declared in the module and checked against the species, population, condition and method of each row; `scripts/08` requires the saved table to equal the classification.

## `results/tables/graphical_abstract_values.json`

Every number drawn on `results/figures/graphical_abstract.*.png`, loaded from `ode_calibrated_results.json`, `mimic_dosing_feasibility.json` and `kinetic_parameters.csv` at draw time rather than typed as a separate string, plus the `param_id`s each figure grounds a claim in. Nothing here is a new computation; it is a record of which already-verified number went where.

## `results/tables/mechanism_animations.json`

The numbers each animation in `results/figures/anim_*.gif` draws: the final AD/control Aβ42 ratio, the two species' washout times, and the two species' peak-to-average ratios and the penalty between them. Every field is recomputed by `scripts/08` from `ode_calibrated_results.json` and `mimic_dosing_feasibility.json` and required to match, so the manifest cannot drift from the tables its own numbers came from.

## `results/figures/anim_*.gif`

Three animated GIFs, each in both languages: `anim_ad_monomer` (Aβ42 accumulating to steady state, no free parameter), `anim_mimic_washout` (a single mimic bolus decaying at its measured half-life), `anim_dosing_sawtooth` (the repeated-dosing peak-to-average penalty as a moving cycle). Each is produced by `scripts/19_mechanism_animations.py`, which imports its equations from `scripts/12` and `scripts/17` rather than re-deriving them. The α-synuclein pH-gate trajectory is not among them, because its acidic-pH rate constant is illustrative rather than measured (K042).

## `results/figures/`

Every figure exists twice, `<name>.en.png` and `<name>.pt-BR.png`. Both are drawn by the same code from the same arrays in the same run; only the label text differs.

## `results/figures/structure_story_panel.*.png`

The four `scripts/13` renders composited into one figure, with the accompanying facts (PDB id, method, resolution, catalytic dyad residues, protofilament count and spacing) pulled from `structure_figure_provenance.json` rather than retyped. No new render, no new claim - only the layout and the grouping are new.

## Conventions

- Proportions are stored as proportions (0.82), not percentages (82%).
- Empty cells mean **not stated in the source**, never zero and never imputed.
- Every enum value that appears in the data is listed above; unlisted values indicate a data error.

## `results/tables/mimic_series_ranking.csv` and `mimic_series_ranking.json`

Written by `scripts/46`. One row per series, anchor and threshold: `axis`, `target`, `target_knockdown`, `species`, `param_id`, `half_life_h`, `anchor` (scan low, midpoint, high, or K065 floor), `rho`, `threshold_x` (3, 10 or 30), `interval_h` (longest dosing interval within the threshold; 0 when none qualifies, 168 when every interval on the grid does), `flag` (`interp`, `none` or `beyond`), `reference` (True for the median-miRNA reference) and `rank_among_studied` (competition ranking among measured species; empty for the reference). The JSON states the criterion in both languages, the featured series of each axis, whether first place changes, and the ties.

## `data/raw/sequences_2026/` and `data/extracted/utr3_sequences.fa`

Inputs of `scripts/49`, with their provenance in `provenance.json` (URL, `Last-Modified` or ETag when the server gives one, size, SHA-256 of the whole file and of each extract). `mirbase_mature_selected.fa` holds the five mature sequences used (hsa-miR-29a-3p, -29b-3p, -29c-3p, hsa-miR-7-5p, hsa-miR-153-3p) as they are written in a copy of miRBase's `mature.fa` kept by nf-core/test-datasets; the miRBase release is not stated in that file, and the sequence of hsa-miR-7-5p there has 23 nt (a 24-nt form with one more 3′ U is also in circulation; the seed and every pair counted here lie in the first 23 nt). `refGene_selected.tsv` holds the two RefSeq rows (NM_012104, NM_000345) of the UCSC refGene table of hg38. `utr3_sequences.fa` holds the 3′UTR of each transcript (exons after the stop codon, read from the hg38 reference by byte range, 5′ to 3′ of the mRNA, written with U); the header gives the genomic intervals.

## `results/tables/mirna_site_pairing.csv` and `.json`

Written by `scripts/49`. One row per seed-match site of a microRNA in a 3′UTR: `gene`, `transcript`, `utr3_length_nt`, `mirna`, `mirna_accession` (MIMAT), `mirna_length_nt`, `site_type` (8mer, 7mer-m8 or 7mer-A1), `utr_start` and `utr_end` (1-based from the first nucleotide after the stop codon), `a1_position` (the 3′UTR nucleotide that faces nucleotide 1 of the microRNA), `seed_pairs_2_8`, `supplementary_pairs_13_16`, `pairs_total` and `wobble_pairs` (ungapped, Watson–Crick plus G–U), `target_5to3` and `pairing_marks` (`|` Watson–Crick, `:` G–U). It is a count of complementarity, not a binding energy and not evidence that a site is functional.
