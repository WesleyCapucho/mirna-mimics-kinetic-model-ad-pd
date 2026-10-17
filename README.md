# Kinetic model of microRNA mimics in Alzheimer's and Parkinson's disease: data, code and results

🇧🇷 Versão em português: [README.pt-BR.md](README.pt-BR.md)

This repository holds the data, scripts, tables and figures of a kinetic study of microRNA (miRNA) mimics for two axes: miR-29, BACE1 and Aβ42 in Alzheimer's disease (AD), and miR-7 and α-synuclein in Parkinson's disease (PD). It holds no manuscript text. Every number a manuscript quotes is read from a result table by a script, and `scripts/08_verify_consistency.py` recomputes the derived numbers from their sources.

## What it contains

- **A source-traced parameter table** (`data/extracted/kinetic_parameters.csv`): rate constants, half-lives and concentrations, each tied to the sentence or supplementary table it was read from, with the system in which it was measured.
- **One ODE model per disease**, calibrated only with measured values, and a scan of the one unmeasured quantity, the strength of repression, over six orders of magnitude.
- **The cost of repeated dosing**: the peak a mimic must reach to hold a target average, and what that does to the target mRNA and protein.
- **A design-space map** of the dose needed over repression strength and dosing interval, for nine series.
- **Structure renders** of four deposited structures (Argonaute 2 with a guide and a target, BACE1, Aβ42 and α-synuclein fibrils), for illustration only.

## Relation to the companion repository

The companion repository, [`-Epigen-tica-e-micrornas-em-doen-as-neurodegenerativas`](https://github.com/WesleyCapucho/-Epigen-tica-e-micrornas-em-doen-as-neurodegenerativas) (Zenodo, [10.5281/zenodo.23438256](https://doi.org/10.5281/zenodo.23438256)), holds a systematic review and meta-analysis of circulating miRNAs as diagnostic biomarkers in AD and PD. The two studies share no data and no results. Until this repository was created, the mechanistic scripts lived in the companion repository, so its archived versions up to 1.2.1 still contain them; this repository is their home from now on. Script numbers continue the numbering of the companion repository, so `scripts/12` here is `scripts/12` there.

## Layout

```
├── data/
│   ├── extracted/kinetic_parameters.csv   # The parameter table, one row per parameter per source
│   └── raw/
│       ├── kinetics_2026/                 # Supplementary tables of published papers (see LICENSE-DATA.md)
│       └── structures_2026/               # Deposited structures: 6N4O, 4D8C, 5OQV, 6CU7
├── scripts/
│   ├── 12_ode_models_calibrated.py        # The two models, calibrated with measured values; free-parameter scans
│   ├── 13_structure_figures.py            # Reads the four structures and records their provenance
│   ├── 17_mimic_dosing_feasibility.py     # Peak-to-average cost of repeated dosing by mimic half-life
│   ├── 36_mimic_target_response.py        # Dosing, mRNA and protein together; repression strength; nucleation scaling
│   ├── 40_mimic_design_space.py           # Dose needed over repression strength and dosing interval, nine series
│   ├── 46_mimic_series_ranking.py         # Which series holds a target down for the longest dosing interval; every threshold and anchor ranked
│   ├── 47_mimic_series_figure.py          # Figure of the series ranking: renders of the targets above the interval-over-rho plots (main-text figure)
│   ├── 48_mechanism_algorithm_scheme.py   # Scheme of the molecular mechanism of the two diseases and of the six steps of the algorithm, with renders and small plots drawn from the result tables
│   ├── 39_structure_renders_article2.py   # PyMOL renders of the four structures (illustration only)
│   ├── 37_manuscript2_numbers.py          # Every number a manuscript quotes, read from the tables (numbers only)
│   ├── 38_mechanistic_figures.py          # The seven main-text figures, fitted to the journal's artwork limits
│   ├── 18_graphical_abstract.py           # Schematic figure of both axes, computed numbers
│   ├── 19_mechanism_animations.py         # Animated GIFs of already-simulated mechanisms
│   ├── 20_structure_story_panel.py        # Composites the structure renders into one figure
│   ├── 45_check_figure_artwork.py         # Lettering contrast, line width, page size and image resolution on figure PDFs
│   ├── 08_verify_consistency.py           # Recomputes every derived number from its source; must pass before committing
│   └── _bilingual.py, _journal_fit.py, _provenance.py, _viz_mech.py, _viz_style.py   # Shared helpers
├── results/
│   ├── tables/                            # Model results, scans, design space, provenance classes, manuscript numbers
│   └── figures/                           # Every figure in English and pt-BR (PNG, PDF, EPS), animations, structure renders
├── docs/                                  # METHODS and DATA_DICTIONARY, in English (docs/en) and Portuguese (docs/pt-BR)
├── CHANGELOG.md, CHANGELOG.pt-BR.md       # What changed between versions
├── CITATION.cff, LICENSE, LICENSE-DATA.md
└── requirements.txt
```

## Reproducing the analysis

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

python scripts/12_ode_models_calibrated.py
python scripts/17_mimic_dosing_feasibility.py
python scripts/36_mimic_target_response.py        # needs 12 and 17; about 3 minutes
python scripts/40_mimic_design_space.py            # needs 12, 17 and 36; about 8 minutes
pip install pymol-open-source                     # only for scripts/13 and 39; it installs for the Python that pip belongs to, so run those two with that Python
python scripts/46_mimic_series_ranking.py          # needs 40
python scripts/13_structure_figures.py
python scripts/39_structure_renders_article2.py   # needs the provenance table of 13; several minutes of headless ray tracing
python scripts/37_manuscript2_numbers.py           # needs 12, 17, 36, 39 and 40
python scripts/38_mechanistic_figures.py           # needs 12, 13, 17, 36, 39 and 40
python scripts/47_mimic_series_figure.py           # needs 39, 40 and 46
python scripts/48_mechanism_algorithm_scheme.py     # needs 12, 13, 17, 36, 39, 40 and 46
python scripts/19_mechanism_animations.py
python scripts/18_graphical_abstract.py            # needs 12, 17 and 13
python scripts/20_structure_story_panel.py
python scripts/08_verify_consistency.py            # must pass before committing
```

`scripts/08` enforces the consistency of parameters, model results, the design-space grid, the renders and the figures. 991 checks currently pass.

## Limits

The model restates what has been measured and scans what has not. It is not a pharmacokinetic prediction for any mimic, and the renders illustrate molecules at the steps of the two axes; they test nothing. See `docs/en/METHODS.md`.

## Citation

See `CITATION.cff`. The DOI of the first archived version is added when it is released on Zenodo.

## Licences

Scripts: MIT (`LICENSE`). Data, tables, figures and documentation: CC BY 4.0 (`LICENSE-DATA.md`), except third-party material, which keeps its own terms; the notice lists what is excluded.
