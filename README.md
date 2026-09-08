# Toxic effects of *Alexandrium pseudogonyaulax* on the pelagic food web

R code accompanying Möller et al. (2024), *Harmful Algae*. Archived, not maintained.

Bioassays testing the lethal and sublethal effects of the harmful dinoflagellate
*Alexandrium pseudogonyaulax* on four marine trophic levels: microalgae,
microzooplankton, mesozooplankton and fish gill cells.

## Citation

> Möller, K., Tillmann, U., Pöchhacker, M., Varga, E., Krock, B., Porreca, F.,
> Koch, F., Harris, T.M., Meunier, C.L. (2024). Toxic effects of the emerging
> *Alexandrium pseudogonyaulax* (Dinophyceae) on multiple trophic levels of the
> pelagic food web. *Harmful Algae* 138, 102705.
> https://doi.org/10.1016/j.hal.2024.102705

Open access, CC BY 4.0.

## Data

The curated, citable datasets are archived on PANGAEA:

| Dataset | DOI |
|---|---|
| *A. pseudogonyaulax* and *A. tonsa* predator-prey interactions (bundle) | [10.1594/PANGAEA.967577](https://doi.org/10.1594/PANGAEA.967577) |
| — Ingestion rates, three *A. tonsa* life-stages | [10.1594/PANGAEA.967353](https://doi.org/10.1594/PANGAEA.967353) |
| — Intracellular goniodomin A content (LC-MS/MS) | [10.1594/PANGAEA.967418](https://doi.org/10.1594/PANGAEA.967418) |
| — *A. tonsa* egg hatching rate | [10.1594/PANGAEA.967422](https://doi.org/10.1594/PANGAEA.967422) |
| — Temporal depletion of goniodomin A in *A. tonsa* | [10.1594/PANGAEA.967425](https://doi.org/10.1594/PANGAEA.967425) |
| *P. kofoidii* predator-prey interactions with three *Alexandrium* species | [10.1594/PANGAEA.967725](https://doi.org/10.1594/PANGAEA.967725) |
| RTgill-W1 viability and membrane integrity bioassays | [10.1594/PANGAEA.968675](https://doi.org/10.1594/PANGAEA.968675) |
| *R. salina* lysis, *A. pseudogonyaulax* and *A. monilatum* supernatants | [10.1594/PANGAEA.968485](https://doi.org/10.1594/PANGAEA.968485) |
| *R. salina* lysis, purified goniodomin congeners | [10.1594/PANGAEA.968492](https://doi.org/10.1594/PANGAEA.968492) |

**On the two versions of the data.** The files on PANGAEA were restructured for
archiving: long format, one row per observation, explicit column names, added
metadata and taxonomic identifiers. The scripts in this repository were written
against the original working files, which used a wide layout with one column per
treatment. Those working files are included under `data/` so that the analysis
runs exactly as it did for the paper. **PANGAEA is the authoritative version and
the one to cite.** Use `data/` only to reproduce these scripts.

## Contents

| Script | Experiment | Figure / table |
|---|---|---|
| `A_tonsa_predator_prey_interactions_ingestion_rates.R` | Ingestion rates of *A. tonsa* on three *A. pseudogonyaulax* strains, across N4-nauplii, C4-copepodites and adults | Fig. 2a |
| `A_tonsa_predator_prey_interactions_goniodomins.R` | Intracellular goniodomin content of *A. pseudogonyaulax* in the grazing experiments | Fig. 2b |
| `A_tonsa_predator_prey_interactions_hatching_rates.R` | Egg hatching success of *A. tonsa* exposed to cell-free supernatants for 48 h | Fig. 3 |
| `copepodamide_toxin_induction.R` | Toxin induction in *A. pseudogonyaulax* exposed to copepodamides from *C. finmarchicus* and *A. tonsa* | Fig. 4 |
| `P_kofoidii_predator_prey_interactions.R` | Grazing and mortality of *Polykrikos kofoidii* on *A. pseudogonyaulax*, *A. catenella* and *A. limii*, monoalgal and mixed prey | Fig. 1, Table 1 |
| `RTgillW1_assays.R` | RTgill-W1 gill cell assays (CellTiter-Blue metabolic activity, LDH release) and *Rhodomonas salina* lysis bioassays. Dose-response curves and EC50s for purified goniodomins and cell-free supernatants | Fig. 5, Fig. 6, Table 3 |
| `A_tonsa_predator_prey_interactions_extra_experiments.R` | Depletion of goniodomin A in *A. tonsa* over eight days after switching from *A. pseudogonyaulax* to non-toxic *R. salina* | none, see below |

Each script is standalone. There is no shared functions file and no driver
script, so run whichever analysis you need.

**A note on `extra_experiments.R`.** This analysis produces no figure or table in
the paper and is not discussed in the text. It is kept because the underlying
data is part of the archived record (PANGAEA 967425), and code for a deposited
dataset should exist somewhere. Treat it as a supplementary side experiment, not
as part of the published analysis.

## Methods implemented

- Four-parameter log-logistic dose-response fitting and EC50 estimation with
  `drc`, applied per treatment and plate.
- Outlier screening with `Routliers` (median absolute deviation) and a Dixon
  test before analysis.
- One-way ANOVA with Tukey HSD, or Kruskal-Wallis with a Conover-Iman post hoc
  test where normality or homoscedasticity failed. Repeated-measures ANOVA for
  factors measured over time. P-values adjusted after Benjamini-Hochberg.
- Ingestion and clearance rates calculated after Frost (1972) from cell count
  time series, converted to ingested carbon using per-cell POC.
- Figure assembly with `ggplot2`, `ggpubr` and `patchwork`, colourblind-safe
  palettes throughout.

Packages are installed on demand via `pacman::p_load`.

## Running it

Each script opens with a configuration block:

```r
data_dir <- "C:/path/to/your/data"
out_dir  <- "C:/path/to/your/output"
```

Set `data_dir` to the `data/` folder of this repository (`"data"` if you open the
repository as an RStudio project), and `out_dir` to wherever figures and exported
tables should go. Then run any script top to bottom.

```r
source("RTgillW1_assays.R")
```

## Known limitations

The scripts were archived as they were used and have not been refactored.

- `RTgillW1_assays.R` calls `extrafont::loadfonts(device = "win")`, which fails
  outside Windows. Comment out that line on Linux or macOS. Plot fonts fall back
  to the system default.
- Grouping is positional in several scripts. Treatments are reconstructed from
  column order and row order rather than from named factor columns, so the input
  files must keep their original layout. This is why the working files are
  included rather than the PANGAEA versions.
- POC per-cell conversion factors in `ingestion_rates.R` are hardcoded literals
  transcribed from a spreadsheet. They are also available as proper columns in
  the PANGAEA datasets.
- Package loading is scattered through the scripts rather than collected at the
  top, so a partial run may hit a missing package mid-script.
- Some exploratory blocks are commented out rather than removed. Those are not
  part of the published analysis.

## Status

Archived on publication. Kept for transparency. Issues are not monitored.

## License

Code: MIT. Data: CC-BY 4.0, as archived on PANGAEA.
