Run from the project root:

```r
shiny::runApp('app_pop-model')
```

The app sources the shared helpers in `R/` and uses
`run_population_model()` to call `evalpheno::phenoflex_population()`.
It needs the current evalpheno package with its population model and
`calculate_chill_dynamic()` / `calculate_heat_gdh_unscaled()` functions,
plus shiny, bslib, patchwork, dplyr, tidyr, ggplot2, lubridate, chillR,
readxl, purrr, scales, ggrepel, ggh4x, yaml and sn (for skewed distributions).
No local C++ compilation is required.

Parameters load from the same YAML files in `parameters/` as the analyses.
The default is `topaz_population_normal.yaml`; select another file or click
**Reload parameters from YAML** to discard edits and read the selected file
again. Loading updates model parameters, standard deviations, distribution
types, skewness, population size, seed and the budbreak heat multiplier.
Edits in the app do not overwrite the YAML files.

Kinetic YAML parameters are used directly and show E0/E1/A0/A1 controls.
Characteristic YAML files show theta_star/theta_c/tau/pi_c controls instead;
temperatures are displayed in Celsius and converted back to Kelvin before
evalpheno converts them to kinetic parameters. Population size and random
seed are editable; a blank seed requests fresh random draws.

Forcing observations use the current files in `data/`. Orchard predictions
use `Ravensburg_hourly_temp_fixed.csv` and Topaz bloom observations.
The heat multiplier applies during forcing; for orchard budbreak it scales
the mean heat requirement, retaining the specified standard deviation,
following the manuscript analysis. Flowering uses the full heat requirement.
Year groups can be combined with individual years.

Integration checks:

```r
source('tests/test_app_population_model.R')
```
