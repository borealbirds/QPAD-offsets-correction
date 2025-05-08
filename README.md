# QPAD offset correction

## Background
In late March 2025, we discovered and fixed a bug in the code for calculating QPAD offsets that dated back approximately ten years.

The bug was within the code used to adjust time zones and therefore affects QPAD offsets used for:
1) [species with time since sunrise in the top model](https://github.com/borealbirds/QPAD-offsets-correction/blob/main/qpad_tssr_species.csv) and,
2) in areas outside the mountain time zone.

Since QPAD offsets only adjust the intercept of model estimates, this bug will only affect model outcomes if:
1) models were built for density or population estimates per se, or
2) models were compared or integrated across time zones.
Relative patterns of density (e.g., habitat coefficients) within time zones should be unaffected.

The bug per se has been fixed in all available code, specifically the (`qpad-offsets repository`)[https://github.com/borealbirds/qpad-offsets] and the (`wildrtrax`)[https://github.com/ABbiodiversity/wildRtrax] R package. For previous analyses such as [V4 of the BAM landbird models](https://borealbirds.github.io/), we have developed correction factors for post-hoc correction. 

## Repository components
This repository is intended to document the bug and support users for post-hoc correction. It contains the [correction factors](https://github.com/borealbirds/QPAD-offsets-correction/blob/main/offset-correctionfactors-2025-04-04.csv) and three markdown documents:
1) (An explanation of the bug and how it was fixed)[https://github.com/borealbirds/QPAD-offsets-correction/blob/main/01_qpadbug_explanation.pdf]
2) (An explanation of how the correction factors were derived)[https://github.com/borealbirds/QPAD-offsets-correction/blob/main/02_correctionfactor_calculation.pdf]
3) An example of how the correction factors can be applied

## Support
Please reach out at bamp@ualberta.ca with questions or for support.
