# `analysis/`

## `refs/`
`.Rev` files containing the set up of variables that specify which revbayes script and data to use, and the true parameter combination of the generated data (to organize the output). `modelComb` refers to the model to use, with 1 being FBD, 2 being BiSSE, and 3 being both. `parComb` refers to the parameter combination used to generate the data, with 1 being base (no trait dependence), 2 being high lambda1, 3 being high mu0, and 4 being low q10. `psiComb` refers to the true value of fossil sampling rate used to generate the data, 1 being low, 2 being medium, and 3 high. See `simulation/pipeline.R` for details on the exact values used. Finally, `traitComb` refers to the trait to be analyzed, 1 being the trait that truly affected diversification, and 2-4 being neutral traits with increasing symmetrical transition rates.

## `array.sh`
Slurm script calling an array of 76 jobs (referring to the 76 parameter combinations, see `refs/`), each calling 100 slurm jobs to run the RevBayes scripts. Note that these are set up to run in the Iowa State HPC clusters, and other clusters might have different default options and requirements.

## `bisse.sh`
Slurm script used to run each RevBayes analysis. First it creates an auxiliary `.Rev` file containing the rep number (1-100, supplied by `array.sh` when it calls this script). Then it sources that auxiliary script, the script in `refs/` that specifies the variables regarding parameter combinations, and `master.Rev`.

## `master.Rev`
Master RevBayes script that parses the parameter combination variables supplied through `bisse.sh`, then chooses to call either `fbd.Rev`, `bisse.Rev`, or `both.Rev`. It also loads the tensorphylo plugin, so note that to reproduce this analysis one would need to replace the file path in that line with the path to the tensorphylo library in their machine.

## `fbd.Rev`
RevBayes script to run analyses using only FBD, without any state-dependence. Note that these analyses were only run for comparison purposes if BiSSE or BiSSE+FBD analyses returned unexpected results, and therefore those results are not discussed in the manuscript. It could be interesting for a future study analyzing how FBD reacts to datasets generated with trait-dependent diversification.

## `bisse.Rev` 
RevBayes script to run analyses using only BiSSE, without fossils. When this is called, the script loads the extant trees as data (i.e. the trees in each `trees/ultrametric/` directory).

## `both.Rev`
RevBayes script to run analyses using BiSSE+FBD. When this is called, FBD trees are used (i.e. the trees in each `trees/fbd/` directory).

## `plotting.R`
R script to read the output of the RevBayes analyses and plot the figures for the paper. Note that we omit the output files from the repo due to space considerations, but if analyses are run as-is the `output` directory should be organized in an analogous way to how the `replicates` directory is.
