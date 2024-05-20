# `simulation/`

## `pipeline.R`
Simulation pipeline. Uses paleobuddy to simulate species evolution following BiSSE and fossil sampling evolution following a Poisson process, then extracts both an extant (i.e. ultrametric) and FBD (i.e. with fossil tips and sampled ancestors) trees. Parameter variation summary present in a comment close to the beginning of the file. Saves simulations in a `replicates` directory

## `replicates/`
Organized first by fossil sampling parameters, then BiSSE parameters. Seeds used can be found in `seeds.RData`, and raw simulation files (from `paleobuddy::bd.sim.traits`) in `sim_list.RData`. Complete trees can be found in `trees/`, with `trees/ultrametric` containing the corresponding extant trees, and `trees/fbd` the trees with fossil tips and sampled ancestors, using the fossil records contained in `fossils/`. `traits/ultrametric` and `traits/fbd` contain the trait information for the corresponding trees in `trees/`. Lastly, `sims/` simply contains files summarizing each birth-death simulation, equivalent to calling `print(sim)`, where `sim` is an output of `bd.sim.traits`.
