This PR mainly provides an example (`docs/src/examples/do-events.jl`) of Critical Slowing Down analysis applied on the $\delta^{18}O$ record from NGRIP. This time series displays DO-events and therefore many abrupt transitions in real-world, paleoclimatic data. Furthtermore, the same analysis is performed on simulation data from CLIMBER-X (less noisy and even time sampling), an Earth Model of Intermediate Complexity able to reproduce DO-like events.

Two things are mildly ugly in this example and are subject to future changes:
- The preprocessing is quite long. Maybe upload preprocessed time series to cut down code of example? On the other hand, maybe nice for people to see a typical preprocessing.
- The transition times of the DO are eye-balled + hand-coded. This should be replaced by an automated procedure like those provided by `ChangePoints.jl`, which however did not give good results so far on the NGRIP data.

This PR additionally comes with minor suggestions:
- Created `docs/src/examples` folder for literate examples.
- Renamed `examples.jl` to `permutation_entropy.jl` and moved to `docs/src/examples`.
- Replaced `docs/src/tutorial.md` by `docs/src/examples/tutorial.jl` such that the tutorial is built by Literate.jl.
- `significant_transitions` was renamed to `infer!` since it mutates the passed `SurrogatesSignificance` (more generally we want to dispatch this function on `<:TransitionSignificance`?). Because of this, the order of the function inputs was changed to comply with the style guidelines `infer!(signif::SurrogatesSignificance, res::WindowedIndicatorResults)`.
- `precompute_ridgematrix` was renamed to `ridgematrix` since it simply computes the matrix.
- `estimate_transitions` is renamed to `transition_metrics` as it describes more accurately what the function does, namely computing the indicators and their change metric (I suggest we generally refer to them as transition metrics). So far, no transition is estimated since it is done in the later step of significance analysis.
- `_cha` is a non-obvious prefix to me



GitHub issue:

Asyncronous significant p-values (maybe bs, just look at your results)


Hi George,

The example with the NGRIP data for the DO-events is finished and I've made a PR for it. Going through the docs page should be quite self-explanatory. I've also made minor changes, which are all listed in the PR. I hope that going through all of that will be smooth for you :)