---
title: The DuMuX adapter
permalink: adapter-dumux.html
url: /adapter-dumux.html
keywords: DuMuX, DUNE, C++
summary: A DuMuX-specific DUNE module for coupling to other codes with preCICE.
---

This is a [DuMuX](https://dumux.org/) adapter (a DUNE module specific to DuMuX) to couple to other codes using [preCICE](https://precice.org/). You can find the source code of this adapter [on GitHub](https://github.com/precice/dumux-adapter).

[Get](adapter-dumux-get.html) and [learn how to use](adapter-dumux-use.html) the adapter.

## Supported features

- Surface and volume coupling as demonstrated in the two tutorial cases [free-flow-over-porous-media](https://precice.org/tutorials-free-flow-over-porous-media-2d.html) and [two-scale-heat-conduction](https://precice.org/tutorials-two-scale-heat-conduction.html).
- Any data field can be exchanged. The adapter does not have any awareness of the physical definition of the exchanged fields.
- In the configuration, one mesh can be set for each interface, and multiple data fields can be read or written on each interface. Multiple interfaces can be configured.
- Implicit coupling and subcycling are supported, while the time step size in the DuMux solver is to be managed by the user.
- Mesh connectivity is not supported.

## Release strategy

Any change that enforces the use of a new DuMuX version on the user side, or any change to the adapter API, is deemed a breaking change and warrants a major release. We follow the semantic versioning scheme: `v{major version}.{minor version}.{patch}`.

## Cite

There is no code-specific publication related to the DuMuX adapter available yet.

### Publications using dumux-precice

You may find more examples and related theory in these publications.

- Jaust A., Weishaupt K., Mehl M., Flemisch B. (2020) Partitioned Coupling Schemes for Free-Flow and Porous-Media Applications with Sharp Interfaces. In: Klöfkorn R., Keilegavlen E., Radu F., Fuhrmann J. (eds) Finite Volumes for Complex Applications IX - Methods, Theoretical Aspects, Examples. FVCA 2020. Springer Proceedings in Mathematics & Statistics, vol 323. Springer, Cham. <https://doi.org/10.1007/978-3-030-43651-3_57>

  - Code can be found at: <https://git.iws.uni-stuttgart.de/dumux-pub/jaust2020a>
