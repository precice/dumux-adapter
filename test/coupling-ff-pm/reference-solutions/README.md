# Short overview over reference solutions

The reference solutions are based on the data used for publication at FVCA IX

## Problem setup

- 40x40 mesh elements
- Stokes flow
- Pressure difference 1e-9
- Permeability 1e-6
- Porosity 0.4
- AlphaBeaversJoseph 1.0

There are two ways to exchange data:

1. Has no extra identifier
  - Pressure: Porous medium -> Free flow
  - Normal velocity/Mass flux: Free flow -> Porous medium
  - **Note** This does not work at the moment and has not been used in the paper.
2. Partially referred to as "reverse" in the current repository
  - Pressure: Free flow -> Porous medium
  - Normal velocity/Mass flux: Porous medium -> Free flow

## Iterative coupling

- Coupling convergence limit: 1e-8

### Serial-implicit coupling

#### Darcy solver first

#### Stokes solver first

## Monolithic coupling

### Erroneous reconstruction of the velocity (published version)

### Correct reconstruction of the velocity (unpublished version)

The error in  the velocity reconstruction has been fixed. For the sake of completeness and to emphasize the correctness of our coupling approaches we add these results here.


## Remarks
