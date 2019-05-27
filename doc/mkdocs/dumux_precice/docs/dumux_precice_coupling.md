# DuMuX-preCICE coupling

## Milestones

### Free flow and poious media flow coupling

We plan to do the following

| Solver | writes | reads |
| ------ | ------ | ----- |
| FreeFlow | Pressure | Velocity |
| Darcy | Velocity | Pressure |

The data is then incorporated into the BJS boundary conditions. Hopefully in a smart way.

- Let's try a `shared_ptr` concept instead of the singleton approach!
    - Does not prevent user to create more than one instance what is somewhat bad!
    - Does not really work. So forget about that now. Cannot pass it properly to our hacky helper functions.

#### Questions

In FreeFlow `ffproblem.hh`

- Should alphaBJ and permeability be transfered? I thought they would be constant?

In Darcy

- Do we receive mass flux or velocity?
- How to reconstruct certain things?!
- What happens when `values.setBJS` is called? Is it callable with preCICE running?


### Conjugate heat transfer problem 

See [test case description](dumux_test_case.md#conjugate-heat-transfer) 


We do the following now:

| Solver | writes | reads |
| ------ | ------ | ----- |
| SolidEnergy | Temperature | Heat flux |
| FreeFlow | Heat flux | Temperature |

**Attention**
This is the opposite way of writing things than in the OpenFOAM example! If we want to couple to OpenFOAM we have to change that.


#### Code related questions

What is the best way to reload a checkpoint?

Saving `sol` and then

1. reloading `sol` via `GridVariables->update(sol)` followed by a `GridVariables->advanceTimeStep()`.
1. reloading `sol` via `GridVariables->init(sol)`;

Codes:


```c++
GridVariables->update(sol);
GridVariables->advanceTimeStep();
```

```c++
GridVariables->init(sol)
```

#### Problems

- [x] We have to be careful with the initialization of the solvers. 
    - **Valid case** Heat solver starts: The heat solver computes a temperature to be send, the heat flux is zero. Thus, the temperature equals the temperature in the cell center.
    - **Invalid case** Flow solver starts: The flow solver computes a heat flux from its temperatures. However, the temperature on the interface (stored in the precice adapter) is zero! Thus, a **very high** heat flux is initialized.
    - **Note**: There were some errors in the preCICE configuration file.
- [x] The implicit coupling creates a solution that develops too fast.
  - **Needs more testing** The checkpointing die only reload parts of the solution. In the end the Newton solver continued from its previous state, but with the initial guess being reset. 
- [x] How do we control the adaptive time stepping? It increases above the suggested time step size which makes it impossible for preCICE to catch the solvers.
    - I have adapted the adaptive solver. I hope it works properly now!
    - Seems to be fine



#### Work 

**Done**:

- Test case with heat transfer is up and running
    - Follows the idea of the preCICE-tutorial test case, but is not exactly the same
        - Compressible flow
        - Slightly different boundary conditions
        - Neumann-Neumann boundary condition
        - Computation of heat flux is a bit different (harmonic mean) `lambda_harmonic = 2 / ( 1 / lambda_solid + 1 / lambda_fluid )`
        - Working with dimensions
- Set up monolithic solver
- Set up two separate solvers
- Set up project with preCICE

**On the way**

- Implementing an adapter/wrapper

**Todo**:

- Write a preCICE-DuMuX-adapter
    - Couple two separate solvers with preCICE **partially**
    - Read/write temperature and heat flux
- Move `PreciceWrapper` into a different directory.
- Check solution of monolithic and partitioned approach

#### Requirements

What functions do we need?

Heat fluxes:

- `readSolidHeatFlux`
- `writeSolidHeatFlux`
- `writeFluidHeatFlux`
- `readFluidHeatFlux`

Are they read/written pointwise? Can we read/write them as vectors/arrays?

Geometry:

- How do we communicate the geometry to preCICE? 
    - **Preferred way**: Hand over a DuMuX data structure to preCICE. This would be easier for DuMuX users. Similar to fracturebenchmark setting?
    - **Bad way**: Collect points in an array in DuMuX and hand it over to the adapter.

Generation of interface data in fracture case:
```c++
using FractureFVGridGeometry = typename FVGridGeometry::Type<1>;
const auto numFracturePoints = fvGridGeometry[onePFractureId].numDofs();
static constexpr int dimWorld = FractureFVGridGeometry::GridView::dimensionworld;
std::vector<double> pointCloud(numFracturePoints*dimWorld);

for (const auto& element : elements(fvGridGeometry[onePFractureId].gridView()))
{
    auto fvGeometry = localView(fvGridGeometry[onePFractureId]);
    fvGeometry.bind(element);

    for (const auto& scv : scvs(fvGeometry))
        for (unsigned int i = 0; i < dimWorld; ++i)
            pointCloud[scv.dofIndex()*dimWorld + i] = scv.dofPosition()[i];
}
```
