# DuMuX-preCICE coupling

## Milestones

### Conjugate heat transfer problem 

See [test case description](dumux_test_case.md#conjugate-heat-transfer) 


We do the following now:

| Solver | writes | reads |
| ------ | ------ | ----- |
| SolidEnergy | Temperature | Heat flux |
| FreeFlow | Heat flux | Temperature |

**Attention**
This is the opposite way of writing things than in the OpenFOAM example! If we want to couple to OpenFOAM we have to change that.

#### Problems

- How do we control the adaptive time stepping? It increases above the suggested time step size which makes it impossible for preCICE to catch the solvers.
    - I have adapted the adaptive solver. I hope it works properly now!

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
