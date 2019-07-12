# DuMuX-preCICE coupling

## Milestones

### Free flow and poious media flow coupling



We plan to do the following

| Solver | writes | reads |
| ------ | ------ | ----- |
| Darcy | Pressure | Velocity |
| FreeFlow | Velocity | Pressure |

The data is then incorporated into the BJS boundary conditions. Hopefully in a smart way.

- Let's try a `shared_ptr` concept instead of the singleton approach!
    - Does not prevent user to create more than one instance what is somewhat bad!
    - **Does not really work.** So forget about that now. Cannot pass it properly to our hacky helper functions.

- We need to set the initial pressure in the Darcy problem back to 0.

#### Questions

Flow is reversed after 1.5 iterations!

**Initial velocity**
![Initial velocity](./figs/initial_data_pm-ff-coupling.png)
**Velocity after 1 ff step and 1 pm step (=1 coupling step)**
![Velocity after 1 coupling step](./figs/1step_pm-ff-coupling.png)
**Velocity after 2 ff steps and 1 pm step (=1.5 coupling step)**
![Velocity after 2 ff steps and 1 pm step](./figs/1andhalf_step_pm-ff-coupling.png)

**FF solver**
```c++
# 1st iteration
velocities to be sent to ff
v[0]=-6.61414e-08
v[1]=-6.19115e-08
v[2]=-5.6572e-08
v[3]=-5.04326e-08
v[4]=-4.36579e-08
v[5]=-3.63742e-08
v[6]=-2.86902e-08
v[7]=-2.07036e-08
v[8]=-1.25056e-08
v[9]=-4.18234e-09
v[10]=4.18234e-09
v[11]=1.25056e-08
v[12]=2.07036e-08
v[13]=2.86902e-08
v[14]=3.63742e-08
v[15]=4.36579e-08
v[16]=5.04326e-08
v[17]=5.6572e-08
v[18]=6.19115e-08
v[19]=6.61414e-08

# 2nd iteration
velocities to be sent to ff
v[0]=-2.41576e-05
v[1]=-2.28153e-05
v[2]=-2.1064e-05
v[3]=-1.897e-05
v[4]=-1.65748e-05
v[5]=-1.39209e-05
v[6]=-1.10525e-05
v[7]=-8.01602e-06
v[8]=-4.8588e-06
v[9]=-1.62892e-06
v[10]=1.62521e-06
v[11]=4.85508e-06
v[12]=8.01227e-06
v[13]=1.10487e-05
v[14]=1.3917e-05
v[15]=1.65709e-05
v[16]=1.8966e-05
v[17]=2.106e-05
v[18]=2.28111e-05
v[19]=2.41533e-05
```

**PM solver**
```c++
# initial data
p[0]=5e-10
p[1]=5e-10
p[2]=5e-10
p[3]=5e-10
p[4]=5e-10
p[5]=5e-10
p[6]=5e-10
p[7]=5e-10
p[8]=5e-10
p[9]=5e-10
p[10]=5e-10
p[11]=5e-10
p[12]=5e-10
p[13]=5e-10
p[14]=5e-10
p[15]=5e-10
p[16]=5e-10
p[17]=5e-10
p[18]=5e-10
p[19]=5e-10

# 1st iteration
p[0]=-1.6348e-06
p[1]=-1.52962e-06
p[2]=-1.39717e-06
p[3]=-1.24513e-06
p[4]=-1.07754e-06
p[5]=-8.97521e-07
p[6]=-7.07721e-07
p[7]=-5.10538e-07
p[8]=-3.08192e-07
p[9]=-1.02794e-07
p[10]=1.03616e-07
p[11]=3.09015e-07
p[12]=5.1136e-07
p[13]=7.08544e-07
p[14]=8.98344e-07
p[15]=1.07837e-06
p[16]=1.24595e-06
p[17]=1.39799e-06
p[18]=1.53044e-06
p[19]=1.63562e-06
```


In FreeFlow `ffproblem.hh`

- Should alphaBJ and permeability be transfered? I thought they would be constant?

In Darcy

- Do we receive mass flux or velocity?
- How to reconstruct certain things?!
- What happens when `values.setBJS` is called? Is it callable with preCICE running?

**NEW**

- I even get negative pressures. Should that be expected? I would guess pressures would be in between 0 and 1e-9 as prescribed on the domain boundaries. Although I understand that is not total pressure so negative values are not forbidden.
- What are the pressures in the Darcy domain normally? Is there a jump between free flow and Darcy flow?
- Change direction of data transfer to change the Dirichlet/Neumann coupling? The idea is to give the Darcy-solver more freedom when solving its problem. 

| Solver | writes | reads |
| ------ | ------ | ----- |
| FreeFlow | Pressure | Velocity |
| Darcy | Velocity | Pressure |

- Monolithic solver is broken?!



### Conjugate heat transfer problem 

See [test case description](dumux_test_case.md#conjugate-heat-transfer) 


We do the following now:
| Solver | writes | reads |
| ------------ | ------ | ----- |
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
