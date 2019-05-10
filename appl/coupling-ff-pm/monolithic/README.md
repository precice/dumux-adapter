# Exercise Coupling free flow/porous medium flow

The aim of this exercise is to get familiar with setting up coupled free flow/porous medium flow problems.

## Problem set-up

The model domain consists of two non-overlapping subdomains.
Free flow is modeled in the upper subdomain, while the lower subdomain models a flow in a porous medium. Both single-phase flow and two-phase flow will be considered in the porous domain.


### 0. Getting familiar with the code

* Navigate to the directory `exercises/exercise-coupling-ff-pm`

There are three sub folders: `interface` (Exercise 1), `models` (Exercise 2) and `turbulence` (Exercise 3).

The problem-related files for this exercise are:
* Three __main files__ for the three sub-tasks :`ex_interface_coupling_ff-pm.cc`, `ex_models_coupling_ff-pm.cc`, `ex_turbulence_coupling_ff-pm.cc`,
* Three __free flow problem files__: `ex_interface_ffproblem.hh`, `ex_models_ffproblem.hh`, `ex_turbulence__ffproblem.hh`
* Three __porous medium flow problem files__: `ex_interface_pmproblem.hh`, `ex_models_pmproblem.hh`, `ex_turbulence_pmproblem.hh`
* The __input files__: `ex_interface_coupling_ff-pm.input`, `ex_models_coupling_ff-pm.input`, `ex_turbulence_coupling_ff-pm.input`,
* The __spatial parameters files__: `1pspatialparams.hh`, `2pspatialparams.hh`


In the main file, `TypeTags` for both submodels are defined.
The same applies for types such as `GridManager`, `FVGridGeometry`, `Problem`, etc...
Since we use a monolithic coupling scheme, there is only one `Assembler` and one `NewtonSolver`.

The problem files very much look like "regular", uncoupled ones with the exception that they hold a pointer to the `CouplingManager` which allows to evaluate the coupling conditions and to exchange information between the coupled models.
The coupling conditions are realized technically in terms of boundary condition. For instance, in lines 178 and 179
in `ex_interface_ffproblem.hh`, `couplingNeumann` boundary conditions are set, which means that the free flow models evaluates the
mass and momentum fluxes coming from the porous domain and uses these values as boundary conditions at the interface.

Note the certain checks are performed when combining different models, e.g., the fluid system has to be the same for both domains
and the sub-control-volume faces at the interface have to match.


We will use a staggered grid for the free flow and a cell-centered finite volume method for the porous medium.
Keep in mind that the staggered grid implementation distinguishes between face variables (velocity components) and
cell center variables (all other variables), therefore in some cases either the `stokesCellCenterIdx`
or the `stokesFaceIdx` is used respectively, while for the porous medium all variables can be accessed with `darcyIdx`.

__Task__:
Take a closer look at the Stokes/Darcy coupling files before moving to the next part of the exercise:


### 1. Changing the interface

In this part of the exercise, a simple coupled system consisting of a one-phase (1p) free flow and a one-phase flow in a porous medium is set up. Both subproblems have no-flow boundaries at the sides.
Currently, a velocity profile is set on the upper free flow boundary, which leads to a vertical flow into the porous medium:

![](../extradoc/ex_ff-pm-vertical-flow.png)

* We will first change the flow direction such that the free flow is parallel to the porous medium.
* Afterwards, the Beavers-Joseph-Saffman condition will be used as an interface condition for the tangential momentum transfer.
* Last, we change the flat interface between the two domains to a wave-shaped one.

__Task A: Change the flow direction__

Open the file `ex_interface_ffproblem.hh` and navigate to line 148, where the types of boundary condition are set.
Instead of applying a fixed velocity profile at the top of the domain, we want to use fixed pressure boundary conditions
at the left and right side of the free flow domain, while the top represents an impermeable wall.

Set a Dirichlet boundary condition for the pressure at the left and right side of the domain:
``` cpp
if(onLeftBoundary_(globalPos) || onRightBoundary_(globalPos))
    values.setDirichlet(Indices::pressureIdx);

```

Set a Dirichlet boundary condition for the velocities at the top:
``` cpp
if(onUpperBoundary_(globalPos))
{
    values.setDirichlet(Indices::velocityXIdx);
    values.setDirichlet(Indices::velocityYIdx);
}
```

Keep the coupling boundary condition:
``` cpp
if(couplingManager().isCoupledEntity(CouplingManager::stokesIdx, scvf))
{
    values.setCouplingNeumann(Indices::conti0EqIdx);
    values.setCouplingNeumann(Indices::momentumYBalanceIdx);
    values.setDirichlet(Indices::velocityXIdx); // assume no slip on interface
}
```

Having changed the types of boundary conditions, we must now assign the correct values for them.

Set a no-slip, no-flow condition for the velocity at the top:
``` cpp
values[Indices::velocityXIdx] = 0.0;
values[Indices::velocityYIdx] = 0.0;
```
Apply a fixed pressure difference between the inlet and outlet, e.g.:
``` cpp
if(onLeftBoundary_(globalPos))
    values[Indices::pressureIdx] = deltaP_;
if(onRightBoundary_(globalPos))
    values[Indices::pressureIdx] = 0.0;
```

For changing the flow direction, the boundary conditions for the porous medium have to be changed as well.

Use Neumann no-flow boundaries everywhere, keep the coupling conditions.
``` cpp
values.setAllNeumann();

if (couplingManager().isCoupledEntity(CouplingManager::darcyIdx, scvf))
    values.setAllCouplingNeumann();
```

This should make the flow go from left to right.

__Task B: Include slip-condition__

However, we are still missing one important feature:
at the moment, the velocity component tangential to the interface gets a no-slip condition.
In the next step we want to implement the Beavers-Joseph-Saffman slip condition at the interface:

$`\frac{\partial v_x}{\partial y} = \frac{\alpha}{\sqrt K} (v_x - q_{pm})\quad`$ at $`\quad y=0`$

with  $`\quad q_{pm}=0`$.

To include this, just replace the no-slip condition at the interface
``` cpp
values.setDirichlet(Indices::velocityXIdx); // assume no slip on interface
```
with a Beavers-Joseph-Saffman (BJS) boundary condition for the respective momentum balance:
``` cpp
values.setBJS(Indices::momentumXBalanceIdx);
```

at the position where the coupling boundary conditions are set in `ex_interface_ffproblem.hh`.

To check if the simulation behaves as expected, we can compare the velocity profile $`v_x(y)`$ with the analytical solution provided by [Beavers and Joseph (1967)](https://doi.org/10.1017/S0022112067001375).
For doing so, we uncomment line 212 in `ex_interface_coupling_ff-pm.cc`.
```cpp
stokesVtkWriter.addField(stokesProblem->getAnalyticalVelocityX(), "analyticalV_x");
```

After re-compiling and re-running the executable, we should be able to see also
the analytical solution of $`v_x`$ on the free flow domain. Play around with the grid resolution to see how that affects the velocity profile.

__Task C: Cange shape of interface__

Now we want to include a non-flat interface between the two domains. We use `dune-subgrid` to construct
two grids for the two domains from one common host grid. Comment out lines 96-106 in `ex_interface_coupling_ff-pm.cc` and comment lines 114-149 in the same file. This will instantiate a host grid and define two helper lambda functions that are used to choose elements from to host grid for the respective sub grid. In the given case,
the domain is split in two haves, separated by a sinusoidal interface.

```cpp
auto elementSelectorStokes = [&](const auto& element)
{
    double interface = params.amplitude * std::sin(( element.geometry().center()[0] -params.offset) / params.scaling * 2.0 * M_PI) + params.baseline;
    return element.geometry().center()[1] > interface;
};

auto elementSelectorDarcy = [&](const auto& element)
{
    double interface  =  params.amplitude * std::sin(( element.geometry().center()[0] - params.offset) / params.scaling * 2.0 * M_PI) + params.baseline;
    return element.geometry().center()[1] < interface;
};
```

Make sure, that you have uncommented the lines including the gridcreators in both problem files
```cpp
#include <dumux/io/grid/subgridgridcreator.hh>
```

and do the changes in the respective lines for the `Grid` property.

The problem should now compile. However, an error occurs due to the coupling conditions.
So far, we assumed a flat interface, therefore the normal momentum coupling condition

 $`[\sigma \cdot \mathbf{n}]^{FF} = p^{PM}`$

 was always set for a fixed $`\mathbf{n} = (0,1)^T`$. We need to account for the curvature of the interface and thus replace
 ```cpp
values.setCouplingNeumann(Indices::momentumYBalanceIdx);
 ```
 with
 ```cpp
values.setCouplingNeumann(scvf.directionIndex());
 ```
The same if true for the BJS condition, however, here we need to consider the tangential direction:
```cpp
values.setBJS(1 - scvf.directionIndex());
```

The final result should look something like this:

![](../extradoc/ex_ff-pm-wave-interface.png)


### 2. Changing the porous medium model

In this part of the exercise, the coupled system will be extended such that the transport of components
in both parts is included and the presence of a second phase in the porous medium is considered.
This enables the simulation of the drying of a porous medium (however no energy balance is included yet).

This part of the exercise consists of the following steps:
* replacing the 1pnc porous-medium model by a 2pnc porous-medium  model,
* adding the output for fluxes and storage terms (having gnuplot is recommended),
* enable a variable switch for a complete drying.

We start with an example in which the transport of water vapor in the gas phase is considered.
The porous medium is filled with gas, initially the mole fraction of water vapor is $`x^w_g = 0.1`$.
Above the porous medium, a dry gas flows and by diffusion, the porous medium dries out.
(Don't mind the compiler warning, we will deal with it in task B.)


__Task A: Change the model__:

In the first task, the porous-medium model will be changed from a 1p2c system to a 2p2c system.
Although a 2p2c system is plugged in, we still want to simulate the same situation as before, i.e., air with water vapor in both domains.
The following changes have to be made in the porous-medium model (`ex_models_pmproblem.hh`):
* Include the 2pnc model: include the respective headers and inherit from the new model `TwoPNC`
* Exchange the spatial parameters for the 1-phase system by those for a 2-phase system (hint: two occurrences).
* Since two phases are involved now, we do not need to use the OnePAdapter anymore. Change to property of the FluidSystem such that `H2OAir` is used directly.  
  Afterwards, set the `transportCompIdx` to `Indices::switchIdx`.

One big difference between the 1p and 2p model is, that the primary variables for which
the problem is solved, are not fixed.
It is possible to use different formulations, e.g. with the gas pressure and the
liquid saturation as primary variables (p_g-S_l -> `p1s0`) or vice versa.
* Set the property

```
template<class TypeTag>
struct Formulation<TypeTag, TTag::DarcyOnePNC>
{ static constexpr auto value = TwoPFormulation::p1s0; };
```
  in the Properties section in the problem file.

In contrast to the formulation, which stays the same during one simulation, the meaning of
the primary variables may vary during one simulation.
In the case we are investigating, we want to use the gas pressure and the liquid saturation as primary variables.
However, if only the gas phase is present, the liquid saturation is always zero.
In this case, the chosen formulation will set the given value as the mole fraction of water vapor in the gas phase.
* To tell to program which phases are present in which parts of the domain at the beginning of the simulation,
  you have to call `values.setState(MY_PHASE_PRESENCE);` in `initialAtPos()`. Have a look at the `indices.hh`
  in the `2p2c` model to figure out which is the correct value of `MY_PHASE_PRESENCE` for the presence of
  a gas-phase only (hint: the numbering of phase indices begins with 0, the numbering of the phase presence states begins
  with 1. Take a look at your formulation to find out which phase index to use for the gas phase.)

__Task B: Add output__:

In the next step, we want to add some output to the simulation.
First we want to know the water mass in the (porous-medium) system.
Therefore, we evaluate the storage term of the water component.
* Have a look at the function `evaluateWaterMassStorageTerm()` in the porous medium subproblem.
  Then implement a calculation of the total water mass:
  $`\sum_{\alpha \in \textrm{g,l}} \left( \phi S_\alpha X^\text{w}_\alpha \varrho_\alpha V_\textrm{scv} \right)`$.
  Afterwards, adapt the method `init()` such that the variable `initialWaterContent_` is initialized correctly using the `evaluateWaterMassStorageTerm()` method and assign that value also to the variable `lastWaterMass_`.

We also want to investigate the temporal evolution of the water mass.
The following steps need to be done to do so.
Check if all instructions are implemented accordingly:
* Calculate the initial water mass at the beginning of the simulation and add the water mass loss to the output.
  Based on the water mass loss you can derive the evaporation rate. Because despite at the interface, all
  boundaries have Neumann no-flow conditions and no sources are present, the water can only leave the system
  via the porous-medium free-flow interface. The evaporation in [mm/d] can be calculated by:
  $`e = \frac{\textrm{d}\, M^\text{w}}{\textrm{d}\, t} / A_\textrm{interface}`$.

Finally we want to know the distribution of the water mass fluxes across the interface.
* Have a look at the function `evaluateInterfaceFluxes()` in the porous medium problem.
  Use the facilities therein to return the values of ...massCouplingCondition... from the `couplingManager`
  for each coupling scvf. Then the fluxes are visualized with gnuplot, when setting `Problem.PlotFluxes = true`.
  If the simulation is too fast, you can have a look at the flux*.png files after the simulation.
* You can use the parameter `Problem.PlotStorage = true` to see the temporal evolution of the evaporation rate
  and the cumulative water mass loss.

Compile and run the simulation and take a look at the results.

__Task C: Let it dry out__:

In this exercise we want to completely dry an initial water-filled porous medium.
- Change the initial condition such that it is started with a two-phase system.
  Set the initial saturation to $`S_\text{l} = 0.1`$.

If one wants to simulate the complete drying of a porous medium, the standard capillary pressure-saturation
relationships have to be regularized. If no regularization is used, then the capillary pressure would
approach infinity when the liquid saturation goes to zero. This means that water can flow to the
interface from any region of the porous medium and of course this poses numerical problems.
Therefore, the capillary pressure-saturation relationship has to be regularized for low saturations.
The regularization has a strong influence on how long liquid water is present at the interface, see
[Mosthaf (2014)](http://dx.doi.org/10.18419/opus-519).
* Take a look at how the regularization is set in the `2pspatialparams.hh` and see how adapating
  the parameter for `pcLowSw` and `pcHighSw` affects the results.

The porous-medium model can now reach a liquid saturation of zero. As already explained above,
for this case the liquid saturation cannot serve as primary variable anymore. However,
manually adapting the primary variable states and values is inconvenient.
[Class et al. (2002)](http://dx.doi.org/10.1016/S0309-1708(02)00014-3)
describe an algorithm to switch the primary variables, if phases should appear or disappear during a simulation.

Now you are able to simulate a complete drying of the porous medium.


### 4. Use a turbulence model in the free flow domain

Several RANS turbulence models are implemented in DuMuX.
This part of the exercise consists of the following steps:
* replacing the Navier-Stokes model by the zero equation turbulence model,
* switching to a symmetry boundary condition,
* applying a grid refinement towards the interface,
* subsequently refining the grid (convergence study).

We will work with a `1p2cni/2p2cni` coupled problem, where `ni` stands for non-isothermal.
All the prepared files can be found in the subfolder `exercise-coupling-ff-pm/turbulence`.

__Task A__:

The file `ex_turbulence_ffproblem.hh` is your free flow problem file within this exercise.

For using the compositional zero equation turbulence model, the following header files need to be included:
```
#include <dumux/freeflow/compositional/zeroeqncmodel.hh>
#include <dumux/freeflow/rans/zeroeq/problem.hh>
```
The includes for the NavierStokesNC model and the NavierStokesProblem are no longer needed and can be removed.

Make sure your free flow problem inherits from the correct parent type:
* Change the entry in the `StokesZeroEq` definition accordingly (non-isothermal zero equation model)
* Adapt the inheritance of the problem class (hint: two occurrences)

Take a look into the two headers included above to see how the correct TypeTag and the Problem class the inherit from are called.

Here, the turbulent free flow is wall bounded which means that the main reason for the development
of turbulent flow is the presence of walls.
The porous medium at the bottom and also the channel wall are such walls.
Inside the problem you have to tell the turbulence model, where these walls are located
by providing an `isOnWallAtPos()` function:
```c++
bool isOnWallAtPos(const GlobalPosition& globalPos) const
{
    return (onLowerBoundary_(globalPos) || onUpperBoundary_(globalPos));
}
```

In addition, especially for the zero-equation models, any element in the free-flow domain interacts with the walls,
e.g. this defines the wall distance which is needed to calculate the eddy viscosity.
To get all these interactions, you have to call  
 ```cpp
 stokesProblem->updateStaticWallProperties()
 ```
in `ex_turbulence_coupling_ff-pm.cc`.
However, there is also a solution-dependent component of these interactions, e.g. for a correct
damping of the eddy viscosity toward the wall, the velocity gradient at the wall and inside the
cells is needed.
These dynamic interactions are to be updated by calling  
```cpp
stokesProblem->updateDynamicWallProperties(stokesSol)
```
in the time loop (after `// update dynamic wall properties`).

Compile and run your new coupled problem and take a look at the results in Paraview.
In addition to the standard variables and parameters, you can now analyze turbulence model specific quantities
(e.g. the turbulent viscosity $`\nu_\textrm{t}`$ or the turbulent diffusivity $`D_\textrm{t}`$) for the free flow domain.
In paraview you may compare the magnitude of $`D`$ and $`D_\textrm{t}`$ to see where the transport is affected by turbulence.
The result for the turbulent viscosity should look like this:

![](../extradoc/ex_ff-pm-turb_diffusivity.png)

__Task B__:

Instead of computing the whole cross-section of a channel, you can use symmetric boundary conditions at the top boundary of your free flow domain by replacing all previous boundary conditions with
```c++
values.setAllSymmetry();
```

In addition, you have to remove the condition `onUpperBoundary_(globalPos)` from the `isOnWallAtPos(globalPos)`
and `initialAtPos(globalPos)` method.

__Task C__:

Choose `Surface With Edges` instead of `Surface` in the Paraview toolbar to see the discretization grid.
We will refine the grid now in order to better resolve the processes at the coupling interface.
Since not much is happening at the upper and lower boundaries of the whole domain, we want to keep the resolution low in these areas to save some computation time.

A grid refinement towards the interface is called _grading_.
Try different gradings by changing the values in the respective sections in the input file:
```c++
Grading0 = 1.0
Grading1 = 1.0
```

__Task D__:

For the grid convergence study, run various simulations with the following grading parameters:
```c++
* [Stokes.Grid] Grading1 = 1.2, [Darcy.Grid] Grading1 = -1.2
* [Stokes.Grid] Grading1 = 1.3, [Darcy.Grid] Grading1 = -1.3
* [Stokes.Grid] Grading1 = 1.4, [Darcy.Grid] Grading1 = -1.4
* [Stokes.Grid] Grading1 = 1.5, [Darcy.Grid] Grading1 = -1.5
* [Stokes.Grid] Grading1 = 1.6, [Darcy.Grid] Grading1 = -1.6
```

By changing the parameter `Problem.Name` for each grading factor, you avoid loosing the `.vtu` and `.pvd` files of the previous simulation runs.
Check the first lines of the output to see how the grading factors change the height of your grid cells.
Compare the velocity fields with Paraview.
