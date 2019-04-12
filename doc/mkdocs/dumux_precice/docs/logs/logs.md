# Overview over the last meetings

**Next meeting**: 2019-??-??

## 2019-04-12

Side note: My class `PreciceWrapper` is in fact a singleton.

- We found out that certain questions have come up regarding the test case

    - Maybe double check wall boundary conditions
    - Inflow: Verify that it is a rectangular flow profile (slip walls)
    - Does the fluid have to be incompressible?
    - Are the flow parameters depending on the solution? This means, are flow parameters like the viscosity depending on the temperature?
    - What are the flow quantities. 
    - Where are the values prescribed in OpenFOAM. Are these actual boundary values or are they cell centered?
    - How is the coupling carried out? What values are used for the heat flux? What temperatures do we use and which conductivity?


Question for us:
**What quantities do we need?**

Known parameters

| Parameter | value |
| --------- | ----- |
| Reynolds  |  500  |
| Prandtl   |  0.1  |
| k         |  1    |

**Note**: k=1 means that k_s = k_f = 1. The heat conductivity of fluid and solid are the same.

**Fluid**

- Inflow velocity: 
- Inflow temperature: 
- Viscosity (given Reynolds number)

**Solid**

- Conductivity: Comes from the Prandtl number. We need to now c_p and density rho of the fluid

## 2019-04-03

Participants: Kilian, Ned, Alex

### Contents

- Discussion on the test case and how to set it up. 
- Some discussion about tools (`MkDocs`, `pandoc`, `sphinx` etc.). 
- Alex should come to one of the DuMuX meetings some time. Maybe present some of the tools.


### Outcome/Plans: 

- Create GitLab repository
- Set up [conjugate heat transfer test case](../dumux_test_case.md#conjugate-heat-transfer) in DuMuX. This will be done twice:
    1. Monolithic setup with DuMuX
    1. Two seperate DuMuX solvers that will be coupled with preCICE afterwards.
- Alex should push his notes about the `dumux-lecture` to the repository. 
- Checkpointing of the solution for preCICE has to be done "by hand". There is no functionality for this currently in DuMuX.
- Integrate the preCICE-DuMuX-adapter in the DuMuX repository. This keeps it close to potential users and it will not be forgotten in the preCICE repository.