# dumux-precice-wrapper

A preCICE adapter toolbox for DuMuX

## Ideas for a tighter integration

- Minimum requirement of what should work
    - Coupling of conjugate heat flow problem (flow over heated plate)
    - Tutorial test case with horizontal flow over a porous medium. See https://git.iws.uni-stuttgart.de/dumux-repositories/dumux-course/tree/master/exercises/exercise-coupling-ff-pm
        - Coupling in both directions should be done
- Have specific function depending on coupling?!
    - `Pressure`
    - `Velocity`
    - `Temperature`
    - `HeatFlux`
    - `NormalVelocity`
    - `TangentialVelocity` <-- Is this necessary?
- Derive from coupling manager in DuMuX.
    - What can I do 1 to 1 as in the coupling manager?
    - How can we iterate several times for a time step? This should be done "internally" so that the user does not see itz
- Template the adapter wrt to dimension

