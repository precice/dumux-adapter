#ifndef HELPERFUNCTIONS_HH
#define HELPERFUNCTIONS_HH

#include <type_traits>
// DuMuX
#include <dumux/common/properties.hh>
// preCICE
#include "preciceadapter.hh"

namespace helperfunctions
{
enum ProblemType { FreeFlow, Heat };

template<class Problem, class GridVariables, class SolutionVector>
void setBoundaryHeatFluxes(const Problem &problem,
                           const GridVariables &gridVars,
                           const SolutionVector &sol)
{
    const auto &fvGridGeometry = problem.fvGridGeometry();
    auto fvGeometry = localView(fvGridGeometry);
    auto elemVolVars = localView(gridVars.curGridVolVars());
    auto elemFaceVars = localView(gridVars.curGridFaceVars());

    auto &couplingInterface = precice_adapter::PreciceAdapter::getInstance();

    for (const auto &element : elements(fvGridGeometry.gridView())) {
        fvGeometry.bindElement(element);
        elemVolVars.bindElement(element, fvGeometry, sol);
        elemFaceVars.bindElement(element, fvGeometry, sol);

        for (const auto &scvf : scvfs(fvGeometry)) {
            if (couplingInterface.isCoupledEntity(scvf.index())) {
                //TODO: Actually writes temperature
                const auto heatFlux = problem.neumann(
                    element, fvGeometry, elemVolVars, elemFaceVars, scvf)[3];
                couplingInterface.writeHeatFluxOnFace(scvf.index(), heatFlux);
            }
        }
    }
}

template<class Problem, class GridVariables, class SolutionVector>
void printCellCenterTemperatures(const Problem &problem,
                                 const GridVariables &gridVars,
                                 const SolutionVector &sol)
{
    const auto &fvGridGeometry = problem.fvGridGeometry();
    auto fvGeometry = localView(fvGridGeometry);
    auto elemVolVars = localView(gridVars.curGridVolVars());
    auto elemFaceVars = localView(gridVars.curGridFaceVars());

    auto &couplingInterface = precice_adapter::PreciceAdapter::getInstance();

    for (const auto &element : elements(fvGridGeometry.gridView())) {
        fvGeometry.bindElement(element);
        elemVolVars.bindElement(element, fvGeometry, sol);
        elemFaceVars.bindElement(element, fvGeometry, sol);

        for (const auto &scvf : scvfs(fvGeometry)) {
            if (couplingInterface.isCoupledEntity(scvf.index())) {
                //TODO: Actually writes temperature
                const auto &scv = fvGeometry.scv(scvf.insideScvIdx());
                const auto &volVars = elemVolVars[scv];

                std::cout << "Temperature on cell center is: "
                          << volVars.temperature() << std::endl;
                //              const auto heatFlux = problem.neumann( element, fvGeometry, elemVolVars, elemFaceVars, scvf )[3];
                //              couplingInterface.writeHeatFluxOnFace( scvf.index(), heatFlux );
            }
        }
    }
}

namespace freeflow
{
template<class ThermalConductivityModel,
         class Problem,
         class FVElementGeometry,
         class ElementVolumeVariables>
auto reconstructBoundaryTemperature(
    const Problem &problem,
    const typename FVElementGeometry::FVGridGeometry::GridView::template Codim<
        0>::Entity &element,
    const FVElementGeometry &fvGeometry,
    const ElementVolumeVariables &elemVolVars,
    const typename FVElementGeometry::SubControlVolumeFace &scvf)
{
    using Scalar = typename ElementVolumeVariables::VolumeVariables::
        PrimaryVariables::value_type;

    const auto &scv = fvGeometry.scv(scvf.insideScvIdx());
    const auto &volVars = elemVolVars[scv];
    const Scalar cellCenterTemperature = volVars.temperature();
    const Scalar distance = (scvf.center() - scv.center()).two_norm();

    const Scalar insideLambda = volVars.effectiveThermalConductivity();

    const Scalar qHeat =
        problem.neumann(element, fvGeometry, elemVolVars, scvf);
    // q = -lambda * (t_face - t_cc) / dx
    // t_face = -q * dx / lambda + t_cc
    return -qHeat * distance / insideLambda + cellCenterTemperature;
}
}  // namespace freeflow

template<class Problem, class GridVariables, class SolutionVector>
void printCellCenterTemperatures(const Problem &problem,
                                 const GridVariables &gridVars,
                                 const SolutionVector &sol)
{
    const auto &fvGridGeometry = problem.fvGridGeometry();
    auto fvGeometry = localView(fvGridGeometry);
    auto elemVolVars = localView(gridVars.curGridVolVars());
    auto elemFaceVars = localView(gridVars.curGridFaceVars());

    auto &couplingInterface = precice_adapter::PreciceAdapter::getInstance();

    for (const auto &element : elements(fvGridGeometry.gridView())) {
        fvGeometry.bindElement(element);
        elemVolVars.bindElement(element, fvGeometry, sol);
        elemFaceVars.bindElement(element, fvGeometry, sol);

        for (const auto &scvf : scvfs(fvGeometry)) {
            if (couplingInterface.isCoupledEntity(scvf.index())) {
                //TODO: Actually writes temperature
                const auto &scv = fvGeometry.scv(scvf.insideScvIdx());
                const auto &volVars = elemVolVars[scv];

                std::cout << "Temperature on cell center is: "
                          << volVars.temperature() << std::endl;
                //              const auto heatFlux = problem.neumann( element, fvGeometry, elemVolVars, elemFaceVars, scvf )[3];
                //              couplingInterface.writeHeatFluxOnFace( scvf.index(), heatFlux );
            }
        }
    }
}

template<class ThermalConductivityModel,
         class Problem,
         class FVElementGeometry,
         class ElementVolumeVariables>
auto getThermalConductivity(
    const Problem &problem,
    const typename FVElementGeometry::FVGridGeometry::GridView::template Codim<
        0>::Entity &element,
    const FVElementGeometry &fvGeometry,
    const ElementVolumeVariables &elemVolVars,
    const typename FVElementGeometry::SubControlVolumeFace &scvf,
    const ProblemType problemType)
{
    const auto &scv = fvGeometry.scv(scvf.insideScvIdx());
    const auto &volVars = elemVolVars[scv];

    switch (problemType) {
        case FreeFlow:
            return ThermalConductivityModel::effectiveThermalConductivity(
                volVars, problem.spatialParams(), element, fvGeometry, scv);
        case Heat:
            return ThermalConductivityModel::thermalConductivity(
                volVars, problem.spatialParams(), element, fvGeometry, scv);
    }
    throw std::runtime_error(
        "helperfunctions.hh: getThermalConductivity was called with wrong "
        "problemType ");
}

namespace freeflow
{
template<class ThermalConductivityModel,
         class Problem,
         class FVElementGeometry,
         class ElementVolumeVariables>
auto reconstructBoundaryTemperature(
    const Problem &problem,
    const typename FVElementGeometry::FVGridGeometry::GridView::template Codim<
        0>::Entity &element,
    const FVElementGeometry &fvGeometry,
    const ElementVolumeVariables &elemVolVars,
    const typename FVElementGeometry::SubControlVolumeFace &scvf)
{
    using Scalar = typename ElementVolumeVariables::VolumeVariables::
        PrimaryVariables::value_type;

    const auto &scv = fvGeometry.scv(scvf.insideScvIdx());
    const auto &volVars = elemVolVars[scv];
    const Scalar cellCenterTemperature = volVars.temperature();
    const Scalar distance = (scvf.center() - scv.center()).two_norm();

    //        const Scalar insideLambda = ThermalConductivityModel::thermalConductivity(volVars,
    //                                                                                           problem.spatialParams(),
    //                                                                                           element,
    //                                                                                           fvGeometry,
    //                                                                                           scv);

    //        const Scalar qHeat = problem.neumann(element, fvGeometry, elemVolVars, scvf);
    //        // q = -lambda * (t_face - t_cc) / dx
    //        // t_face = -q * dx / lambda + t_cc
    //        return -qHeat * distance / insideLambda + cellCenterTemperature;
    return 0;
}

template<class ThermalConductivityModel,
         class Problem,
         class GridVariables,
         class SolutionVector>
void setBoundaryTemperatures(const Problem &problem,
                             const GridVariables &gridVars,
                             const SolutionVector &sol)
{
    const auto &fvGridGeometry = problem.fvGridGeometry();
    auto fvGeometry = localView(fvGridGeometry);
    auto elemVolVars = localView(gridVars.curGridVolVars());

    auto &couplingInterface = precice_adapter::PreciceAdapter::getInstance();

    for (const auto &element : elements(fvGridGeometry.gridView())) {
        fvGeometry.bindElement(element);
        elemVolVars.bindElement(element, fvGeometry, sol);

        for (const auto &scvf : scvfs(fvGeometry)) {
            if (couplingInterface.isCoupledEntity(scvf.index())) {
                //TODO: Actually writes temperature
                couplingInterface.writeTemperatureOnFace(
                    scvf.index(),
                    reconstructBoundaryTemperature<ThermalConductivityModel>(
                        problem, element, fvGeometry, elemVolVars, scvf));
            }
        }
    }
}

}  // namespace freeflow

namespace heat
{
template<class ThermalConductivityModel,
         class Problem,
         class FVElementGeometry,
         class ElementVolumeVariables>
auto reconstructBoundaryTemperature(
    const Problem &problem,
    const typename FVElementGeometry::FVGridGeometry::GridView::template Codim<
        0>::Entity &element,
    const FVElementGeometry &fvGeometry,
    const ElementVolumeVariables &elemVolVars,
    const typename FVElementGeometry::SubControlVolumeFace &scvf)
{
    using Scalar = typename ElementVolumeVariables::VolumeVariables::
        PrimaryVariables::value_type;

    const auto &scv = fvGeometry.scv(scvf.insideScvIdx());
    const auto &volVars = elemVolVars[scv];
    const Scalar cellCenterTemperature = volVars.temperature();
    const Scalar distance = (scvf.center() - scv.center()).two_norm();

    const Scalar insideLambda =
        ThermalConductivityModel::effectiveThermalConductivity(
            volVars, problem.spatialParams(), element, fvGeometry, scv);

    const Scalar qHeat =
        problem.neumann(element, fvGeometry, elemVolVars, scvf);
    // q = -lambda * (t_face - t_cc) / dx
    // t_face = -q * dx / lambda + t_cc
    return -qHeat * distance / insideLambda + cellCenterTemperature;
}

template<class ThermalConductivityModel,
         class Problem,
         class GridVariables,
         class SolutionVector>
void setBoundaryTemperatures(const Problem &problem,
                             const GridVariables &gridVars,
                             const SolutionVector &sol)
{
    const auto &fvGridGeometry = problem.fvGridGeometry();
    auto fvGeometry = localView(fvGridGeometry);
    auto elemVolVars = localView(gridVars.curGridVolVars());

    auto &couplingInterface = precice_adapter::PreciceAdapter::getInstance();

    for (const auto &element : elements(fvGridGeometry.gridView())) {
        fvGeometry.bindElement(element);
        elemVolVars.bindElement(element, fvGeometry, sol);

        for (const auto &scvf : scvfs(fvGeometry)) {
            if (couplingInterface.isCoupledEntity(scvf.index())) {
                //TODO: Actually writes temperature
                couplingInterface.writeTemperatureOnFace(
                    scvf.index(),
                    reconstructBoundaryTemperature<ThermalConductivityModel>(
                        problem, element, fvGeometry, elemVolVars, scvf));
            }
        }
    }
}
}  // namespace heat

}  // namespace helperfunctions

#endif  // HELPERFUNCTIONS_HH
