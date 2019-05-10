// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief A simple Darcy test problem (cell-centered finite volume method).
 */
#ifndef DUMUX_DARCY_SUBPROBLEM_HH
#define DUMUX_DARCY_SUBPROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/io/gnuplotinterface.hh>
#include <dumux/material/fluidsystems/1padapter.hh>
#include <dumux/material/fluidsystems/h2oair.hh>
#include <dumux/material/fluidmatrixinteractions/diffusivityconstanttortuosity.hh>

#include <dumux/porousmediumflow/problem.hh>

#include <dumux/porousmediumflow/1pnc/model.hh>
#include "../1pspatialparams.hh"

namespace Dumux
{
template <class TypeTag>
class DarcySubProblem;

namespace Properties
{
// Create new type tags
namespace TTag {
struct DarcyOnePNC { using InheritsFrom = std::tuple<OnePNC, CCTpfaModel>; };
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::DarcyOnePNC> { using type = Dumux::DarcySubProblem<TypeTag>; };

// The fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::DarcyOnePNC>
{
    using H2OAir = FluidSystems::H2OAir<GetPropType<TypeTag, Properties::Scalar>>;
    using type = FluidSystems::OnePAdapter<H2OAir, H2OAir::gasPhaseIdx>;
};

// Use moles
template<class TypeTag>
struct UseMoles<TypeTag, TTag::DarcyOnePNC> { static constexpr bool value = true; };

// Do not replace one equation with a total mass balance
template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::DarcyOnePNC> { static constexpr int value = 3; };

//! Use a model with constant tortuosity for the effective diffusivity
template<class TypeTag>
struct EffectiveDiffusivityModel<TypeTag, TTag::DarcyOnePNC>
{ using type = DiffusivityConstantTortuosity<GetPropType<TypeTag, Properties::Scalar>>; };
// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::DarcyOnePNC> { using type = Dune::YaspGrid<2>; };

// Set the spatial paramaters type
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::DarcyOnePNC> {
    using type = OnePSpatialParams<GetPropType<TypeTag, FVGridGeometry>, GetPropType<TypeTag, Scalar>>;
};

} // end namespace Properties

template <class TypeTag>
class DarcySubProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::FVGridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;

    // copy some indices for convenience
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    enum {
        // grid and world dimension
        dim = GridView::dimension,
        dimworld = GridView::dimensionworld,

        // primary variable indices
        conti0EqIdx = Indices::conti0EqIdx,
        pressureIdx = Indices::pressureIdx,
        phaseIdx = 0,
        transportCompIdx = 1
    };

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = Dune::FieldVector<Scalar, dimworld>;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using TimeLoopPtr = std::shared_ptr<TimeLoop<Scalar>>;

public:
    DarcySubProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                   std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(fvGridGeometry, "Darcy"), eps_(1e-7), couplingManager_(couplingManager)
    {
        moleFraction_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.MoleFraction");

        // initialize output file
        plotFluxes_ = getParamFromGroup<bool>(this->paramGroup(), "Problem.PlotFluxes", false);
        plotStorage_ = getParamFromGroup<bool>(this->paramGroup(), "Problem.PlotStorage", false);
        storageFileName_ = "storage_" + getParam<std::string>("Problem.Name") + "_" + this->name() + ".csv";
        storageFile_.open(storageFileName_);
        storageFile_ << "#Time[s]" << ";"
                     << "WaterMass[kg]" << ";"
                     << "WaterMassLoss[kg]" << ";"
                     << "EvaporationRate[mm/d]"
                     << std::endl;
    }

    /*!
     * \name Simulation steering
     */
    // \{

    /*!
     * \brief Initialize the problem.
     */
    template<class SolutionVector, class GridVariables>
    void init(const SolutionVector& curSol,
              const GridVariables& gridVariables)
    { }

    template<class SolutionVector, class GridVariables>
    void postTimeStep(const SolutionVector& curSol,
                      const GridVariables& gridVariables)

    {
        evaluateWaterMassStorageTerm(curSol, gridVariables);
        evaluateInterfaceFluxes(curSol, gridVariables);

        gnuplotStorage_.resetPlot();
        gnuplotStorage_.setDatafileSeparator(';');
        gnuplotStorage_.setXlabel("time [d]");
        gnuplotStorage_.setXRange(0.0, getParam<Scalar>("TimeLoop.TEnd"));
        gnuplotStorage_.setYlabel("evaporation rate [mm/d]");
        gnuplotStorage_.setOption("set yrange [0.0:]");
        gnuplotStorage_.setOption("set y2label 'cumulative mass loss'");
        gnuplotStorage_.setOption("set y2range [0.0:0.5]");
        gnuplotStorage_.setOption("set y2range [0.0:0.5]");
        gnuplotStorage_.addFileToPlot(storageFileName_, "using 1:4 with lines title 'evaporation rate'");
        gnuplotStorage_.addFileToPlot(storageFileName_, "using 1:3 axes x1y2 with lines title 'cumulative mass loss'");
        if (plotStorage_)
            gnuplotStorage_.plot("temp");
    }

    template<class SolutionVector, class GridVariables>
    Scalar evaluateWaterMassStorageTerm(const SolutionVector& curSol,
                                        const GridVariables& gridVariables)

    {
        // compute the mass in the entire domain
        Scalar waterMass = 0.0;

        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);

            auto elemVolVars = localView(gridVariables.curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, curSol);

            for (auto&& scv : scvs(fvGeometry))
            {
                // const auto& volVars = elemVolVars[scv];
                // insert calculation of the water mass here
                waterMass += 0.0;
            }
        }

        Scalar cumMassLoss = initialWaterContent_ - waterMass;
        Scalar evaporationRate = (lastWaterMass_ - waterMass) * 86400
                                 / (this->fvGridGeometry().bBoxMax()[0] - this->fvGridGeometry().bBoxMin()[0])
                                 / timeLoop_->timeStepSize();
        lastWaterMass_ = waterMass;

        storageFile_ << timeLoop_->time() << ";"
                     << waterMass << ";"
                     << cumMassLoss << ";"
                     << evaporationRate
                     << std::endl;

        return waterMass;
    }

    template<class SolutionVector, class GridVariables>
    void evaluateInterfaceFluxes(const SolutionVector& curSol,
                                 const GridVariables& gridVariables)

    {
        std::vector<Scalar> x;
        std::vector<Scalar> y;

        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);

            auto elemVolVars = localView(gridVariables.curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, curSol);

            for (auto&& scvf : scvfs(fvGeometry))
            {
                if (!couplingManager().isCoupledEntity(CouplingManager::darcyIdx, scvf))
                    continue;

                NumEqVector flux(0.0); // use "massCouplingCondition" from the couplingManager here

                x.push_back(scvf.center()[0]);
                y.push_back(flux[transportCompIdx]);
            }
        }

        gnuplotInterfaceFluxes_.resetPlot();
        gnuplotInterfaceFluxes_.setXlabel("x-position [m]");
        gnuplotInterfaceFluxes_.setXRange(this->fvGridGeometry().bBoxMin()[0], this->fvGridGeometry().bBoxMax()[0]);
        gnuplotInterfaceFluxes_.setYlabel("flux [kg/(m^2 s)]");
        gnuplotInterfaceFluxes_.setYRange(-5e-4, 0.0);
        gnuplotInterfaceFluxes_.setOption("set label 'time: " + std::to_string(timeLoop_->time()/86400.) + "d' at graph 0.8,0.8 ");
        std::string fluxFileName = "flux_" + std::to_string(timeLoop_->timeStepIndex()) +
                                   "_" + getParam<std::string>("Problem.Name") + "_" + this->name() + ".csv";
        gnuplotInterfaceFluxes_.addDataSetToPlot(x, y, fluxFileName, "with lines title 'water mass flux'");
        if (plotFluxes_)
            gnuplotInterfaceFluxes_.plot("flux_" + std::to_string(timeLoop_->timeStepIndex()));
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Return the temperature within the domain in [K].
     *
     */
    Scalar temperature() const
    { return 293.15; }
    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
      * \brief Specifies which kind of boundary condition should be
      *        used for which equation on a given boundary control volume.
      *
      * \param element The element
      * \param scvf The boundary sub control volume face
      */
    BoundaryTypes boundaryTypes(const Element& element, const SubControlVolumeFace& scvf) const
    {
        BoundaryTypes values;
        values.setAllNeumann();

        if (couplingManager().isCoupledEntity(CouplingManager::darcyIdx, scvf))
            values.setAllCouplingNeumann();

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a Neumann control volume.
     *
     * \param element The element for which the Neumann boundary condition is set
     * \param fvGeomentry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param scvf The boundary sub control volume face
     *
     * For this method, the \a values variable stores primary variables.
     */
    template<class ElementVolumeVariables>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);

        if (couplingManager().isCoupledEntity(CouplingManager::darcyIdx, scvf))
            values = couplingManager().couplingData().massCouplingCondition(element, fvGeometry, elemVolVars, scvf);

        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{
    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * \param element The element for which the source term is set
     * \param fvGeomentry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param scv The subcontrolvolume
     */
    template<class ElementVolumeVariables>
    NumEqVector source(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume& scv) const
    { return NumEqVector(0.0); }

    // \}

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param element The element
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        static const Scalar stokesPressure = getParamFromGroup<Scalar>("Stokes", "Problem.Pressure");

        PrimaryVariables values(0.0);
        values[Indices::pressureIdx] = stokesPressure;
        values[transportCompIdx] = moleFraction_;
        return values;
    }

    // \}

    //! Set the coupling manager
    void setCouplingManager(std::shared_ptr<CouplingManager> cm)
    { couplingManager_ = cm; }

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

    void setTimeLoop(TimeLoopPtr timeLoop)
    { timeLoop_ = timeLoop; }

private:
    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < this->fvGridGeometry().bBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > this->fvGridGeometry().bBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < this->fvGridGeometry().bBoxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->fvGridGeometry().bBoxMax()[1] - eps_; }

    Scalar eps_;
    Scalar moleFraction_;

    Scalar initialWaterContent_ = 0.0;
    Scalar lastWaterMass_ = 0.0;

    TimeLoopPtr timeLoop_;
    std::shared_ptr<CouplingManager> couplingManager_;

    std::string storageFileName_;
    std::ofstream storageFile_;
    bool plotFluxes_;
    bool plotStorage_;
    Dumux::GnuplotInterface<Scalar> gnuplotInterfaceFluxes_;
    Dumux::GnuplotInterface<Scalar> gnuplotStorage_;
};
} //end namespace Dumux

#endif //DUMUX_DARCY_SUBPROBLEM_HH
