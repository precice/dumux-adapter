// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup StokesDarcyCoupling
 * \copydoc Dumux::StokesDarcyCouplingData
 */

#ifndef DUMUX_FREEFLOW_SOLIDENERGY_COUPLINGDATA_HH
#define DUMUX_FREEFLOW_SOLIDENERGY_COUPLINGDATA_HH

#include <numeric>

#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>
#include <dumux/multidomain/couplingmanager.hh>

namespace Dumux
{
/*!
 * \ingroup StokesDarcyCoupling
 * \brief This structs holds a set of options which allow to modify the Stokes-Darcy
 *        coupling mechanism during runtime.
 */
struct FreeFlowHeatCouplingOptions {
    /*!
     * \brief Defines which kind of averanging of diffusion coefficiencients
     *        (moleculat diffusion or thermal conductance) at the interface
     *        between free flow and porous medium shall be used.
     */
    enum class DiffusionCoefficientAveragingType {
        harmonic,
        arithmethic,
        ffOnly,
        pmOnly
    };

    /*!
     * \brief Convenience function to convert user input given as std::string to the corresponding enum class used for chosing the type
     *        of averaging of the diffusion/conduction parameter at the interface between the two domains.
     */
    static DiffusionCoefficientAveragingType stringToEnum(
        const std::string &diffusionCoefficientAveragingType)
    {
        if (diffusionCoefficientAveragingType == "Harmonic")
            return DiffusionCoefficientAveragingType::harmonic;
        else if (diffusionCoefficientAveragingType == "Arithmethic")
            return DiffusionCoefficientAveragingType::arithmethic;
        else if (diffusionCoefficientAveragingType == "FreeFlowOnly")
            return DiffusionCoefficientAveragingType::ffOnly;
        else if (diffusionCoefficientAveragingType == "SolidOnly")
            return DiffusionCoefficientAveragingType::pmOnly;
        else
            DUNE_THROW(Dune::IOError,
                       "Unknown DiffusionCoefficientAveragingType");
    }
};

/*!
 * \ingroup StokesDarcyCoupling
 * \brief A base class which provides some common methods used for Stokes-Darcy coupling.
 */
template<class MDTraits, class CouplingManager>
class FreeFlowSolidEnergyCouplingData
{
    using Scalar = typename MDTraits::Scalar;

    template<std::size_t id>
    using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;
    template<std::size_t id>
    using FVGridGeometry =
        GetPropType<SubDomainTypeTag<id>, Properties::FVGridGeometry>;
    template<std::size_t id>
    using Element =
        typename FVGridGeometry<id>::GridView::template Codim<0>::Entity;
    template<std::size_t id>
    using FVElementGeometry = typename FVGridGeometry<id>::LocalView;
    template<std::size_t id>
    using SubControlVolumeFace =
        typename FVGridGeometry<id>::LocalView::SubControlVolumeFace;
    template<std::size_t id>
    using SubControlVolume =
        typename FVGridGeometry<id>::LocalView::SubControlVolume;
    template<std::size_t id>
    using Indices = typename GetPropType<SubDomainTypeTag<id>,
                                         Properties::ModelTraits>::Indices;
    template<std::size_t id>
    using ElementVolumeVariables =
        typename GetPropType<SubDomainTypeTag<id>,
                             Properties::GridVolumeVariables>::LocalView;
    template<std::size_t id>
    using ElementFaceVariables =
        typename GetPropType<SubDomainTypeTag<id>,
                             Properties::GridFaceVariables>::LocalView;
    template<std::size_t id>
    using VolumeVariables =
        typename GetPropType<SubDomainTypeTag<id>,
                             Properties::GridVolumeVariables>::VolumeVariables;
    template<std::size_t id>
    using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id>
    using FluidSystem =
        GetPropType<SubDomainTypeTag<id>, Properties::FluidSystem>;
    template<std::size_t id>
    using ModelTraits =
        GetPropType<SubDomainTypeTag<id>, Properties::ModelTraits>;

    static constexpr auto freeFlowIdx = CouplingManager::freeFlowIdx;
    static constexpr auto solidEnergyIdx = CouplingManager::solidEnergyIdx;

    static constexpr int enableEnergyBalance =
        GetPropType<SubDomainTypeTag<freeFlowIdx>,
                    Properties::ModelTraits>::enableEnergyBalance();
    static_assert(
        GetPropType<SubDomainTypeTag<solidEnergyIdx>,
                    Properties::ModelTraits>::enableEnergyBalance() ==
            enableEnergyBalance,
        "All submodels must both be either isothermal or non-isothermal");

public:
    using DiffusionCoefficientAveragingType =
        typename FreeFlowHeatCouplingOptions::DiffusionCoefficientAveragingType;

    FreeFlowSolidEnergyCouplingData(const CouplingManager &couplingmanager)
        : couplingManager_(couplingmanager)
    {
    }

    /*!
     * \brief Returns a reference to the coupling manager.
     */
    const CouplingManager &couplingManager() const { return couplingManager_; }

    /*!
     * \brief Returns the energy flux across the coupling boundary as seen from the Darcy domain.
     */
    template<bool isNI = enableEnergyBalance,
             typename std::enable_if_t<isNI, int> = 0>
    Scalar energyCouplingCondition(
        const Element<solidEnergyIdx> &element,
        const FVElementGeometry<solidEnergyIdx> &fvGeometry,
        const ElementVolumeVariables<solidEnergyIdx> &darcyElemVolVars,
        const SubControlVolumeFace<solidEnergyIdx> &scvf,
        const DiffusionCoefficientAveragingType diffCoeffAvgType =
            DiffusionCoefficientAveragingType::ffOnly) const
    {
        const auto &darcyContext =
            this->couplingManager().darcyCouplingContext(element, scvf);
        const auto &darcyVolVars = darcyElemVolVars[scvf.insideScvIdx()];
        const auto &stokesVolVars = darcyContext.volVars;

        return energyFlux_(solidEnergyIdx, freeFlowIdx, fvGeometry,
                           darcyContext.fvGeometry, scvf, darcyVolVars,
                           stokesVolVars, diffCoeffAvgType);
    }

    /*!
     * \brief Returns the energy flux across the coupling boundary as seen from the free-flow domain.
     */
    template<bool isNI = enableEnergyBalance,
             typename std::enable_if_t<isNI, int> = 0>
    Scalar energyCouplingCondition(
        const Element<freeFlowIdx> &element,
        const FVElementGeometry<freeFlowIdx> &fvGeometry,
        const ElementVolumeVariables<freeFlowIdx> &stokesElemVolVars,
        const ElementFaceVariables<freeFlowIdx> &stokesElemFaceVars,
        const SubControlVolumeFace<freeFlowIdx> &scvf,
        const DiffusionCoefficientAveragingType diffCoeffAvgType =
            DiffusionCoefficientAveragingType::ffOnly) const
    {
        const auto &stokesContext =
            this->couplingManager().stokesCouplingContext(element, scvf);
        const auto &stokesVolVars = stokesElemVolVars[scvf.insideScvIdx()];
        const auto &darcyVolVars = stokesContext.volVars;

        return energyFlux_(freeFlowIdx, solidEnergyIdx, fvGeometry,
                           stokesContext.fvGeometry, scvf, stokesVolVars,
                           darcyVolVars, diffCoeffAvgType);
    }

private:
    /*!
     * \brief Evaluate the energy flux across the interface.
     */
    template<std::size_t i,
             std::size_t j,
             bool isNI = enableEnergyBalance,
             typename std::enable_if_t<isNI, int> = 0>
    Scalar energyFlux_(
        Dune::index_constant<i> domainI,
        Dune::index_constant<j> domainJ,
        const FVElementGeometry<i> &insideFvGeometry,
        const FVElementGeometry<j> &outsideFvGeometry,
        const SubControlVolumeFace<i> &scvf,
        const VolumeVariables<i> &insideVolVars,
        const VolumeVariables<j> &outsideVolVars,
        const DiffusionCoefficientAveragingType diffCoeffAvgType) const
    {
        Scalar flux(0.0);

        const auto &insideScv = (*scvs(insideFvGeometry).begin());
        const auto &outsideScv = (*scvs(outsideFvGeometry).begin());

        flux += this->conductiveEnergyFlux_(domainI, domainJ, insideFvGeometry,
                                            outsideFvGeometry, scvf, insideScv,
                                            outsideScv, insideVolVars,
                                            outsideVolVars, diffCoeffAvgType);

        return flux;
    }

    /*!
     * \brief Returns the transmissibility used for either molecular diffusion or thermal conductivity.
     */
    template<std::size_t i, std::size_t j>
    Scalar transmissibility_(
        Dune::index_constant<i> domainI,
        Dune::index_constant<j> domainJ,
        const Scalar insideDistance,
        const Scalar outsideDistance,
        const Scalar avgQuantityI,
        const Scalar avgQuantityJ,
        const DiffusionCoefficientAveragingType diffCoeffAvgType) const
    {
        const Scalar totalDistance = insideDistance + outsideDistance;
        if (diffCoeffAvgType == DiffusionCoefficientAveragingType::harmonic) {
            return harmonicMean(avgQuantityI, avgQuantityJ, insideDistance,
                                outsideDistance) /
                   totalDistance;
        } else if (diffCoeffAvgType ==
                   DiffusionCoefficientAveragingType::arithmethic) {
            return arithmeticMean(avgQuantityI, avgQuantityJ, insideDistance,
                                  outsideDistance) /
                   totalDistance;
        } else if (diffCoeffAvgType ==
                   DiffusionCoefficientAveragingType::ffOnly)
            return domainI == freeFlowIdx ? avgQuantityI / totalDistance
                                          : avgQuantityJ / totalDistance;

        else  // diffCoeffAvgType == DiffusionCoefficientAveragingType::pmOnly)
            return domainI == solidEnergyIdx ? avgQuantityI / totalDistance
                                             : avgQuantityJ / totalDistance;
    }

    /*!
     * \brief Returns the distance between an scvf and the corresponding scv center.
     */
    template<class Scv, class Scvf>
    Scalar getDistance_(const Scv &scv, const Scvf &scvf) const
    {
        return (scv.dofPosition() - scvf.ipGlobal()).two_norm();
    }

    /*!
     * \brief Returns the conductive energy flux acorss the interface.
     */
    template<std::size_t i,
             std::size_t j,
             bool isNI = enableEnergyBalance,
             typename std::enable_if_t<isNI, int> = 0>
    Scalar conductiveEnergyFlux_(
        Dune::index_constant<i> domainI,
        Dune::index_constant<j> domainJ,
        const FVElementGeometry<i> &fvGeometryI,
        const FVElementGeometry<j> &fvGeometryJ,
        const SubControlVolumeFace<i> &scvfI,
        const SubControlVolume<i> &scvI,
        const SubControlVolume<j> &scvJ,
        const VolumeVariables<i> &volVarsI,
        const VolumeVariables<j> &volVarsJ,
        const DiffusionCoefficientAveragingType diffCoeffAvgType) const
    {
        const Scalar insideDistance = getDistance_(scvI, scvfI);
        const Scalar outsideDistance = getDistance_(scvJ, scvfI);

        const Scalar deltaT = volVarsJ.temperature() - volVarsI.temperature();
        const Scalar tij =
            transmissibility_(domainI, domainJ, insideDistance, outsideDistance,
                              thermalConductivity_(volVarsI, fvGeometryI, scvI),
                              thermalConductivity_(volVarsJ, fvGeometryJ, scvJ),
                              diffCoeffAvgType);

        return -tij * deltaT;
    }

    /*!
     * \brief Returns the effective thermal conductivity (lumped parameter) within the porous medium.
     */
    template<bool isNI = enableEnergyBalance,
             typename std::enable_if_t<isNI, int> = 0>
    Scalar thermalConductivity_(
        const VolumeVariables<solidEnergyIdx> &volVars,
        const FVElementGeometry<solidEnergyIdx> &fvGeometry,
        const SubControlVolume<solidEnergyIdx> &scv) const
    {
        using ThermalConductivityModel =
            GetPropType<SubDomainTypeTag<solidEnergyIdx>,
                        Properties::ThermalConductivityModel>;
        const auto &problem = this->couplingManager().problem(solidEnergyIdx);
        return ThermalConductivityModel::effectiveThermalConductivity(
            volVars, problem.spatialParams(),
            fvGeometry.fvGridGeometry().element(scv), fvGeometry, scv);
    }

    /*!
     * \brief Returns the thermal conductivity of the fluid phase within the free flow domain.
     */
    template<bool isNI = enableEnergyBalance,
             typename std::enable_if_t<isNI, int> = 0>
    Scalar thermalConductivity_(
        const VolumeVariables<freeFlowIdx> &volVars,
        const FVElementGeometry<freeFlowIdx> &fvGeometry,
        const SubControlVolume<freeFlowIdx> &scv) const
    {
        return volVars.effectiveThermalConductivity();
    }

    const CouplingManager &couplingManager_;
};

}  // end namespace Dumux

#endif  // DUMUX_STOKES_DARCY_COUPLINGDATA_HH
