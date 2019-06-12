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
 * \brief A test problem for the coupled Stokes/Darcy problem (1p)
 */
#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/istl/io.hh>

#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/valgrind.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/defaultusagemessage.hh>

#include <dumux/linear/amgbackend.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>

#include <dumux/discretization/method.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager.hh>

#include "../monolithic/pmproblem.hh"

#include "../../precice-adapter/include/preciceadapter.hh"

 /*!
  * \brief Returns the pressure at the interface using Darcy's law for reconstruction
  */
 template<class Problem, class Element, class FVElementGeometry, class ElementVolumeVariables, class SubControlVolumeFace>
 auto pressureAtInterface(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf)
 {
     using Scalar = typename ElementVolumeVariables::VolumeVariables::PrimaryVariables::value_type;
     const auto& volVars = elemVolVars[scvf.insideScvIdx()];

     const Scalar boundaryFlux = problem.neumann(element, fvGeometry, elemVolVars, scvf);

     // v = -K/mu * (gradP + rho*g)
     auto velocity = scvf.unitOuterNormal();
     velocity *= boundaryFlux; // TODO check sign

    const Scalar ccPressure = volVars.pressure();
    const Scalar mobility = volVars.mobility();
    const Scalar density = volVars.density();
    const auto K = volVars.permeability();

    // v = -kr/mu*K * (gradP + rho*g) = -mobility*K * (gradP + rho*g)
    const auto alpha = Dumux::vtmv( scvf.unitOuterNormal(), K, problem.gravity() );

    auto distanceVector = scvf.center() - element.geometry().center();
    distanceVector /= distanceVector.two_norm2();
    const Scalar ti = Dumux::vtmv(distanceVector, K, scvf.unitOuterNormal());

    return (1/mobility * (scvf.unitOuterNormal() * velocity) + density * alpha)/ti
           + ccPressure;
 }

 template<class Problem, class GridVariables, class SolutionVector>
 void setInterfacePressures(const Problem& problem,
                            const GridVariables& gridVars,
                            const SolutionVector& sol)
 {
   const auto& fvGridGeometry = problem.fvGridGeometry();
   auto fvGeometry = localView(fvGridGeometry);
   auto elemVolVars = localView(gridVars.curGridVolVars());

   auto& couplingInterface = precice_adapter::PreciceAdapter::getInstance();
   const auto pressureId = couplingInterface.getIdFromName( "Pressure" );

   for (const auto& element : elements(fvGridGeometry.gridView()))
   {
     fvGeometry.bindElement(element);
     elemVolVars.bindElement(element, fvGeometry, sol);

     for (const auto& scvf : scvfs(fvGeometry))
     {

       if ( couplingInterface.isCoupledEntity( scvf.index() ) )
       {
         //TODO: What to do here?
         const double p = pressureAtInterface(problem, element, fvGridGeometry, elemVolVars, scvf);
         couplingInterface.writeQuantityOnFace( pressureId, scvf.index(), p );
       }
     }
   }

 }

int main(int argc, char** argv) try
{
    using namespace Dumux;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    using DarcyTypeTag = Properties::TTag::DarcyOneP;

    using DarcyGridManager = Dumux::GridManager<GetPropType<DarcyTypeTag, Properties::Grid>>;
    DarcyGridManager darcyGridManager;
    darcyGridManager.init("Darcy"); // pass parameter group

    // we compute on the leaf grid view
    const auto& darcyGridView = darcyGridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using DarcyFVGridGeometry = GetPropType<DarcyTypeTag, Properties::FVGridGeometry>;
    auto darcyFvGridGeometry = std::make_shared<DarcyFVGridGeometry>(darcyGridView);
    darcyFvGridGeometry->update();

    using DarcyProblem = GetPropType<DarcyTypeTag, Properties::Problem>;
    auto darcyProblem = std::make_shared<DarcyProblem>(darcyFvGridGeometry);

    // the solution vector
    GetPropType<DarcyTypeTag, Properties::SolutionVector> sol;
    sol.resize(darcyFvGridGeometry->numDofs());

    // Initialize preCICE.Tell preCICE about:
    // - Name of solver
    // - What rank of how many ranks this instance is
    // Configure preCICE. For now the config file is hardcoded.
    //couplingInterface.createInstance( "darcy", mpiHelper.rank(), mpiHelper.size() );
    std::string preciceConfigFilename = "precice-config.xml";
    if (argc == 3)
      preciceConfigFilename = argv[2];

    auto& couplingInterface =
        precice_adapter::PreciceAdapter::getInstance();
    couplingInterface.announceSolver( "Darcy", preciceConfigFilename,
                                      mpiHelper.rank(), mpiHelper.size() );

    const int dim = couplingInterface.getDimensions();
    std::cout << dim << "  " << int(DarcyFVGridGeometry::GridView::dimension) << std::endl;
    if (dim != int(DarcyFVGridGeometry::GridView::dimension))
        DUNE_THROW(Dune::InvalidStateException, "Dimensions do not match");

    // GET mesh corodinates
    const double xMin = getParamFromGroup<std::vector<double>>("Darcy", "Grid.Positions0")[0];
    const double xMax = getParamFromGroup<std::vector<double>>("Darcy", "Grid.Positions0").back();
    std::vector<double> coords; //( dim * vertexSize );
    std::vector<int> coupledScvfIndices;

    for (const auto& element : elements(darcyGridView))
    {
        auto fvGeometry = localView(*darcyFvGridGeometry);
        fvGeometry.bindElement(element);

        for (const auto& scvf : scvfs(fvGeometry))
        {
            static constexpr auto eps = 1e-7;
            const auto& pos = scvf.center();
            if (pos[1] > darcyFvGridGeometry->bBoxMax()[1] - eps)
            {
                if (pos[0] > xMin - eps && pos[0] < xMax + eps)
                {
                  coupledScvfIndices.push_back(scvf.index());
                    for (const auto p : pos)
                        coords.push_back(p);
                }
            }
        }
    }

    const auto numberOfPoints = coords.size() / dim;
    const double preciceDt = couplingInterface.setMeshAndInitialize( "DarcyMesh",
                                                                     numberOfPoints,
                                                                     coords,
                                                                     coupledScvfIndices );

    const auto velocityId = couplingInterface.announceQuantity( "Velocity" );
    const auto pressureId = couplingInterface.announceQuantity( "Pressure" );

    darcyProblem->updatePreciceDataIds();

    darcyProblem->applyInitialSolution(sol);

    // the grid variables
    using DarcyGridVariables = GetPropType<DarcyTypeTag, Properties::GridVariables>;
    auto darcyGridVariables = std::make_shared<DarcyGridVariables>(darcyProblem, darcyFvGridGeometry);
    darcyGridVariables->init(sol);

    // intialize the vtk output module
    const auto darcyName = getParam<std::string>("Problem.Name") + "_" + darcyProblem->name();

    VtkOutputModule<DarcyGridVariables, GetPropType<DarcyTypeTag, Properties::SolutionVector>> darcyVtkWriter(*darcyGridVariables, sol, darcyProblem->name());
    using DarcyVelocityOutput = GetPropType<DarcyTypeTag, Properties::VelocityOutput>;
    darcyVtkWriter.addVelocityOutput(std::make_shared<DarcyVelocityOutput>(*darcyGridVariables));
    GetPropType<DarcyTypeTag, Properties::IOFields>::initOutputModule(darcyVtkWriter);
    darcyVtkWriter.write(0.0);


    if ( couplingInterface.hasToWriteInitialData() )
    {
      //TODO
      //couplingInterface.writeQuantityVector(velocityId);
      setInterfacePressures( *darcyProblem, *darcyGridVariables, sol );
      couplingInterface.writeQuantityToOtherSolver( pressureId );
      couplingInterface.announceInitialDataWritten();
    }
    couplingInterface.initializeData();

    // the assembler for a stationary problem
    using Assembler = FVAssembler<DarcyTypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(darcyProblem, darcyFvGridGeometry, darcyGridVariables);

    // the linear solver
    using LinearSolver = UMFPackBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    auto dt = preciceDt;
    auto sol_checkpoint = sol;

    while ( couplingInterface.isCouplingOngoing() )
    {
        if ( couplingInterface.hasToWriteIterationCheckpoint() )
        {
            //DO CHECKPOINTING
            sol_checkpoint = sol;
            couplingInterface.announceIterationCheckpointWritten();
        }

        // TODO
        couplingInterface.readQuantityFromOtherSolver( velocityId );

        // solve the non-linear system
        nonLinearSolver.solve(sol);
        setInterfacePressures( *darcyProblem, *darcyGridVariables, sol );
        couplingInterface.writeQuantityToOtherSolver( pressureId );

        const double preciceDt = couplingInterface.advance( dt );
        dt = std::min( preciceDt, dt );

        if ( couplingInterface.hasToReadIterationCheckpoint() )
        {
            //Read checkpoint
            sol = sol_checkpoint;
            darcyGridVariables->update(sol);
            darcyGridVariables->advanceTimeStep();
            //darcyGridVariables->init(sol);
            couplingInterface.announceIterationCheckpointRead();
        }
        else // coupling successful
        {
          // write vtk output
          darcyVtkWriter.write(1.0);
        }

    }
    // write vtk output
    darcyVtkWriter.write(1.0);

    couplingInterface.finalize();

    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////

    // print dumux end message
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
} // end main
catch (Dumux::ParameterException &e)
{
    std::cerr << std::endl << e << " ---> Abort!" << std::endl;
    return 1;
}
catch (Dune::DGFException & e)
{
    std::cerr << "DGF exception thrown (" << e <<
                 "). Most likely, the DGF file name is wrong "
                 "or the DGF file is corrupted, "
                 "e.g. missing hash at end of file or wrong number (dimensions) of entries."
                 << " ---> Abort!" << std::endl;
    return 2;
}
catch (Dune::Exception &e)
{
    std::cerr << "Dune reported error: " << e << " ---> Abort!" << std::endl;
    return 3;
}
catch (...)
{
    std::cerr << "Unknown exception thrown! ---> Abort!" << std::endl;
    return 4;
}
