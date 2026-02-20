#include <config.h>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <vector>
#include <algorithm>
#include <iostream>

#include "dumux-precice/couplingadapter.hh"

int main(int argc, char** argv)
{
    auto& mpi = Dune::MPIHelper::instance(argc, argv);
    Dune::TestSuite t;

    static constexpr int dim = 2;

    // Build a 8*8 grid;
    Dune::FieldVector<double, dim> L(1.0);
    std::array<int, dim> N{8, 8};
    std::bitset<dim> periodic(false);
    int overlap = 1;

    Dune::YaspGrid<dim> grid(L, N, periodic, overlap);
    auto gv = grid.leafGridView();
    const auto& iset = gv.indexSet();

    std::vector<int> ids;
    std::vector<double> pos;
    for (const auto& e : elements(gv)) {
        ids.push_back(static_cast<int>(iset.index(e)));
        const auto c = e.geometry().center();
        pos.push_back(c[0]);
        pos.push_back(c[1]);
    }

    // Reference: local interior element indices
    std::vector<int> refInterior;
    for (const auto& e : elements(gv)) {
        if (e.partitionType() == Dune::InteriorEntity)
            refInterior.push_back(static_cast<int>(iset.index(e)));
    }

    std::cout << "\n=== Rank " << mpi.rank() << " ===\n";
    std::cout << "Total local elements: " << gv.size(0) << "\n";
    std::cout << "Initial ids.size(): " << ids.size() << "\n";
    std::cout << "Initial pos.size(): " << pos.size() << "\n";
    std::cout << "Reference interior count: " << refInterior.size() << "\n";

    std::cout << "Initial local IDs on rank "<< mpi.rank()<<":" << "\n";
    for (auto id : ids)
        std::cout << id << " ";
    std::cout << "\n";

    
    auto &couplingParticipant = Dumux::Precice::CouplingAdapter::getInstance();
    couplingParticipant.filterInteriorEntities(gv, ids, pos);

    std::cout << "After filtering ids.size(): " << ids.size() << " on rank "<< mpi.rank() << "\n";
    std::cout << "After filtering pos.size(): " << pos.size() << " on rank "<< mpi.rank() << "\n";

    std::cout << "Local IDs after filtering on rank "<< mpi.rank()<<":" << "\n";
    for (auto id : ids)
        std::cout << id << " ";
    std::cout << "\n";


    t.check(pos.size() == ids.size() * dim)
        << "pos size inconsistent after filtering";

    t.check(ids.size() == refInterior.size())
        << "filtered id count mismatch";

    for (auto id : ids) {
        t.check(std::find(refInterior.begin(),
                          refInterior.end(),
                          id) != refInterior.end())
            << "filtered id not interior: " << id;
    }

    std::cout << "Rank " << mpi.rank() << " test finished.\n";

    return t.exit();
}