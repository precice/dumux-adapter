#ifndef OUTPUTHELPER_HH
#define OUTPUTHELPER_HH

#include <fstream>
#include <iomanip>

namespace outputhelper
{
namespace monolithic
{
template<class FluxVariables,
         class CouplingManager,
         class Problem,
         class GridVariables,
         class SolutionVector>
std::tuple<double, double, double> writeVelocitiesOnInterfaceToFile(
    const std::string &filename,
    const CouplingManager &couplingManager,
    const Problem &problem,
    const GridVariables &gridVars,
    const SolutionVector &sol)
{
    const auto &gridGeometry = problem.gridGeometry();
    auto fvGeometry = localView(gridGeometry);
    auto elemVolVars = localView(gridVars.curGridVolVars());
    auto elemFluxVarsCache = localView(gridVars.gridFluxVarsCache());

    std::ofstream ofs(filename + ".csv",
                      std::ofstream::out | std::ofstream::trunc);
    ofs << "x,y,";
    ofs << "velocityY"
        << "\n";
    double min = std::numeric_limits<double>::max();
    double max = std::numeric_limits<double>::min();
    double sum = 0.;
    for (const auto &element : elements(gridGeometry.gridView())) {
        fvGeometry.bind(element);
        elemVolVars.bind(element, fvGeometry, sol);
        elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);

        for (const auto &scvf : scvfs(fvGeometry)) {
            if (couplingManager.isCoupledEntity(CouplingManager::darcyIdx,
                                                scvf)) {
                const auto &pos = scvf.center();
                for (int i = 0; i < 2; ++i) {
                    ofs << pos[i] << ",";
                }
                const double v =
                    couplingManager.couplingData().massCouplingCondition(
                        element, fvGeometry, elemVolVars, scvf) /
                    1e3;
                max = std::max(v, max);
                min = std::min(v, min);
                sum += v;
                const int prec = ofs.precision();
                ofs << std::setprecision(std::numeric_limits<double>::digits10 +
                                         1)
                    << v << "\n";
                ofs.precision(prec);
                //        ofs << v / 1e3 << "\n";
            }
        }
    }

    ofs.close();
    return std::make_tuple(min, max, sum);
}

template<class CouplingManager, class Problem, class SolutionVector>
std::tuple<double, double, double> writeStokesVelocitiesOnInterfaceToFile(
    const std::string &filename,
    const CouplingManager &couplingManager,
    const Problem &problem,
    const SolutionVector &sol)
{
    const auto &gridGeometry = problem.gridGeometry();
    auto fvGeometry = localView(gridGeometry);
    using GridGeometry = std::decay_t<decltype(gridGeometry)>;

    std::ofstream ofs(filename + ".csv",
                      std::ofstream::out | std::ofstream::trunc);
    ofs << "x,y,";
    ofs << "velocityY"
        << "\n";

    double min = std::numeric_limits<double>::max();
    double max = std::numeric_limits<double>::min();
    double sum = 0.;
    for (const auto &element : elements(gridGeometry.gridView())) {
        fvGeometry.bind(element);

        for (const auto &scvf : scvfs(fvGeometry)) {
            if (couplingManager.isCoupledEntity(CouplingManager::stokesIdx,
                                                scvf)) {
                const auto &pos = scvf.center();
                for (int i = 0; i < 2; ++i) {
                    ofs << pos[i] << ",";
                }
                const double v = sol[scvf.dofIndex()];
                max = std::max(v, max);
                min = std::min(v, min);
                sum += v;
                const int prec = ofs.precision();
                ofs << std::setprecision(std::numeric_limits<double>::digits10 +
                                         1)
                    << v << "\n";
                ofs.precision(prec);
            }
        }
    }

    ofs.close();
    return std::make_tuple(min, max, sum);
}

}  // namespace monolithic

}  // namespace outputhelper

#endif  // OUTPUTHELPER_HH
