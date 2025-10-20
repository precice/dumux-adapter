#ifndef SOLVERSTATE_HH
#define SOLVERSTATE_HH

#include <dumux/common/properties.hh>
#include <ostream>
#include <precice/precice.hpp>
#include <string>

/*!
 * @brief Namespace of dumux-precice
 *
 */
namespace Dumux::Precice
{
struct SolverStateBase {
    virtual ~SolverStateBase() = default;
    virtual void writeState() = 0;
    virtual void readState() = 0;
};
/*!
    * @brief A class to store and provide the state of the solver while checkpointing, one SolutionVector object is supported.
    */
template<class SolutionVector>
class SolverStateOnly : public SolverStateBase
{
private:
    SolutionVector *x_;
    SolutionVector xCheckpoint_;

public:
    SolverStateOnly(SolutionVector &x) : x_(&x), xCheckpoint_(*x_) {}

    void writeState() override { xCheckpoint_ = *x_; }

    void readState() override { *x_ = xCheckpoint_; }
};
/*!
    * @brief A class to store and provide the state of the solver while checkpointing, one SolutionVector object is supported.
    */
template<class SolutionVector, class GridVariables>
class SolverStateGridVar : public SolverStateBase
{
private:
    SolutionVector *x_;
    SolutionVector xCheckpoint_;
    GridVariables *gv_;

public:
    SolverStateGridVar(SolutionVector &x, GridVariables &gv)
        : x_(&x), xCheckpoint_(*x_), gv_(&gv)
    {
    }

    void writeState() override { xCheckpoint_ = *x_; }

    void readState() override
    {
        *x_ = xCheckpoint_;
        gv_->update(*x_);
    }
};
}  // namespace Dumux::Precice
#endif
