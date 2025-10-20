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
    virtual void readState(double dt) = 0;
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

    void readState(double dt) override { *x_ = xCheckpoint_; }
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

    void readState(double dt) override
    {
        *x_ = xCheckpoint_;
        gv_->update(*x_);
    }
};
/*!
    * @brief A class to store and provide the state of the solver while checkpointing, one SolutionVector object is supported.
    */
template<class SolutionVector, class TimeLoop, class GridVariables>
class SolverStateGridVarTime : public SolverStateBase
{
private:
    SolutionVector *x_;
    SolutionVector xCheckpoint_;
    TimeLoop *tl_;
    double timeCheckpoint_{0.0};
    long timeStepCheckpoint_{0};
    GridVariables *gv_;

public:
    SolverStateGridVarTime(SolutionVector &x, TimeLoop &tl, GridVariables &gv)
        : x_(&x), xCheckpoint_(*x_), tl_(&tl), gv_(&gv)
    {
    }

    void writeState() override
    {
        xCheckpoint_ = *x_;
        timeCheckpoint_ = tl_->time();
        timeStepCheckpoint_ = tl_->timeStepIndex();
    }

    void readState(double dt) override
    {
        *x_ = xCheckpoint_;
        tl_->setTime(timeCheckpoint_, timeStepCheckpoint_);
        tl_->setTimeStepSize(dt);
        gv_->update(*x_);
    }
};
}  // namespace Dumux::Precice
#endif
