#ifndef SOLVERSTATE_HH
#define SOLVERSTATE_HH

#include <precice/precice.hpp>

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
/*!
 * @brief A class to store and provide the state of the solver while checkpointing, one SolutionVector object is supported.
 */
template<class SolutionVector, class GridVariables, class TimeLoop>
class SolverStateGridVarTimeLoop : public SolverStateBase
{
private:
    SolutionVector *x_;
    SolutionVector xCheckpoint_;
    GridVariables *gv_;
    TimeLoop *tl_;
    double timeCheckpoint_ = 0.0;
    double dtCheckpoint_ = 0.0;
    long timeStepCheckpoint_ = 0;

public:
    SolverStateGridVarTimeLoop(SolutionVector &x,
                               GridVariables &gv,
                               TimeLoop &tl)
        : x_(&x), xCheckpoint_(*x_), gv_(&gv), tl_(&tl)
    {
    }

    void writeState() override
    {
        xCheckpoint_ = *x_;
        timeCheckpoint_ = tl_->time();
        timeStepCheckpoint_ = tl_->timeStepIndex();
        dtCheckpoint_ = tl_->timeStepSize();
    }

    void readState() override
    {
        *x_ = xCheckpoint_;
        gv_->update(*x_);
        tl_->setTime(timeCheckpoint_, timeStepCheckpoint_);
        tl_->setTimeStepSize(dtCheckpoint_);
    }
};
}  // namespace Dumux::Precice
#endif
