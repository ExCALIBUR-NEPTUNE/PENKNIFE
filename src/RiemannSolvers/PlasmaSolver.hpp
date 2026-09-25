#ifndef PLASMASOLVER_HPP
#define PLASMASOLVER_HPP

#include <SolverUtils/RiemannSolvers/RiemannSolver.h>
#include <nektar_interface/solver_base/neso_reader.hpp>

namespace PENKNIFE
{
using namespace Nektar;
namespace LU = Nektar::LibUtilities;
namespace SU = Nektar::SolverUtils;

class PlasmaSolver : public SU::RiemannSolver
{
public:
    int omega_idx;
    std::weak_ptr<PlasmaSystem> m_system;

protected:
    bool m_pointSolve;

    PlasmaSolver(const LU::SessionReaderSharedPtr &pSession);

    void v_Solve(const int nDim,
                 const Array<OneD, const Array<OneD, NekDouble>> &Fwd,
                 const Array<OneD, const Array<OneD, NekDouble>> &Bwd,
                 Array<OneD, Array<OneD, NekDouble>> &flux) override;

    virtual void v_ArraySolve(
        [[maybe_unused]] const Array<OneD, const Array<OneD, NekDouble>> &Fwd,
        [[maybe_unused]] const Array<OneD, const Array<OneD, NekDouble>> &Bwd,
        [[maybe_unused]] Array<OneD, Array<OneD, NekDouble>> &flux)
    {
        NEKERROR(ErrorUtil::efatal,
                 "This function should be defined by subclasses.");
    }

    virtual void v_PointSolve(
        [[maybe_unused]] NekDouble rhoL, [[maybe_unused]] NekDouble rhouL,
        [[maybe_unused]] NekDouble rhovL, [[maybe_unused]] NekDouble rhowL,
        [[maybe_unused]] NekDouble EL, [[maybe_unused]] NekDouble rhoR,
        [[maybe_unused]] NekDouble rhouR, [[maybe_unused]] NekDouble rhovR,
        [[maybe_unused]] NekDouble rhowR, [[maybe_unused]] NekDouble ER,
        [[maybe_unused]] NekDouble &rhof, [[maybe_unused]] NekDouble &rhouf,
        [[maybe_unused]] NekDouble &rhovf, [[maybe_unused]] NekDouble &rhowf,
        [[maybe_unused]] NekDouble &Ef)
    {
        NEKERROR(ErrorUtil::efatal,
                 "This function should be defined by subclasses.");
    }
};
} // namespace PENKNIFE

#endif
