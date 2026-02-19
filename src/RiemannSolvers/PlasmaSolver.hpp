#ifndef PLASMASOLVER_HPP
#define PLASMASOLVER_HPP

#include <SolverUtils/RiemannSolvers/RiemannSolver.h>
#include <nektar_interface/solver_base/neso_reader.hpp>

using namespace Nektar;
namespace LU = Nektar::LibUtilities;
namespace SU = Nektar::SolverUtils;
namespace PENKNIFE
{
class PlasmaSolver : public SU::RiemannSolver
{
public:
    int omega_idx;

protected:
    bool m_pointSolve;

    PlasmaSolver(const LU::SessionReaderSharedPtr &pSession);

    using ND = NekDouble;
    void v_Solve(const int nDim, const Array<OneD, const Array<OneD, ND>> &Fwd,
                 const Array<OneD, const Array<OneD, ND>> &Bwd,
                 Array<OneD, Array<OneD, ND>> &flux) override;

    virtual void v_ArraySolve(
        [[maybe_unused]] const Array<OneD, const Array<OneD, ND>> &Fwd,
        [[maybe_unused]] const Array<OneD, const Array<OneD, ND>> &Bwd,
        [[maybe_unused]] Array<OneD, Array<OneD, ND>> &flux)
    {
        NEKERROR(ErrorUtil::efatal,
                 "This function should be defined by subclasses.");
    }

    virtual void v_PointSolve(
        [[maybe_unused]] ND rhoL, [[maybe_unused]] ND rhouL,
        [[maybe_unused]] ND rhovL, [[maybe_unused]] ND rhowL,
        [[maybe_unused]] ND EL, [[maybe_unused]] ND rhoR,
        [[maybe_unused]] ND rhouR, [[maybe_unused]] ND rhovR,
        [[maybe_unused]] ND rhowR, [[maybe_unused]] ND ER,
        [[maybe_unused]] ND &rhof, [[maybe_unused]] ND &rhouf,
        [[maybe_unused]] ND &rhovf, [[maybe_unused]] ND &rhowf,
        [[maybe_unused]] ND &Ef)
    {
        NEKERROR(ErrorUtil::efatal,
                 "This function should be defined by subclasses.");
    }
};
} // namespace PENKNIFE

#endif
