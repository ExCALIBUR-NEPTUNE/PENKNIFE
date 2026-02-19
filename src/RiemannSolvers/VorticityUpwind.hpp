#ifndef VORTICITY_UPWINDSOLVER_HPP
#define VORTICITY_UPWINDSOLVER_HPP

#include "PlasmaSolver.hpp"

using namespace Nektar;
namespace LU = Nektar::LibUtilities;
namespace SU = Nektar::SolverUtils;
namespace PENKNIFE
{
class VorticityUpwindSolver : public PlasmaSolver
{
public:
    static SU::RiemannSolverSharedPtr create(
        const LU::SessionReaderSharedPtr &pSession)
    {
        return SU::RiemannSolverSharedPtr(new VorticityUpwindSolver(pSession));
    }

    static std::string solverName;

protected:
    VorticityUpwindSolver(const LU::SessionReaderSharedPtr &pSession);

    void v_ArraySolve(const Array<OneD, const Array<OneD, NekDouble>> &Fwd,
                      const Array<OneD, const Array<OneD, NekDouble>> &Bwd,
                      Array<OneD, Array<OneD, NekDouble>> &flux) final;
};
} // namespace PENKNIFE

#endif
