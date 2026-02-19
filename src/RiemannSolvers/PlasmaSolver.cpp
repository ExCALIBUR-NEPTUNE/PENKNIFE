#include "PlasmaSolver.hpp"

namespace PENKNIFE
{

/**
 * @brief Default constructor.
 */
PlasmaSolver::PlasmaSolver(
    const LU::SessionReaderSharedPtr &pSession)
    : RiemannSolver(pSession)
{
}

void PlasmaSolver::v_Solve(
    const int nDim, const Array<OneD, const Array<OneD, NekDouble>> &Fwd,
    const Array<OneD, const Array<OneD, NekDouble>> &Bwd,
    Array<OneD, Array<OneD, NekDouble>> &flux)
{
    if (m_pointSolve)
    {
        size_t expDim = nDim;
        NekDouble rhouf{}, rhovf{};

        if (expDim == 1)
        {
            for (size_t i = 0; i < Fwd[0].size(); ++i)
            {
                v_PointSolve(Fwd[0][i], Fwd[1][i], 0.0, 0.0, Fwd[2][i],
                             Bwd[0][i], Bwd[1][i], 0.0, 0.0, Bwd[2][i],
                             flux[0][i], flux[1][i], rhouf, rhovf, flux[2][i]);
            }
        }
        else if (expDim == 2)
        {
            for (size_t i = 0; i < Fwd[0].size(); ++i)
            {
                v_PointSolve(Fwd[0][i], Fwd[1][i], Fwd[2][i], 0.0, Fwd[3][i],
                             Bwd[0][i], Bwd[1][i], Bwd[2][i], 0.0, Bwd[3][i],
                             flux[0][i], flux[1][i], flux[2][i], rhovf,
                             flux[3][i]);
            }
        }
        else if (expDim == 3)
        {
            for (size_t i = 0; i < Fwd[0].size(); ++i)
            {
                v_PointSolve(Fwd[0][i], Fwd[1][i], Fwd[2][i], Fwd[3][i],
                             Fwd[4][i], Bwd[0][i], Bwd[1][i], Bwd[2][i],
                             Bwd[3][i], Bwd[4][i], flux[0][i], flux[1][i],
                             flux[2][i], flux[3][i], flux[4][i]);
            }
        }
    }
    else
    {
        v_ArraySolve(Fwd, Bwd, flux);
    }
}

} // namespace PENKNIFE
