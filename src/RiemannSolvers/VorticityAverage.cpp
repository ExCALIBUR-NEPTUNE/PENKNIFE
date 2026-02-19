#include "VorticityAverage.hpp"

namespace PENKNIFE
{
std::string VorticityAverageSolver::solverName =
    SU::GetRiemannSolverFactory().RegisterCreatorFunction(
        "VorticityAverage", VorticityAverageSolver::create, "Average solver");

/**
 * @class AverageSolver
 *
 * @brief Average scheme Riemann solver.
 *
 * The Average solver determines the flux based upon an advection
 * velocity \f$\mathbf{V}\f$ and trace normal \f$\mathbf{n}\f$. In
 * particular, the flux for each component of the velocity field is
 * deterined by:
 *
 * \f[ \mathbf{f}(u^+,u^-) = \begin{cases} \mathbf{V}u^+, &
 * \mathbf{V}\cdot\mathbf{n} \geq 0,\\ \mathbf{V}u^-, &
 * \mathbf{V}\cdot\mathbf{n} < 0.\end{cases} \f]
 *
 * Here the superscript + and - denotes forwards and backwards spaces
 * respectively.
 */

/**
 * @brief Default constructor.
 */
VorticityAverageSolver::VorticityAverageSolver(
    const LU::SessionReaderSharedPtr &pSession)
    : PlasmaSolver(pSession)
{
    m_pointSolve = false;
}

/**
 * @brief Implementation of the Average solver.
 *
 * The Average solver assumes that a scalar field Vn is defined, which
 * corresponds with the dot product \f$\mathbf{V}\cdot\mathbf{n}\f$,
 * where \f$\mathbf{V}\f$ is the advection velocity and \f$\mathbf{n}\f$
 * defines the normal of a vertex, edge or face at each quadrature point
 * of the trace space.
 *
 * @param Fwd   Forwards trace space.
 * @param Bwd   Backwards trace space.
 * @param flux  Resulting flux.
 */
void VorticityAverageSolver::v_ArraySolve(
    const Array<OneD, const Array<OneD, NekDouble>> &Fwd,
    const Array<OneD, const Array<OneD, NekDouble>> &Bwd,
    Array<OneD, Array<OneD, NekDouble>> &flux)
{
    ASSERTL1(CheckVectors("Vn"), "Vn not defined.");
    const Array<OneD, Array<OneD, NekDouble>> &traceVel = m_vectors["Vn"]();
    const Array<OneD, NekDouble> &flux_omega            = m_scalars["wf"]();

    for (int p = 0; p < traceVel[0].size(); ++p)
    {
        for (int i = 0; i < traceVel.size(); ++i)
        {
            flux[i][p] = traceVel[i][p] * 0.5 * (Fwd[i][p] + Bwd[i][p]);
        }

        flux[omega_idx][p] = flux_omega[p];
    }
}
} // namespace PENKNIFE
