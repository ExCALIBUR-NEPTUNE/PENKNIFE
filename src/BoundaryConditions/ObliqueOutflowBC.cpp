#include "ObliqueOutflowBC.hpp"
#include "../EquationSystems/PlasmaSystem.hpp"

using namespace std;

namespace PENKNIFE
{

std::string ObliqueOutflowBC::className =
    GetPlasmaBaseBndCondFactory().RegisterCreatorFunction(
        "ObliqueOutflow", ObliqueOutflowBC::create,
        "ObliqueOutflow boundary condition.");

ObliqueOutflowBC::ObliqueOutflowBC(
    const LU::SessionReaderSharedPtr &pSession,
    const std::weak_ptr<PlasmaSystem> &pSystem,
    const Array<OneD, MR::ExpListSharedPtr> &pFields,
    const Array<OneD, MR::DisContFieldSharedPtr> &pB,
    const Array<OneD, MR::DisContFieldSharedPtr> &pE,
    Array<OneD, SpatialDomains::BoundaryConditionShPtr> cond,
    Array<OneD, MultiRegions::ExpListSharedPtr> exp, const int pSpaceDim,
    const int bcRegion)
    : PlasmaBaseBndCond(pSession, pSystem, pFields, pB, pE, cond, exp,
                         pSpaceDim, bcRegion)
{
    m_session->LoadParameter("gamma", gamma, 7.0);
    m_session->LoadParameter("m_i", m_i);
    m_session->LoadParameter("k_B", k_B);
    m_session->LoadParameter("k_par", this->k_par, 100.0);
    m_session->LoadParameter("k_perp", this->k_perp, 1.0);
    m_session->LoadParameter("kappa_par", this->kappa_par, 100.0);
    m_session->LoadParameter("kappa_perp", this->kappa_perp, 1.0);

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            m_D[i][j]     = Array<OneD, NekDouble>(m_nEdgePts, 0.0);
            m_kappa[i][j] = Array<OneD, NekDouble>(m_nEdgePts, 0.0);
        }
    }
}

void ObliqueOutflowBC::v_Apply(
    const Array<OneD, const Array<OneD, NekDouble>> &Fwd,
    const Array<OneD, const Array<OneD, NekDouble>> &physarray,
    [[maybe_unused]] const NekDouble &time)
{
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> bndGrad(m_fields.size());

    for (int i = 0; i < m_fields.size(); ++i)
    {
        // Obtain field on boundary elmts
        Array<OneD, NekDouble> elmt;
        m_fields[i]->ExtractPhysToBndElmt(m_bcRegion, physarray[i], elmt);

        // Obtain field derivative
        Array<OneD, Array<OneD, NekDouble>> grad(m_spacedim);
        for (int k = 0; k < m_spacedim; k++)
        {
            grad[k] = Array<OneD, NekDouble>(elmt.size(), 0.0);
        }
        if (m_spacedim == 2)
        {
            m_bndElmtExp->PhysDeriv(elmt, grad[0], grad[1]);
        }
        else if (m_spacedim == 3)
        {
            m_bndElmtExp->PhysDeriv(elmt, grad[0], grad[1], grad[2]);
        }
        // Obtain gradient on boundary of elements
        bndGrad[i] = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
        for (int k = 0; k < m_spacedim; k++)
        {
            bndGrad[i][k] = Array<OneD, NekDouble>(m_nEdgePts, 0.0);
            m_fields[i]->ExtractElmtToBndPhys(m_bcRegion, grad[k],
                                              bndGrad[i][k]);
        }
    }
    // Obtain anisotropic tensors
    CalcDTensor();
    CalcKappaTensor();
    // Obtain Fluxes
    Array<OneD, Array<OneD, NekDouble>> density_flux(m_spacedim);
    Array<OneD, Array<OneD, NekDouble>> energy_flux(m_spacedim);

    // Obtain rhs of bcs
    Array<OneD, NekDouble> density_bc(m_nEdgePts, 0.0);
    Array<OneD, NekDouble> energy_bc(m_nEdgePts, 0.0);

    for (int p = 0; p < m_nEdgePts; ++p)
    {
        NekDouble cs  = -std::sqrt(k_B * m_bndExp[2]->GetPhys()[p] / m_i);
        density_bc[p] = m_bndExp[0]->GetPhys()[p] * cs;
        energy_bc[p]  = gamma * k_B * density_bc[p] * m_bndExp[2]->GetPhys()[p];
    }

    for (int k = 0; k < m_spacedim; k++)
    {
        density_flux[k] = Array<OneD, NekDouble>(m_nEdgePts, 0.0);
        energy_flux[k]  = Array<OneD, NekDouble>(m_nEdgePts, 0.0);

        for (int j = 0; j < m_spacedim; j++)
        {
            Vmath::Vvtvp(m_nEdgePts, m_D[k][j], 1, bndGrad[0][j], 1,
                         density_flux[k], 1, density_flux[k], 1);
            Vmath::Vvtvp(m_nEdgePts, m_kappa[k][j], 1, bndGrad[1][j], 1,
                         energy_flux[k], 1, energy_flux[k], 1);
        }
        Vmath::Vvtvp(m_nEdgePts, density_flux[k], 1, B_bnd[k], 1, density_bc, 1,
                     density_bc, 1);
        Vmath::Vvtvp(m_nEdgePts, energy_flux[k], 1, B_bnd[k], 1, energy_bc, 1,
                     energy_bc, 1);
    }

    m_bndExp[0]->FwdTrans(density_bc, m_bndExp[0]->UpdateCoeffs());
    m_bndExp[1]->FwdTrans(energy_bc, m_bndExp[1]->UpdateCoeffs());
}

void ObliqueOutflowBC::CalcKPar()
{
    NekDouble k_par = this->k_par;
    kpar            = Array<OneD, NekDouble>(m_nEdgePts, k_par);
}

void ObliqueOutflowBC::CalcKPerp()
{
    NekDouble k_perp = this->k_perp;
    m_session->LoadParameter("k_perp", k_perp, 1.0);
    kperp = Array<OneD, NekDouble>(m_nEdgePts, k_perp);
}

void ObliqueOutflowBC::CalcDTensor()
{
    CalcKPar();
    CalcKPerp();
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            for (int k = 0; k < m_nEdgePts; k++)
            {
                m_D[i][j][k] =
                    (kpar[k] - kperp[k]) * B_bnd[i][k] * B_bnd[j][k] / mag_B[k];
                if (i == j)
                {
                    m_D[i][j][k] += kperp[k];
                }
            }
        }
    }
}

void ObliqueOutflowBC::CalcKappaPar()
{
    NekDouble kappa_par = this->kappa_par;
    kpar                = Array<OneD, NekDouble>(m_nEdgePts, kappa_par);
}

void ObliqueOutflowBC::CalcKappaPerp()
{
    NekDouble kappa_perp = this->kappa_perp;
    kperp                = Array<OneD, NekDouble>(m_nEdgePts, kappa_perp);
}

void ObliqueOutflowBC::CalcKappaTensor()
{
    CalcKappaPar();
    CalcKappaPerp();
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            for (int k = 0; k < m_nEdgePts; k++)
            {
                m_kappa[i][j][k] = (kappapar[k] - kappaperp[k]) * B_bnd[i][k] *
                                   B_bnd[j][k] / mag_B[k];
                if (i == j)
                {
                    m_kappa[i][j][k] += kappaperp[k];
                }
            }
        }
    }
}

} // namespace PENKNIFE
