#include "BohmBC.hpp"
#include "../EquationSystems/PlasmaSystem.hpp"

using namespace std;

namespace PENKNIFE
{

std::string BohmBC::className =
    GetPlasmaBaseBndCondFactory().RegisterCreatorFunction(
        "Bohm", BohmBC::create, "Bohm boundary condition.");

BohmBC::BohmBC(const LU::SessionReaderSharedPtr &pSession,
               const std::weak_ptr<PlasmaSystem> &pSystem,
               const Array<OneD, MR::ExpListSharedPtr> &pFields,
               const Array<OneD, MR::DisContFieldSharedPtr> &pB,
               const Array<OneD, MR::DisContFieldSharedPtr> &pE,
               Array<OneD, SpatialDomains::BoundaryConditionShPtr> cond,
               Array<OneD, MultiRegions::ExpListSharedPtr> exp,
               const int pSpaceDim, const int bcRegion)
    : PlasmaBaseBndCond(pSession, pSystem, pFields, pB, pE, cond, exp,
                         pSpaceDim, bcRegion)
{
    m_session->LoadParameter("Ge", Ge, 0.0);
    m_session->LoadParameter("lambda", lambda, 0.0);
    this->v_ExB = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
    for (int d = 0; d < m_spacedim; ++d)
    {
        v_ExB[d] = Array<OneD, NekDouble>(m_nEdgePts);
    }
    this->mag_B = Array<OneD, NekDouble>(m_nEdgePts, 0.0);
}

void BohmBC::v_Apply(const Array<OneD, const Array<OneD, NekDouble>> &Fwd,
                     const Array<OneD, const Array<OneD, NekDouble>> &physarray,
                     [[maybe_unused]] const NekDouble &time)
{
    // Get ExB velocity on boundary
    Vmath::Vvtvvtm(m_nEdgePts, this->E_bnd[1], 1, this->B_bnd[2], 1,
                   this->E_bnd[2], 1, this->B_bnd[1], 1, this->v_ExB[0], 1);
    Vmath::Vdiv(m_nEdgePts, this->v_ExB[0], 1, this->mag_B, 1, this->v_ExB[0],
                1);
    Vmath::Vvtvvtm(m_nEdgePts, this->E_bnd[2], 1, this->B_bnd[0], 1, E_bnd[0],
                   1, B_bnd[2], 1, this->v_ExB[1], 1);
    Vmath::Vdiv(m_nEdgePts, this->v_ExB[1], 1, this->mag_B, 1, this->v_ExB[1],
                1);

    if (m_spacedim == 3)
    {
        Vmath::Vvtvvtm(m_nEdgePts, this->E_bnd[0], 1, this->B_bnd[1], 1,
                       this->E_bnd[1], 1, this->B_bnd[0], 1, this->v_ExB[2], 1);
        Vmath::Vdiv(m_nEdgePts, this->v_ExB[2], 1, this->mag_B, 1,
                    this->v_ExB[2], 1);
    }
    // Get electron density on boundary
    Array<OneD, NekDouble> ne(m_nEdgePts, 0.0);

    for (const auto &[k, v] : m_system.lock()->GetIons())
    {
        int ni_idx = v.fields.at(field_to_index.at("n"));
        Vmath::Svtvp(m_nEdgePts, v.charge, Fwd[ni_idx], 1, ne, 1, ne, 1);
    }

    Array<OneD, NekDouble> ion_sum(m_nEdgePts, 0.0);

    for (const auto &[k, v] : m_system.lock()->GetIons())
    {
        int ni_idx = v.fields.at(field_to_index.at("n"));
        int vi_idx = v.fields.at(field_to_index.at("v"));
        int pi_idx = v.fields.at(field_to_index.at("p"));

        // Obtain pressure gradient in boundary elements
        Array<OneD, NekDouble> pi_bndelmt;
        m_fields[0]->ExtractPhysToBndElmt(m_bcRegion, physarray[pi_idx],
                                          pi_bndelmt);

        Array<OneD, Array<OneD, NekDouble>> pi_grad(m_spacedim);
        for (int d = 0; d < m_spacedim; d++)
        {
            pi_grad[d] = Array<OneD, NekDouble>(pi_bndelmt.size(), 0.0);
        }
        if (m_spacedim == 2)
        {
            m_bndElmtExp->PhysDeriv(pi_bndelmt, pi_grad[0], pi_grad[1]);
        }
        else if (m_spacedim == 3)
        {
            m_bndElmtExp->PhysDeriv(pi_bndelmt, pi_grad[0], pi_grad[1],
                                    pi_grad[2]);
        }
        // Obtain pressure gradient on boundary
        Array<OneD, Array<OneD, NekDouble>> pi_bndGrad(m_spacedim);
        Array<OneD, Array<OneD, NekDouble>> v_di(m_spacedim);

        for (int d = 0; d < m_spacedim; d++)
        {
            pi_bndGrad[d] = Array<OneD, NekDouble>(m_nEdgePts, 0.0);
            v_di[d]       = Array<OneD, NekDouble>(m_nEdgePts, 0.0);
            m_fields[0]->ExtractElmtToBndPhys(m_bcRegion, pi_grad[d],
                                              pi_bndGrad[d]);
        }

        Vmath::Vvtvvtm(m_nEdgePts, pi_bndGrad[1], 1, B_bnd[2], 1, pi_bndGrad[2],
                       1, B_bnd[1], 1, v_di[0], 1);
        Vmath::Vdiv(m_nEdgePts, v_di[0], 1, this->mag_B, 1, v_di[0], 1);
        Vmath::Vvtvvtm(m_nEdgePts, pi_bndGrad[2], 1, B_bnd[0], 1, pi_bndGrad[0],
                       1, B_bnd[2], 1, v_di[1], 1);
        Vmath::Vdiv(m_nEdgePts, v_di[1], 1, this->mag_B, 1, v_di[1], 1);

        if (m_spacedim == 3)
        {
            Vmath::Vvtvvtm(m_nEdgePts, pi_bndGrad[0], 1, B_bnd[1], 1,
                           pi_bndGrad[1], 1, B_bnd[0], 1, v_di[2], 1);
            Vmath::Vdiv(m_nEdgePts, v_di[2], 1, this->mag_B, 1, v_di[2], 1);
        }

        // Loop over points to set parallel ion momentum boundary condition
        Array<OneD, NekDouble> v_par(m_nEdgePts);

        for (int p = 0; p < m_nEdgePts; ++p)
        {
            NekDouble v_en = 0;
            NekDouble v_dn = 0;
            NekDouble bn   = 0;
            for (int d = 0; d < m_spacedim; ++d)
            {
                v_en += m_normals[d][p] * v_ExB[d][p];
                v_dn += m_normals[d][p] * v_di[d][p];
                bn += m_normals[d][p] * B_bnd[d][p];
            }
            bn /= std::sqrt(mag_B[p]);

            NekDouble c_i_sq = (gamma_i * Fwd[pi_idx][p] / Fwd[ni_idx][p] +
                                v.charge * Fwd[pe_idx][p] / ne[p]) /
                               v.mass;
            NekDouble c_i = bn > 0 ? std::sqrt(c_i_sq) : -std::sqrt(c_i_sq);

            v_par[p] = std::max(v.mass * Fwd[ni_idx][p] *
                                    (c_i * bn - v_en - v_dn) / bn,
                                Fwd[vi_idx][p]);
        }
        m_bndExp[vi_idx]->FwdTransBndConstrained(
            v_par, m_bndExp[vi_idx]->UpdateCoeffs());
    }

    Array<OneD, NekDouble> phi(m_nEdgePts, 0.0);
    for (int i = 0; i < m_nEdgePts; ++i)
    {
        phi[i] = m_bndExp[pe_idx]->GetPhys()[i] *
                 log(sqrt(m_bndExp[pe_idx]->GetPhys()[i] / (me * 2 * M_PI)) *
                     (1. - Ge) / ion_sum[i]);
    }
    Array<OneD, NekDouble> wall_potential(m_nEdgePts, lambda);
    Vmath::Vadd(m_nEdgePts, wall_potential, 1, phi, 1, phi, 1);

    m_bndExp[0]->FwdTrans(phi, m_bndExp[phi_idx]->UpdateCoeffs());
}

} // namespace PENKNIFE
