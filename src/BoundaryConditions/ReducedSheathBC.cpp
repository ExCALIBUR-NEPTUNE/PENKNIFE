#include "ReducedSheathBC.hpp"
#include "../EquationSystems/PlasmaSystem.hpp"

using namespace std;

namespace PENKNIFE
{

std::string ReducedSheathBC::className =
    GetPlasmaBaseBndCondFactory().RegisterCreatorFunction(
        "ReducedSheath", ReducedSheathBC::create,
        "Reduced Sheath boundary condition.");

ReducedSheathBC::ReducedSheathBC(
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
    className = "Sheath";
    for (size_t i = 0; i < m_spacedim; ++i)
    {
        m_bndConds[i] = cond[i];
        if (cond[i]->GetUserDefined() == className)
        {
            m_bndExp[i] = exp[i];
        }
    }

    m_session->LoadParameter("Ge", Ge, 1.0);
    m_session->LoadParameter("wall", wall, 0.0);
    // this->v_ExB = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
    // for (int d = 0; d < m_spacedim; ++d)
    // {
    //     v_ExB[d] = Array<OneD, NekDouble>(m_nEdgePts);
    // }
    this->mag_B = Array<OneD, NekDouble>(m_nEdgePts, 0.0);
}

void ReducedSheathBC::v_Apply(
    const Array<OneD, const Array<OneD, NekDouble>> &Fwd,
    const Array<OneD, const Array<OneD, NekDouble>> &physarray,
    [[maybe_unused]] const NekDouble &time)
{

    // Get ExB velocity on boundary
    // Vmath::Vvtvvtm(m_nEdgePts, this->E_bnd[1], 1, this->B_bnd[2], 1,
    //                this->E_bnd[2], 1, this->B_bnd[1], 1, this->v_ExB[0], 1);
    // Vmath::Vdiv(m_nEdgePts, this->v_ExB[0], 1, this->mag_B, 1,
    // this->v_ExB[0],
    //             1);
    // Vmath::Vvtvvtm(m_nEdgePts, this->E_bnd[2], 1, this->B_bnd[0], 1,
    // E_bnd[0],
    //                1, B_bnd[2], 1, this->v_ExB[1], 1);
    // Vmath::Vdiv(m_nEdgePts, this->v_ExB[1], 1, this->mag_B, 1,
    // this->v_ExB[1],
    //             1);

    // if (m_spacedim == 3)
    // {
    //     Vmath::Vvtvvtm(m_nEdgePts, this->E_bnd[0], 1, this->B_bnd[1], 1,
    //                    this->E_bnd[1], 1, this->B_bnd[0], 1, this->v_ExB[2],
    //                    1);
    //     Vmath::Vdiv(m_nEdgePts, this->v_ExB[2], 1, this->mag_B, 1,
    //                 this->v_ExB[2], 1);
    // }

    Array<OneD, NekDouble> bn(m_nEdgePts, 0.0);
    for (int d = 0; d < m_spacedim; ++d)
    {
        Vmath::Vvtvp(m_nEdgePts, B_bnd[d], 1, m_normals[d], 1, bn, 1, bn, 1);
    }
    for (int p = 0; p < m_nEdgePts; ++p)
    {
        bn[p] = mag_B[p] == 0 ? 0 : bn[p] / std::sqrt(mag_B[p]);
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

        // Loop over points to set parallel ion momentum boundary condition
        Array<OneD, NekDouble> vi_bc(m_nEdgePts);
        Array<OneD, NekDouble> pi_bc(m_nEdgePts);
        Array<OneD, NekDouble> ni_bc(m_nEdgePts);

        for (int p = 0; p < m_nEdgePts; ++p)
        {
            NekDouble c_i_sq = (gamma_i * Fwd[pi_idx][p] / Fwd[ni_idx][p] +
                                v.charge * Fwd[pe_idx][p] / ne[p]) /
                               v.mass;
            NekDouble c_i = bn[p] > 0 ? std::sqrt(c_i_sq) : -std::sqrt(c_i_sq);
            NekDouble v_trial =
                Fwd[vi_idx][p] * bn[p] / (v.mass * Fwd[ni_idx][p]);
            NekDouble v_sheath = std::max(v_trial, c_i * bn[p]);

            vi_bc[p] = bn[p] == 0 ? Fwd[vi_idx][p]
                                  : v.mass * Fwd[ni_idx][p] * v_sheath / bn[p];
            pi_bc[p] = ((2 / 3) * gamma_i * Fwd[pi_idx][p] +
                        0.5 * vi_bc[p] * vi_bc[p] / (v.mass * Fwd[ni_idx][p])) *
                       v_sheath;
            ni_bc[p] = Fwd[ni_idx][p] * v_sheath;
        }

        Vmath::Vadd(m_nEdgePts, ni_bc, 1, ion_sum, 1, ion_sum, 1);

        // Momentum
        m_bndExp[vi_idx]->FwdTransBndConstrained(
            vi_bc, m_bndExp[vi_idx]->UpdateCoeffs());
        // Density
        m_bndExp[ni_idx]->IProductWRTBase(ni_bc,
                                          m_bndExp[ni_idx]->UpdateCoeffs());
        // Ion Pressure
        m_bndExp[pi_idx]->IProductWRTBase(pi_bc,
                                          m_bndExp[pi_idx]->UpdateCoeffs());
    }

    Array<OneD, NekDouble> phi_bc(m_nEdgePts, 0.0);
    Array<OneD, NekDouble> pe_bc(m_nEdgePts, 0.0);
    for (int p = 0; p < m_nEdgePts; ++p)
    {
        NekDouble Te = Fwd[pe_idx][p] / ne[p];
        phi_bc[p]    = wall + Te * std::log(std::sqrt(Te / (me * 2 * M_PI)) *
                                            (1. - Ge) * ne[p] / ion_sum[p]);
        NekDouble ve = -std::sqrt(Te / (me * 2 * M_PI)) * (1. - Ge) *
                       std::exp(-(phi_bc[p] - wall) / Te);

        pe_bc[p] = (2 / 3) * gamma_e * Fwd[pe_idx][p] * ve * bn[p];
        pe_bc[p] = 0.0;
    }

    // Electron Pressure
    m_bndExp[pe_idx]->IProductWRTBase(pe_bc, m_bndExp[pe_idx]->UpdateCoeffs());
    // Potential
    m_bndExp[phi_idx]->FwdTransBndConstrained(
        phi_bc, m_bndExp[phi_idx]->UpdateCoeffs());
}

} // namespace PENKNIFE
