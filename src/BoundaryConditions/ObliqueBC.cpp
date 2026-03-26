#include "ObliqueBC.hpp"
#include "../EquationSystems/PlasmaSystem.hpp"

using namespace std;

namespace PENKNIFE
{

std::string ObliqueBC::className =
    GetPlasmaBaseBndCondFactory().RegisterCreatorFunction(
        "Oblique", ObliqueBC::create, "Oblique boundary condition.");

ObliqueBC::ObliqueBC(const LU::SessionReaderSharedPtr &pSession,
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
}

void ObliqueBC::v_Apply(
    const Array<OneD, const Array<OneD, NekDouble>> &Fwd,
    const Array<OneD, const Array<OneD, NekDouble>> &physarray,
    [[maybe_unused]] const NekDouble &time)
{

    for (const auto &[k, v] : m_system.lock()->GetIons())
    {
        int ni_idx = v.fields.at(field_to_index.at("n"));

        Array<OneD, NekDouble> ni_bndelmt;
        m_fields[ni_idx]->ExtractPhysToBndElmt(m_bcRegion, physarray[ni_idx],
                                               ni_bndelmt);
        Array<OneD, Array<OneD, NekDouble>> ni_grad(m_spacedim);
        for (int d = 0; d < m_spacedim; d++)
        {
            ni_grad[d] = Array<OneD, NekDouble>(ni_bndelmt.size(), 0.0);
        }
        if (m_spacedim == 2)
        {
            m_bndElmtExp->PhysDeriv(ni_bndelmt, ni_grad[0], ni_grad[1]);
        }
        else if (m_spacedim == 3)
        {
            m_bndElmtExp->PhysDeriv(ni_bndelmt, ni_grad[0], ni_grad[1],
                                    ni_grad[2]);
        }
        Array<OneD, Array<OneD, NekDouble>> ni_bndGrad(m_spacedim);

        for (int d = 0; d < m_spacedim; d++)
        {
            ni_bndGrad[d] = Array<OneD, NekDouble>(m_nEdgePts, 0.0);
            m_fields[ni_idx]->ExtractElmtToBndPhys(m_bcRegion, ni_grad[d],
                                                      ni_bndGrad[d]);
        }

        Array<OneD, NekDouble> ni_par(m_nEdgePts, 0.0);
        for (int d = 0; d < m_spacedim; ++d)
        {
            Vmath::Vvtvp(m_nEdgePts, ni_bndGrad[d], 1, B_bnd[d], 1, ni_par, 1,
                         ni_par, 1);
            Vmath::Vdiv(m_nEdgePts, ni_par, 1, mag_B, 1, ni_par, 1);
        }
        Array<OneD, Array<OneD, NekDouble>> ni_perp(m_spacedim);
        Array<OneD, NekDouble> result(m_nEdgePts, 0.0);

        for (int d = 0; d < m_spacedim; ++d)
        {
            Vmath::Vvtvm(m_nEdgePts, ni_par, 1, B_bnd[d], 1, ni_bndGrad[d], 1,
                         ni_perp[d], 1);
            Vmath::Vvtvp(m_nEdgePts, ni_perp[d], 1, m_normals[d], 1, result, 1,
                         result, 1);
        }

        LibUtilities::EquationSharedPtr cond =
            std::static_pointer_cast<SpatialDomains::NeumannBoundaryCondition>(
                m_fields[ni_idx]->GetBndConditions()[m_bcRegion])
                ->m_neumannCondition;

        Array<OneD, NekDouble> x0(m_nEdgePts, 0.0);
        Array<OneD, NekDouble> x1(m_nEdgePts, 0.0);
        Array<OneD, NekDouble> x2(m_nEdgePts, 0.0);

        m_bndExp[ni_idx]->GetCoords(x0, x1, x2);

        Array<OneD, NekDouble> a(m_nEdgePts, 0.0);
        cond->Evaluate(x0, x1, x2, time, a);

        Vmath::Vsub(m_nEdgePts, a, 1, result, 1, result, 1);

        m_bndExp[ni_idx]->IProductWRTBase(result,
                                          m_bndExp[ni_idx]->UpdateCoeffs());
    }
}

} // namespace PENKNIFE
