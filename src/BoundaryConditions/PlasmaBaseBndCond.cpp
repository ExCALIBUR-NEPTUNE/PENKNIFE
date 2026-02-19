#include "PlasmaBaseBndCond.hpp"
#include "../EquationSystems/PlasmaSystem.hpp"

using namespace std;

namespace PENKNIFE
{
PlasmaBaseBndCondFactory &GetPlasmaBaseBndCondFactory()
{
    static PlasmaBaseBndCondFactory instance;
    return instance;
}

PlasmaBaseBndCond::PlasmaBaseBndCond(
    const LU::SessionReaderSharedPtr &pSession,
    const std::weak_ptr<PlasmaSystem> &pSystem,
    const Array<OneD, MR::ExpListSharedPtr> &pFields,
    const Array<OneD, MR::DisContFieldSharedPtr> &pB,
    const Array<OneD, MR::DisContFieldSharedPtr> &pE,
    Array<OneD, SpatialDomains::BoundaryConditionShPtr> cond,
    Array<OneD, MultiRegions::ExpListSharedPtr> exp, const int pSpaceDim,
    const int bcRegion)
    : m_session(pSession), m_system(pSystem), m_fields(pFields),
      m_spacedim(pSpaceDim), m_bcRegion(bcRegion), B(pB), E(pE),
      m_varConv(pSystem.lock()->m_varConv),
      field_to_index(pSystem.lock()->field_to_index)
{
    for (int v = 0; v < m_fields.size(); ++v)
    {
        m_bndExp[v] = m_fields[v]->GetBndCondExpansions()[m_bcRegion];
    }
    m_nEdgePts    = m_bndExp[0]->GetTotPoints();
    m_nEdgeCoeffs = m_bndExp[0]->GetNcoeffs();

    m_fields[0]->GetBndElmtExpansion(m_bcRegion, m_bndElmtExp, false);

    m_fields[0]->GetBoundaryNormals(m_bcRegion, m_normals);

    B_bnd = Array<OneD, Array<OneD, NekDouble>>(3);
    E_bnd = Array<OneD, Array<OneD, NekDouble>>(3);
    for (int d = 0; d < 3; d++)
    {
        B_bnd[d] = Array<OneD, NekDouble>(m_nEdgePts);
        E_bnd[d] = Array<OneD, NekDouble>(m_nEdgePts);
    }

    m_diffusionAveWeight = 1.0;
}

/**
 * @param   bcRegion      id of the boundary region
 * @param   cnt
 * @param   Fwd
 * @param   physarray
 * @param   time
 */
void PlasmaBaseBndCond::Apply(
    const Array<OneD, const Array<OneD, NekDouble>> &physarray,
    const NekDouble &time)
{
    // Get EM fields on boundary
    Vmath::Zero(m_nEdgePts, mag_B, 1);
    for (int d = 0; d < m_spacedim; ++d)
    {
        this->B[d]->ExtractPhysToBnd(m_bcRegion, this->B[d]->GetPhys(),
                                     this->B_bnd[d]);
        this->E[d]->ExtractPhysToBnd(m_bcRegion, this->E[d]->GetPhys(),
                                     this->E_bnd[d]);
        Vmath::Vvtvp(m_nEdgePts, B_bnd[d], 1, B_bnd[d], 1, mag_B, 1, mag_B, 1);
    }
    Array<OneD, Array<OneD, NekDouble>> Fwd(physarray.size());

    for (int i = 0; i < physarray.size(); ++i)
    {
        m_fields[i]->ExtractPhysToBnd(m_bcRegion, physarray[i], Fwd[i]);
    }

    v_Apply(Fwd, physarray, time);
}

/**
 * @ brief Newly added bc should specify this virtual function
 * if the Bwd/value in m_bndCondExpansions is the target value like Direchlet
 * bc weight should be 1.0.
 * if some average Fwd and Bwd/value in m_bndCondExpansions
 * is the target value like WallViscousBC weight should be 0.5.
 */
void PlasmaBaseBndCond::v_ApplyBwdWeight()
{
    m_fields[0]->SetBndCondBwdWeight(m_bcRegion, m_diffusionAveWeight);
}

} // namespace PENKNIFE
