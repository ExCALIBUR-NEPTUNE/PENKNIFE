#ifndef PLASMA_BASEBNDCOND_HPP
#define PLASMA_BASEBNDCOND_HPP

#include <string>

#include "../Misc/VariableConverter.hpp"
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>

using namespace Nektar;
namespace LU = Nektar::LibUtilities;
namespace MR = Nektar::MultiRegions;
namespace SD = Nektar::SpatialDomains;

namespace PENKNIFE
{
class PlasmaSystem;
class PlasmaBaseBndCond;

/// A shared pointer to a boundary condition object
typedef std::shared_ptr<PlasmaBaseBndCond> PlasmaBaseBndCondSharedPtr;

/// Declaration of the boundary condition factory
typedef LU::NekFactory<
    std::string, PlasmaBaseBndCond, const LU::SessionReaderSharedPtr &,
    const std::weak_ptr<PlasmaSystem> &,
    const Array<OneD, MR::ExpListSharedPtr> &,
    const Array<OneD, MR::DisContFieldSharedPtr> &,
    const Array<OneD, MR::DisContFieldSharedPtr> &,
    Array<OneD, SpatialDomains::BoundaryConditionShPtr>,
    Array<OneD, MultiRegions::ExpListSharedPtr>, const int, const int>
    PlasmaBaseBndCondFactory;

/// Declaration of the boundary condition factory singleton
PlasmaBaseBndCondFactory &GetPlasmaBaseBndCondFactory();

/**
 * @class PlasmaBaseBndCond
 * @brief Encapsulates the user-defined boundary conditions for the
 *        plasma solver.
 */
class PlasmaBaseBndCond
{
public:
    virtual ~PlasmaBaseBndCond()
    {
    }

    /// Apply the boundary condition
    void Apply(const Array<OneD, const Array<OneD, NekDouble>> &physarray,
               const NekDouble &time = 0);

    /// Apply the Weight of boundary condition
    void ApplyBwdWeight()
    {
        v_ApplyBwdWeight();
    }

    int ee_idx;
    int omega_idx;
    int phi_idx;

protected:
    /// Session reader
    LU::SessionReaderSharedPtr m_session;

    const std::weak_ptr<PlasmaSystem> m_system;

    std::map<int, SpatialDomains::BoundaryConditionShPtr> m_bndConds;
    std::map<int, MultiRegions::ExpListSharedPtr> m_bndExp;
    /// Expansion of boundary adjacent elements
    MultiRegions::ExpListSharedPtr m_bndElmtExp;

    /// Array of fields
    Array<OneD, MR::ExpListSharedPtr> m_fields;
    /// Trace normals
    Array<OneD, Array<OneD, NekDouble>> m_normals;

    // EM fields
    Array<OneD, MR::DisContFieldSharedPtr> B;
    Array<OneD, MR::DisContFieldSharedPtr> E;

    /// Trace EM fields
    Array<OneD, Array<OneD, NekDouble>> B_bnd;
    Array<OneD, Array<OneD, NekDouble>> E_bnd;

    Array<OneD, NekDouble> mag_B;

    /// Space dimension
    int m_spacedim;
    /// Auxiliary object to convert variables
    VariableConverterSharedPtr m_varConv;
    const NektarFieldIndexMap &field_to_index;

    /// Weight for average calculation of diffusion term
    NekDouble m_diffusionAveWeight;

    int m_variable;
    /// Id of the boundary region
    int m_bcRegion;

    int m_nEdgePts;
    int m_nEdgeCoeffs;

    /// Constructor
    PlasmaBaseBndCond(const LU::SessionReaderSharedPtr &pSession,
                       const std::weak_ptr<PlasmaSystem> &pSystem,
                       const Array<OneD, MR::ExpListSharedPtr> &pFields,
                       const Array<OneD, MR::DisContFieldSharedPtr> &pB,
                       const Array<OneD, MR::DisContFieldSharedPtr> &pE,
                       Array<OneD, SpatialDomains::BoundaryConditionShPtr> cond,
                       Array<OneD, MultiRegions::ExpListSharedPtr> exp,
                       const int pSpaceDim, const int bcRegion);

    virtual void v_Apply(
        const Array<OneD, const Array<OneD, NekDouble>> &Fwd,
        const Array<OneD, const Array<OneD, NekDouble>> &physarray,
        const NekDouble &time) = 0;

    virtual void v_ApplyBwdWeight();
};
} // namespace PENKNIFE
#endif