#ifndef PLASMA_BNDCONDS_HPP
#define PLASMA_BNDCONDS_HPP

#include "PlasmaBaseBndCond.hpp"

using namespace Nektar;
namespace LU = Nektar::LibUtilities;
namespace MR = Nektar::MultiRegions;
namespace SD = Nektar::SpatialDomains;

namespace PENKNIFE
{

class PlasmaSystem;
class PlasmaBoundaryConditions;
typedef std::shared_ptr<PlasmaBoundaryConditions>
    IncBoundaryConditionsSharedPtr;

class PlasmaBoundaryConditions
{
public:
    PlasmaBoundaryConditions();

    void Initialize(const LibUtilities::SessionReaderSharedPtr pSession,
                    const std::weak_ptr<PlasmaSystem> &pSystem,
                    Array<OneD, MultiRegions::ExpListSharedPtr> pFields,
                    const Array<OneD, MR::DisContFieldSharedPtr> &pB,
                    const Array<OneD, MR::DisContFieldSharedPtr> &pE,
                    int pSpacedim);

    void Update(const Array<OneD, const Array<OneD, NekDouble>> &physarray,
                NekDouble time);

    std::map<int, PlasmaBaseBndCondSharedPtr> &GetBounds()
    {
        return m_bounds;
    };

protected:
    std::map<int, PlasmaBaseBndCondSharedPtr> m_bounds;
    static std::set<std::string> m_BndType;
    int m_spacedim;
    int m_bnd_dim;
};

} // namespace PENKNIFE
#endif