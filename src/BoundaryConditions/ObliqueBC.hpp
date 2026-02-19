#ifndef PLASMA_OBLIQUE_HPP
#define PLASMA_OBLIQUE_HPP

#include "PlasmaBaseBndCond.hpp"

namespace PENKNIFE
{

class ObliqueBC : public PlasmaBaseBndCond
{
public:
    friend class MemoryManager<ObliqueBC>;

    static PlasmaBaseBndCondSharedPtr create(
        const LU::SessionReaderSharedPtr &pSession,
        const std::weak_ptr<PlasmaSystem> &pSystem,
        const Array<OneD, MR::ExpListSharedPtr> &pFields,
        const Array<OneD, MR::DisContFieldSharedPtr> &pB,
        const Array<OneD, MR::DisContFieldSharedPtr> &pE,
        Array<OneD, SpatialDomains::BoundaryConditionShPtr> cond,
        Array<OneD, MultiRegions::ExpListSharedPtr> exp, const int pSpaceDim,
        const int bcRegion)
    {
        PlasmaBaseBndCondSharedPtr p =
            MemoryManager<ObliqueBC>::AllocateSharedPtr(
                pSession, pSystem, pFields, pB, pE, cond, exp, pSpaceDim,
                bcRegion);
        return p;
    }

    static std::string className;

protected:
    void v_Apply(const Array<OneD, const Array<OneD, NekDouble>> &Fwd,
                 const Array<OneD, const Array<OneD, NekDouble>> &physarray,
                 const NekDouble &time) override;

private:
    ObliqueBC(const LU::SessionReaderSharedPtr &pSession,
              const std::weak_ptr<PlasmaSystem> &pSystem,
              const Array<OneD, MR::ExpListSharedPtr> &pFields,
              const Array<OneD, MR::DisContFieldSharedPtr> &pB,
              const Array<OneD, MR::DisContFieldSharedPtr> &pE,
              Array<OneD, SpatialDomains::BoundaryConditionShPtr> cond,
              Array<OneD, MultiRegions::ExpListSharedPtr> exp,
              const int pSpaceDim, const int bcRegion);
    ~ObliqueBC() override {};
};

} // namespace PENKNIFE
#endif