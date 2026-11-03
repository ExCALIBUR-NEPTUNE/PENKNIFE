#ifndef PLASMA_WALL_HPP
#define PLASMA_WALL_HPP

#include "PlasmaBaseBndCond.hpp"

namespace PENKNIFE
{

class WallBC : public PlasmaBaseBndCond
{
public:
    friend class MemoryManager<WallBC>;

    static PlasmaBaseBndCondSharedPtr create(
        const LU::SessionReaderSharedPtr &pSession,
        const Array<OneD, MR::ExpListSharedPtr> &pFields,
        const Array<OneD, MR::DisContFieldSharedPtr> &pB,
        const Array<OneD, MR::DisContFieldSharedPtr> &pE,
        Array<OneD, SpatialDomains::BoundaryConditionShPtr> cond,
        Array<OneD, MultiRegions::ExpListSharedPtr> exp, const int pSpaceDim,
        const int bcRegion)

    {
        PlasmaBaseBndCondSharedPtr p =
            MemoryManager<WallBC>::AllocateSharedPtr(
                pSession, pFields, pB, pE, cond, exp, pSpaceDim, bcRegion);
        return p;
    }

    static std::string className;

protected:
    void v_Apply(const Array<OneD, const Array<OneD, NekDouble>> &Fwd,
                 const Array<OneD, const Array<OneD, NekDouble>> &physarray,
                 const NekDouble &time) override;

private:
    WallBC(const LU::SessionReaderSharedPtr &pSession,
             const Array<OneD, MR::ExpListSharedPtr> &pFields,
             const Array<OneD, MR::DisContFieldSharedPtr> &pB,
             const Array<OneD, MR::DisContFieldSharedPtr> &pE,
             Array<OneD, SpatialDomains::BoundaryConditionShPtr> cond,
             Array<OneD, MultiRegions::ExpListSharedPtr> exp,
             const int pSpaceDim, const int bcRegion);
    ~WallBC() override {};

    // Hardcoded for now
    NekDouble gamma_i = 5 / 2;
    NekDouble gamma_e = 9 / 2;
    NekDouble wall    = 0;
    NekDouble Ge      = 1;
    NekDouble me      = 1 / 2000;

    Array<OneD, Array<OneD, NekDouble>> v_ExB;
};

} // namespace PENKNIFE
#endif