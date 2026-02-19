#ifndef PLASMA_OBLIQUEOUTFLOW_HPP
#define PLASMA_OBLIQUEOUTFLOW_HPP

#include "PlasmaBaseBndCond.hpp"

namespace PENKNIFE
{

class ObliqueOutflowBC : public PlasmaBaseBndCond
{
public:
    friend class MemoryManager<ObliqueOutflowBC>;

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
            MemoryManager<ObliqueOutflowBC>::AllocateSharedPtr(
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
    ObliqueOutflowBC(const LU::SessionReaderSharedPtr &pSession,
                     const std::weak_ptr<PlasmaSystem> &pSystem,
                     const Array<OneD, MR::ExpListSharedPtr> &pFields,
                     const Array<OneD, MR::DisContFieldSharedPtr> &pB,
                     const Array<OneD, MR::DisContFieldSharedPtr> &pE,
                     Array<OneD, SpatialDomains::BoundaryConditionShPtr> cond,
                     Array<OneD, MultiRegions::ExpListSharedPtr> exp,
                     const int pSpaceDim, const int bcRegion);
    ~ObliqueOutflowBC() override {};

    void CalcKPar();
    void CalcKPerp();
    void CalcDTensor();

    void CalcKappaPar();
    void CalcKappaPerp();
    void CalcKappaTensor();

    Array<OneD, NekDouble> m_D[3][3];
    Array<OneD, NekDouble> m_kappa[3][3];
    Array<OneD, NekDouble> kpar;
    Array<OneD, NekDouble> kperp;
    Array<OneD, NekDouble> kappapar;
    Array<OneD, NekDouble> kappaperp;
    NekDouble gamma;
    NekDouble m_i;
    NekDouble k_B;
    NekDouble k_par;
    NekDouble k_perp;
    NekDouble kappa_par;
    NekDouble kappa_perp;
};

} // namespace PENKNIFE
#endif