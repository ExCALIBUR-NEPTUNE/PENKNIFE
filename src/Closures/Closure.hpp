#ifndef CLOSURE_HPP
#define CLOSURE_HPP

#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

#include "nektar_interface/utilities.hpp"
#include <MultiRegions/ContField.h>
#include <SolverUtils/UnsteadySystem.h>

namespace SD = Nektar::SpatialDomains;
namespace LU = Nektar::LibUtilities;

namespace PENKNIFE
{
// Forward declarations
class PlasmaSystem;
class Closure;

/// A shared pointer to an equation of state object
typedef std::shared_ptr<Closure> ClosureSharedPtr;

/// Declaration of the equation of state factory
typedef LU::NekFactory<std::string, Closure,
                       const std::weak_ptr<PlasmaSystem> &, const int>
    ClosureFactory;

/// Declaration of the equation of state factory singleton
ClosureFactory &GetClosureFactory();

class Closure
{

public:
    Closure(const std::weak_ptr<PlasmaSystem> &pSystem, const int spaceDim);

    virtual ~Closure() = default;

    void EvaluateClosure(
        const Array<OneD, Array<OneD, NekDouble>> &values,
        const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &grads,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fluxes,
        Array<OneD, Array<OneD, NekDouble>> &frictions,
        const Array<OneD, NekDouble> &ne, const Array<OneD, NekDouble> &ve)
    {
        return v_EvaluateClosure(values, grads, fluxes, frictions, ne, ve);
    };

protected:
    virtual void v_EvaluateClosure(
        const Array<OneD, Array<OneD, NekDouble>> &values,
        const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &grads,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fluxes,
        Array<OneD, Array<OneD, NekDouble>> &frictions,
        const Array<OneD, NekDouble> &ne, const Array<OneD, NekDouble> &ve) = 0;

    const std::weak_ptr<PlasmaSystem> m_system;
    const NektarFieldIndexMap &field_to_index;
    const Array<OneD, Array<OneD, NekDouble>> &b_unit;
    const Array<OneD, NekDouble> &mag_B;
    int omega_idx;
    int ee_idx;
    int n_pts;
    size_t m_spacedim;

    double Nnorm;
    double Tnorm;
    double omega_c;
};

} // namespace PENKNIFE

#endif