#ifndef CLOSURE_HPP
#define CLOSURE_HPP

#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <MultiRegions/ContField.h>
#include <SolverUtils/UnsteadySystem.h>

#include "nektar_interface/utilities.hpp"

namespace PENKNIFE
{
namespace SD = Nektar::SpatialDomains;
namespace LU = Nektar::LibUtilities;
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

    void CollisionFrequencies(const Array<OneD, Array<OneD, NekDouble>> &values,
                              const Array<OneD, NekDouble> &ne)
    {
        return v_CollisionFrequencies(values, ne);
    }

    void EvaluateHeatFlux(
        const Array<OneD, Array<OneD, NekDouble>> &values,
        const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &grads,
        const Array<OneD, NekDouble> &ne,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fluxes)
    {
        return v_EvaluateHeatFlux(values, grads, ne, fluxes);
    }

    void EvaluateThermalForce(
        const Array<OneD, Array<OneD, NekDouble>> &values,
        const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &grads,
        const Array<OneD, NekDouble> &ne,
        Array<OneD, Array<OneD, NekDouble>> &force)
    {
        return v_EvaluateThermalForce(values, grads, ne, force);
    }

    void EvaluateFrictionHeating(
        const Array<OneD, Array<OneD, NekDouble>> &values,
        const Array<OneD, NekDouble> &ne, const Array<OneD, NekDouble> &ve,
        Array<OneD, Array<OneD, NekDouble>> &frictions,
        Array<OneD, Array<OneD, NekDouble>> &heats)
    {
        return v_EvaluateFrictionHeating(values, ne, ve, frictions, heats);
    }

    void EvaluateConductivity(const Array<OneD, NekDouble> &ne,
                              Array<OneD, NekDouble> &sigma)
    {
        return v_EvaluateConductivity(ne, sigma);
    }

    int ee_idx;

protected:
    virtual void v_CollisionFrequencies(
        const Array<OneD, Array<OneD, NekDouble>> &values,
        const Array<OneD, NekDouble> &ne) = 0;

    virtual void v_EvaluateHeatFlux(
        const Array<OneD, Array<OneD, NekDouble>> &values,
        const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &grads,
        const Array<OneD, NekDouble> &ne,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fluxes) = 0;

    virtual void v_EvaluateThermalForce(
        const Array<OneD, Array<OneD, NekDouble>> &values,
        const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &grads,
        const Array<OneD, NekDouble> &ne,
        Array<OneD, Array<OneD, NekDouble>> &force) = 0;

    virtual void v_EvaluateFrictionHeating(
        const Array<OneD, Array<OneD, NekDouble>> &values,
        const Array<OneD, NekDouble> &ne, const Array<OneD, NekDouble> &ve,
        Array<OneD, Array<OneD, NekDouble>> &frictions,
        Array<OneD, Array<OneD, NekDouble>> &heats) = 0;

    virtual void v_EvaluateConductivity(const Array<OneD, NekDouble> &ne,
                                        Array<OneD, NekDouble> &sigma) = 0;

    const std::weak_ptr<PlasmaSystem> m_system;
    const NektarFieldIndexMap &field_to_index;
    const Array<OneD, Array<OneD, NekDouble>> &b_unit;
    const Array<OneD, NekDouble> &mag_B;
    int omega_idx;

    int n_pts;
    size_t m_spacedim;

    double Nnorm;
    double Tnorm;
    double Bnorm;
    double omega_c;
    double mesh_length;
    double scaling;
};

} // namespace PENKNIFE

#endif