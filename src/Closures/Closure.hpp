#ifndef CLOSURE_HPP
#define CLOSURE_HPP

#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

#include "nektar_interface/utilities.hpp"
#include <MultiRegions/ContField.h>
#include <SolverUtils/UnsteadySystem.h>
#include "../Misc/Constants.hpp"

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

    inline double CoulombLog_ii(double ni1, double ni2, double Ti1, double Ti2,
                                double A1, double A2, double Z1, double Z2)
    {
        return 29.91 - log(sqrt(Nnorm)) -
               log((Z1 * Z2 * (A1 + A2)) / (A1 * Ti2 + A2 * Ti1) *
                   sqrt(ni1 * Z1 * Z1 / Ti1 + ni2 * Z2 * Z2 / Ti2)) +
               1.5 * log(Tnorm);
    }

    inline double CoulombLog_ee(double ne, double Te)
    {
        double logTe = log(Tnorm * Te);
        return 30.4 - 0.5 * log(ne) - 0.5 * log(Nnorm) + (5. / 4) * logTe -
               sqrt(1e-5 + (logTe - 2) * (logTe - 2) / 16.);
    }

    inline double CoulombLog_ei(double ni, double ne, double Ti, double Te,
                                double Ai, double Zi)
    {
        if ((Te * Tnorm < 0.1) || (ni * Nnorm < 1e10) || (ne * Nnorm < 1e10))
            return 10;
        else if (Te * Ai < Ti * constants::m_e_m_p)
            return 23 - 0.5 * log(ni) + 1.5 * log(Ti) - log(Zi * Zi * Ai) -
                   0.5 * log(Nnorm) + 1.5 * log(Tnorm);
        else if (Te * Tnorm < exp(2) * Zi * Zi)
            return 30.0 - 0.5 * log(ne) - log(Zi) + 1.5 * log(Te) -
                   0.5 * log(Nnorm) + 1.5 * log(Tnorm);
        else
            return 31.0 - 0.5 * log(ne) + log(Te) - 0.5 * log(Nnorm) +
                   log(Tnorm);
    }

    inline double CoulombLog(double n1, double n2, double T1, double T2,
                             double A1, double A2, double Z1, double Z2)
    {
        if (A1 == constants::m_e_m_p && A2 == constants::m_e_m_p)
        {
            return CoulombLog_ee(n1, T1);
        }
        else if (A2 == constants::m_e_m_p)
        {
            return CoulombLog_ei(n1, n2, T1, T2, A1, Z1);
        }
        else if (A1 == constants::m_e_m_p)
        {
            return CoulombLog_ei(n2, n1, T2, T1, A2, Z2);
        }
        else
        {
            return CoulombLog_ii(n1, n2, T1, T2, A1, A2, Z1, Z2);
        }
    }

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
    double mesh_length;
};

} // namespace PENKNIFE

#endif