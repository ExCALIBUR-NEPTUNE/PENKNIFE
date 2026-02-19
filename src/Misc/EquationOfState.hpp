#ifndef EQUATIONOFSTATE_HPP
#define EQUATIONOFSTATE_HPP

#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/SimdLib/tinysimd.hpp>

using namespace Nektar;
namespace LU = Nektar::LibUtilities;

namespace PENKNIFE
{
//  Forward declaration
class EquationOfState;

/// A shared pointer to an equation of state object
typedef std::shared_ptr<EquationOfState> EquationOfStateSharedPtr;

/// Declaration of the equation of state factory
typedef LU::NekFactory<std::string, EquationOfState> EquationOfStateFactory;

/// Declaration of the equation of state factory singleton
EquationOfStateFactory &GetEquationOfStateFactory();

using vec_t = tinysimd::simd<NekDouble>;

/**
 * @class EquationOfState
 * @brief Encapsulates equations of state allowing us to obtain thermodynamic
 *        properties: most relations are in the form X(rho,e)
 */
class EquationOfState
{
public:
    virtual ~EquationOfState() = default;

    /// Calculate the temperature
    NekDouble GetTemperature(const NekDouble &rho, const NekDouble &e)
    {
        return v_GetTemperature(rho, e);
    }

    vec_t GetTemperature(const vec_t &rho, const vec_t &e)
    {
        return v_GetTemperature(rho, e);
    }

    /// Calculate the pressure
    NekDouble GetPressure(const NekDouble &rho, const NekDouble &e)
    {
        return v_GetPressure(rho, e);
    }

    vec_t GetPressure(const vec_t &rho, const vec_t &e)
    {
        return v_GetPressure(rho, e);
    }

    /// Calculate the sound speed
    NekDouble GetSoundSpeed(const NekDouble &rho, const NekDouble &e)
    {
        return v_GetSoundSpeed(rho, e);
    }

    /// Calculate the entropy
    NekDouble GetEntropy(const NekDouble &rho, const NekDouble &e)
    {
        return v_GetEntropy(rho, e);
    }

    /// Calculate the partial derivative of P(rho,e) with respect to rho
    NekDouble GetDPDrho_e(const NekDouble &rho, const NekDouble &e)
    {
        return v_GetDPDrho_e(rho, e);
    }

    /// Calculate the partial derivative of P(rho,e) with respect to e
    NekDouble GetDPDe_rho(const NekDouble &rho, const NekDouble &e)
    {
        return v_GetDPDe_rho(rho, e);
    }

    /// Obtain the internal energy from rho and P
    NekDouble GetEFromRhoP(const NekDouble &rho, const NekDouble &p)
    {
        return v_GetEFromRhoP(rho, p);
    }

    /// Obtain the density from P and T
    NekDouble GetRhoFromPT(const NekDouble &p, const NekDouble &T)
    {
        return v_GetRhoFromPT(p, T);
    }

protected:
    NekDouble m_gamma;
    NekDouble m_gasConstant;
    NekDouble m_gammaMone;
    NekDouble m_gammaMoneOgasConst;

    /// Constructor
    EquationOfState();

    /// Programmatic Constructor
    EquationOfState(const NekDouble &gamma, const NekDouble &gasConstant);

    virtual NekDouble v_GetTemperature(const NekDouble &rho,
                                       const NekDouble &e) = 0;

    virtual vec_t v_GetTemperature(const vec_t &rho, const vec_t &e) = 0;

    virtual NekDouble v_GetPressure(const NekDouble &rho,
                                    const NekDouble &e) = 0;

    virtual vec_t v_GetPressure(const vec_t &rho, const vec_t &e) = 0;

    virtual NekDouble v_GetSoundSpeed(const NekDouble &rho, const NekDouble &e);

    virtual NekDouble v_GetEntropy(const NekDouble &rho,
                                   const NekDouble &e) = 0;

    virtual NekDouble v_GetDPDrho_e(const NekDouble &rho,
                                    const NekDouble &e) = 0;

    virtual NekDouble v_GetDPDe_rho(const NekDouble &rho,
                                    const NekDouble &e) = 0;

    virtual NekDouble v_GetEFromRhoP(const NekDouble &rho,
                                     const NekDouble &p) = 0;

    virtual NekDouble v_GetRhoFromPT(const NekDouble &rho,
                                     const NekDouble &p) = 0;
};
} // namespace PENKNIFE

#endif