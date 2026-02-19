#include "IdealGasEoS.hpp"

namespace PENKNIFE
{

std::string IdealGasEoS::className =
    GetEquationOfStateFactory().RegisterCreatorFunction(
        "IdealGas", IdealGasEoS::create, "Ideal gas equation of state.");

IdealGasEoS::IdealGasEoS() : EquationOfState()
{
}

NekDouble IdealGasEoS::v_GetTemperature([[maybe_unused]] const NekDouble &rho,
                                        const NekDouble &e)
{
    return GetTemperatureKernel(rho, e);
}

vec_t IdealGasEoS::v_GetTemperature([[maybe_unused]] const vec_t &rho,
                                    const vec_t &e)
{
    return GetTemperatureKernel(rho, e);
}

NekDouble IdealGasEoS::v_GetPressure(const NekDouble &rho, const NekDouble &e)
{
    return GetPressureKernel(rho, e);
}

vec_t IdealGasEoS::v_GetPressure(const vec_t &rho, const vec_t &e)
{
    return GetPressureKernel(rho, e);
}

NekDouble IdealGasEoS::v_GetSoundSpeed(const NekDouble &rho, const NekDouble &e)
{
    NekDouble T = GetTemperature(rho, e);
    return sqrt(m_gamma * m_gasConstant * T);
}

NekDouble IdealGasEoS::v_GetEntropy(const NekDouble &rho, const NekDouble &e)
{
    NekDouble T = GetTemperature(rho, e);
    return m_gasConstant / (m_gamma - 1) * log(T) - m_gasConstant * log(rho);
}

NekDouble IdealGasEoS::v_GetDPDrho_e([[maybe_unused]] const NekDouble &rho,
                                     const NekDouble &e)
{
    return e * (m_gamma - 1);
}

NekDouble IdealGasEoS::v_GetDPDe_rho(const NekDouble &rho,
                                     [[maybe_unused]] const NekDouble &e)
{
    return rho * (m_gamma - 1);
}

NekDouble IdealGasEoS::v_GetEFromRhoP(const NekDouble &rho, const NekDouble &p)
{
    return p / (rho * (m_gamma - 1));
}

NekDouble IdealGasEoS::v_GetRhoFromPT(const NekDouble &p, const NekDouble &T)
{
    return p / (m_gasConstant * T);
}

} // namespace PENKNIFE
