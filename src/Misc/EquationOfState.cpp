#include "EquationOfState.hpp"

namespace PENKNIFE
{
EquationOfStateFactory &GetEquationOfStateFactory()
{
    static EquationOfStateFactory instance;
    return instance;
}

EquationOfState::EquationOfState()
{
    m_gamma = 1.4;
    m_gasConstant= 287.058;
    m_gammaMone          = m_gamma - 1.0;
    m_gammaMoneOgasConst = m_gammaMone / m_gasConstant;
}

EquationOfState::EquationOfState(const NekDouble &gamma,
                                 const NekDouble &gasConstant)
    : m_gamma{gamma}, m_gasConstant{gasConstant}
{
}

// General implementation for v_GetSoundSpeed: c^2 = xi + kappa * h
//    where xi = dpdrho - e/rho * dp/de    and  kappa = dp/de / rho
NekDouble EquationOfState::v_GetSoundSpeed(const NekDouble &rho,
                                           const NekDouble &e)
{
    NekDouble p      = GetPressure(rho, e);
    NekDouble dpde   = GetDPDe_rho(rho, e);
    NekDouble dpdrho = GetDPDrho_e(rho, e);

    NekDouble enthalpy = e + p / rho;

    NekDouble chi   = dpdrho - e / rho * dpde;
    NekDouble kappa = dpde / rho;

    return std::sqrt(chi + kappa * enthalpy);
}

} // namespace PENKNIFE
