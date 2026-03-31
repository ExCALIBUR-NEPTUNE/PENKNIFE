#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

namespace constants
{
static constexpr double m_e     = 0.51099895069e+6;
static constexpr double m_p     = 938.272089e+6;
static constexpr double m_e_m_p = m_e / m_p;
static constexpr double m_p_si  = 1.6726219259552e-27;
static constexpr double m_e_si  = m_p_si * m_e_m_p;
static constexpr double c       = 299792458;

// epsilon_0 in e^2 eV^-1 um^-1
static constexpr double e            = 1.602176634e-19;
static constexpr double epsilon_0    = 55.26349406;
static constexpr double epsilon_0_si = 8.8541878188e-12;

static constexpr double qeomp = e / m_p_si;

static constexpr double k_B     = 1.380649e-23;
static constexpr double temp_SI = 11604.51812;

static constexpr double mass_amu    = 1.0;
static constexpr double mass_amu_SI = 1.66053904e-27;

} // namespace constants

#endif