#ifndef PLASMA_IDEALGASEOS_HPP
#define PLASMA_IDEALGASEOS_HPP

#include "EquationOfState.hpp"

namespace PENKNIFE
{

/**
 * @brief Ideal gas equation of state:
 *       p = rho * R * T
 */
class IdealGasEoS : public EquationOfState
{
public:
    friend class MemoryManager<IdealGasEoS>;

    /// Creates an instance of this class
    static EquationOfStateSharedPtr create()
    {
        EquationOfStateSharedPtr p =
            MemoryManager<IdealGasEoS>::AllocateSharedPtr();
        return p;
    }

    /// Name of the class
    static std::string className;

protected:
    NekDouble v_GetTemperature(const NekDouble &rho, const NekDouble &e) final;

    vec_t v_GetTemperature(const vec_t &rho, const vec_t &e) final;

    NekDouble v_GetPressure(const NekDouble &rho, const NekDouble &e) final;

    vec_t v_GetPressure(const vec_t &rho, const vec_t &e) final;

    NekDouble v_GetSoundSpeed(const NekDouble &rho, const NekDouble &e) final;

    NekDouble v_GetEntropy(const NekDouble &rho, const NekDouble &e) final;

    NekDouble v_GetDPDrho_e(const NekDouble &rho, const NekDouble &e) final;

    NekDouble v_GetDPDe_rho(const NekDouble &rho, const NekDouble &e) final;

    NekDouble v_GetEFromRhoP(const NekDouble &rho, const NekDouble &p) final;

    NekDouble v_GetRhoFromPT(const NekDouble &rho, const NekDouble &p) final;

private:
    IdealGasEoS();

    ~IdealGasEoS(void) override = default;

    // type agnostic kernels
    template <class T, typename = typename std::enable_if<
                           std::is_floating_point<T>::value ||
                           tinysimd::is_vector_floating_point<T>::value>::type>
    inline T GetTemperatureKernel(const T &rho, const T &e)
    {
        return (2.0 / 3.0) * e / rho;
    }

    template <class T, typename = typename std::enable_if<
                           std::is_floating_point<T>::value ||
                           tinysimd::is_vector_floating_point<T>::value>::type>
    inline T GetPressureKernel(const T &rho, const T &e)
    {
        return rho * e * m_gammaMone;
    }
};

} // namespace PENKNIFE

#endif
