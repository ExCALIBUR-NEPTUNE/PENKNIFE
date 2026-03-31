#ifndef PLASMA_VARIABLECONVERTER_HPP
#define PLASMA_VARIABLECONVERTER_HPP

#include "EquationOfState.hpp"
#include "nektar_interface/utilities.hpp"
#include <MultiRegions/ContField.h>
#include <SolverUtils/UnsteadySystem.h>

namespace SD = Nektar::SpatialDomains;

namespace PENKNIFE
{
// Forward declarations
class PlasmaSystem;
class VariableConverter;
typedef std::shared_ptr<VariableConverter> VariableConverterSharedPtr;
/**
 *
 */
class VariableConverter
{
public:
    VariableConverter(const std::weak_ptr<PlasmaSystem> &pSystem,
                      const int spaceDim);

    ~VariableConverter() = default;

    void GetElectronDensity(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, NekDouble> &density);

    void GetElectronVelocity(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        const Array<OneD, NekDouble> &current,
        const Array<OneD, NekDouble> &density,
        Array<OneD, NekDouble> &velocity);

    void GetElectronDynamicEnergy(const Array<OneD, NekDouble> &velocity,
                                  const Array<OneD, NekDouble> &density,
                                  Array<OneD, NekDouble> &energy);

    // Transformations depending on the equation of state
    void GetElectronTemperature(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, NekDouble> &temperature);
    template <class T, typename = typename std::enable_if<
                           std::is_floating_point<T>::value ||
                           tinysimd::is_vector_floating_point<T>::value>::type>
    inline T GetElectronTemperature(T *physfield)
    {
        T energy = GetElectronInternalEnergy(physfield);
        return m_eos->GetTemperature(physfield[0], energy);
    }
    //
    void GetElectronPressure(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, NekDouble> &pressure);
    template <class T, typename = typename std::enable_if<
                           std::is_floating_point<T>::value ||
                           tinysimd::is_vector_floating_point<T>::value>::type>
    inline T GetElectronPressure(T *physfield)
    {
        T energy = GetElectronInternalEnergy(physfield);
        return m_eos->GetPressure(physfield[0], energy);
    }

    void GetIonDynamicEnergy(
        const int s, const double mass,
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, NekDouble> &energy);
    void GetIonInternalEnergy(
        const int s, const double mass,
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, NekDouble> &energy);

    void GetIonParallelVelocity(
        const int s, const double mass,
        const Array<OneD, Array<OneD, NekDouble>> &physfield,
        Array<OneD, NekDouble> &velocity);

    // Transformations depending on the equation of state
    void GetIonTemperature(
        const int s, const double mass,
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, NekDouble> &temperature);
    template <class T, typename = typename std::enable_if<
                           std::is_floating_point<T>::value ||
                           tinysimd::is_vector_floating_point<T>::value>::type>
    inline T GetIonTemperature(const int s, const double mass, T *physfield)
    {
        T energy = GetIonInternalEnergy(s, mass, physfield);
        return m_eos->GetTemperature(physfield[0], energy);
    }
    //
    void GetIonPressure(
        const int s, const double mass,
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, NekDouble> &pressure);
    template <class T, typename = typename std::enable_if<
                           std::is_floating_point<T>::value ||
                           tinysimd::is_vector_floating_point<T>::value>::type>
    inline T GetIonPressure(const int s, const double mass, T *physfield)
    {
        T energy = GetIonInternalEnergy(s, mass, physfield);
        return m_eos->GetPressure(physfield[0], energy);
    }
    void GetIonSoundSpeed(
        const int s, const double mass, const int Z,
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, NekDouble> &soundspeed);
    void GetInternalEnergy(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, NekDouble> &energy);
    void GetSystemSoundSpeed(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, NekDouble> &soundspeed);
    const EquationOfStateSharedPtr Geteos()
    {
        return m_eos;
    }
    int omega_idx;
    int ee_idx;

protected:
    const std::weak_ptr<PlasmaSystem> m_system;
    const NektarFieldIndexMap &field_to_index;

    EquationOfStateSharedPtr m_eos;
    size_t m_spacedim;
};

} // namespace PENKNIFE

#endif