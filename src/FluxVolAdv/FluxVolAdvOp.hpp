#pragma once
#include <SolverCore/Advection/AdvectionVolumeFluxOp.hpp>

namespace PENKNIFE
{
using namespace Nektar;
template <typename TData>
class FluxVolAdvOp : public SolverCore::AdvectionVolumeFluxOp<TData>
{
public:
    static std::shared_ptr<FluxVolAdvOp<TData>> Create(
        const MultiRegions::ExpListSharedPtr &expansionList,
        const std::vector<std::string> &components,
        const std::string &execStr = "")
    {
        return std::dynamic_pointer_cast<FluxVolAdvOp<TData>>(
            SolverCore::AdvectionVolumeFluxOp<TData>::Create(
                expansionList, components, execStr));
    }
    static inline const std::string name = "FluxVolAdv";

    // Non-owning: fields must outlive the operator.
    void SetMagneticField(LibUtilities::Field<TData, FieldState::Phys> &b_unit)
    {
        m_b_unit = &b_unit;
    }

    void SetExtraFlux(LibUtilities::Field<TData, FieldState::Phys> &velocity,
                      LibUtilities::Field<TData, FieldState::Phys> &flux,
                      unsigned omega)
    {
        extraVelocity = &velocity;
        omegaFlux     = &flux;
        omegaId       = omega;
    }
    void SetElectronFields(LibUtilities::Field<TData, FieldState::Phys> &ne,
                           LibUtilities::Field<TData, FieldState::Phys> &Me)
    {
        m_ne = &ne;
        m_Me = &Me;
    }
    void SetSpecies(const std::vector<unsigned> &density,
                    const std::vector<unsigned> &momentum,
                    const std::vector<unsigned> &energy,
                    const std::vector<TData> &masses,
                    const std::vector<TData> &charges)
    {
        this->v_SetSpecies(density, momentum, energy, masses, charges);
    }
    void SetElectronEnergyId(unsigned id)
    {
        this->v_SetElectronEnergyId(id);
    }

protected:
    virtual void v_SetSpecies(const std::vector<unsigned> &density,
                              const std::vector<unsigned> &momentum,
                              const std::vector<unsigned> &energy,
                              const std::vector<TData> &masses,
                              const std::vector<TData> &charges) = 0;
    virtual void v_SetElectronEnergyId(unsigned id)              = 0;

    FluxVolAdvOp(const MultiRegions::ExpListSharedPtr &expansionList,
                 const std::vector<std::string> &components)
        : SolverCore::AdvectionVolumeFluxOp<TData>(expansionList, components)
    {
    }
    ~FluxVolAdvOp() override = default;

    LibUtilities::Field<TData, FieldState::Phys> *m_ne = nullptr;
    LibUtilities::Field<TData, FieldState::Phys> *m_Me = nullptr;

    LibUtilities::Field<TData, FieldState::Phys> *extraVelocity = nullptr;
    LibUtilities::Field<TData, FieldState::Phys> *omegaFlux     = nullptr;

    LibUtilities::Field<TData, FieldState::Phys> *m_b_unit   = nullptr;
    LibUtilities::Field<TData, FieldState::Phys> *m_friction = nullptr;
    unsigned omegaId;
};
} // namespace PENKNIFE
