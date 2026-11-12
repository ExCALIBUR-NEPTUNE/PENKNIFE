#pragma once

#include <SolverCore/AdvDiffusion/AdvDiffVolumeFluxOp.hpp>
#include <cmath>

namespace PENKNIFE
{
using namespace Nektar;

template <typename TData>
class FluxVolAdvDiffOp : public SolverCore::AdvDiffVolumeFluxOp<TData>
{
public:
    static std::shared_ptr<FluxVolAdvDiffOp<TData>> Create(
        const MultiRegions::ExpListSharedPtr &expansionList,
        const std::vector<std::string> &components,
        const std::string &execStr = "")
    {
        std::string closureName = "Braginskii";
        expansionList->GetSession()->LoadSolverInfo("Closure", closureName,
                                                    "Braginskii");
        ASSERTL0(
            closureName == "Braginskii",
            "Only the Braginskii parallel heat-flux closure is implemented.");
        return std::dynamic_pointer_cast<FluxVolAdvDiffOp<TData>>(
            SolverCore::AdvDiffVolumeFluxOp<TData>::Create(
                expansionList, components, name + closureName, execStr));
    }
    static inline const std::string name = "FluxVolAdvDiff";

    // Non-owning: fields must outlive the operator.
    void SetMagneticField(LibUtilities::Field<TData, FieldState::Phys> &b_unit)
    {
        m_b_unit = &b_unit;
    }
    // Reserved for the separate volume-source calculation. Heat conduction
    // neither reads nor writes this field.
    void SetFriction(LibUtilities::Field<TData, FieldState::Phys> &friction)
    {
        m_friction = &friction;
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

    FluxVolAdvDiffOp(const MultiRegions::ExpListSharedPtr &expansionList,
                     const std::vector<std::string> &components)
        : SolverCore::AdvDiffVolumeFluxOp<TData>(expansionList, components)
    {
    }
    ~FluxVolAdvDiffOp() override = default;

    LibUtilities::Field<TData, FieldState::Phys> *m_ne = nullptr;
    LibUtilities::Field<TData, FieldState::Phys> *m_Me = nullptr;

    LibUtilities::Field<TData, FieldState::Phys> *m_b_unit   = nullptr;
    LibUtilities::Field<TData, FieldState::Phys> *m_friction = nullptr;
};
} // namespace PENKNIFE
