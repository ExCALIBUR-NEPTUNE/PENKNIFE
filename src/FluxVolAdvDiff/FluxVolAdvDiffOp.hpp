#pragma once

#include <SolverCore/AdvDiffusion/AdvDiffVolumeFluxOp.hpp>
#include <cmath>
#include <limits>

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
    static constexpr unsigned InvalidId  = std::numeric_limits<unsigned>::max();

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
                           LibUtilities::Field<TData, FieldState::Phys> &Me,
                           TData electronMass)
    {
        ASSERTL0(std::isfinite(electronMass) && electronMass > 0,
                 "Electron mass must be finite and positive.");
        m_ne           = &ne;
        m_Me           = &Me;
        m_electronMass = electronMass;
    }
    void SetElectronEnergyId(unsigned id)
    {
        ASSERTL0(id < this->m_components.size(),
                 "Invalid electron energy component.");
        m_electronEnergyId = id;
    }
    void SetIndices(const std::vector<unsigned> &n_ids,
                    const std::vector<unsigned> &v_ids,
                    const std::vector<unsigned> &e_ids,
                    const std::vector<TData> &masses,
                    const std::vector<TData> &charges)
    {
        ASSERTL0(n_ids.size() == v_ids.size() && n_ids.size() == e_ids.size() &&
                     n_ids.size() == masses.size() &&
                     n_ids.size() == charges.size(),
                 "Diffusion species metadata vectors must be aligned.");
        for (size_t s = 0; s < n_ids.size(); ++s)
        {
            ASSERTL0(std::isfinite(masses[s]) && masses[s] > 0 &&
                         std::isfinite(charges[s]) && charges[s] >= 0,
                     "Diffusion requires positive finite masses and "
                     "nonnegative finite charges.");
            for (const unsigned id : {n_ids[s], v_ids[s], e_ids[s]})
                ASSERTL0(id == InvalidId || id < this->m_components.size(),
                         "Diffusion species index is outside the state.");
            ASSERTL0(charges[s] == 0 ||
                         (n_ids[s] != InvalidId && e_ids[s] != InvalidId),
                     "Braginskii conduction requires density and energy for "
                     "every ion.");
        }
        v_SetIndices(n_ids, v_ids, e_ids, masses, charges);
        m_indicesConfigured = true;
    }

protected:
    FluxVolAdvDiffOp(const MultiRegions::ExpListSharedPtr &expansionList,
                     const std::vector<std::string> &components)
        : SolverCore::AdvDiffVolumeFluxOp<TData>(expansionList, components)
    {
    }
    ~FluxVolAdvDiffOp() override                             = default;
    virtual void v_SetIndices(const std::vector<unsigned> &,
                              const std::vector<unsigned> &,
                              const std::vector<unsigned> &,
                              const std::vector<TData> &,
                              const std::vector<TData> &)    = 0;
    unsigned m_electronEnergyId                              = InvalidId;
    bool m_indicesConfigured                                 = false;
    LibUtilities::Field<TData, FieldState::Phys> *m_ne       = nullptr;
    LibUtilities::Field<TData, FieldState::Phys> *m_Me       = nullptr;
    TData m_electronMass                                     = 0;
    LibUtilities::Field<TData, FieldState::Phys> *m_b_unit   = nullptr;
    LibUtilities::Field<TData, FieldState::Phys> *m_friction = nullptr;
};
} // namespace PENKNIFE
