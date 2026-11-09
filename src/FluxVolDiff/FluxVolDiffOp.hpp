#pragma once

#include <SolverCore/Diffusion/DiffusionVolumeFluxOp.hpp>
#include <cmath>
#include <limits>

namespace PENKNIFE
{
using namespace Nektar;

template <typename TData>
class FluxVolDiffOp : public SolverCore::DiffusionVolumeFluxOp<TData>
{
public:
    static std::shared_ptr<FluxVolDiffOp<TData>> Create(
        const MultiRegions::ExpListSharedPtr &expansionList,
        const std::vector<std::string> &components,
        const std::string &execStr = "", const std::string &model = "")
    {
        std::string closureName = "Braginskii";
        expansionList->GetSession()->LoadSolverInfo("Closure", closureName,
                                                    "Braginskii");
        if (!model.empty())
            closureName = model;
        ASSERTL0(
            closureName == "Braginskii" || closureName == "SingleDiffusive" ||
                closureName == "DoubleDiffusive",
            "Only the Braginskii parallel heat-flux closure is implemented.");
        auto op = std::dynamic_pointer_cast<FluxVolDiffOp<TData>>(
            SolverCore::DiffusionVolumeFluxOp<TData>::Create(
                expansionList, components, name + closureName, execStr));
        op->m_requireEnergy = closureName != "SingleDiffusive";
        return op;
    }
    static inline const std::string name = "FluxVolDiff";
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
            ASSERTL0(!m_requireEnergy || charges[s] == 0 ||
                         (n_ids[s] != InvalidId && e_ids[s] != InvalidId),
                     "Braginskii conduction requires density and energy for "
                     "every ion.");
        }
        v_SetIndices(n_ids, v_ids, e_ids, masses, charges);
        m_indicesConfigured = true;
    }

protected:
    FluxVolDiffOp(const MultiRegions::ExpListSharedPtr &expansionList,
                  const std::vector<std::string> &components)
        : SolverCore::DiffusionVolumeFluxOp<TData>(expansionList, components)
    {
    }
    ~FluxVolDiffOp() override                                = default;
    virtual void v_SetIndices(const std::vector<unsigned> &,
                              const std::vector<unsigned> &,
                              const std::vector<unsigned> &,
                              const std::vector<TData> &,
                              const std::vector<TData> &)    = 0;
    bool m_requireEnergy                                     = true;
    unsigned m_electronEnergyId                              = InvalidId;
    bool m_indicesConfigured                                 = false;
    LibUtilities::Field<TData, FieldState::Phys> *m_b_unit   = nullptr;
    LibUtilities::Field<TData, FieldState::Phys> *m_friction = nullptr;
};
} // namespace PENKNIFE
