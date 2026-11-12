#pragma once

#include <SolverCore/Diffusion/DiffusionVolumeFluxOp.hpp>
#include <cmath>

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
    
    FluxVolDiffOp(const MultiRegions::ExpListSharedPtr &expansionList,
                  const std::vector<std::string> &components)
        : SolverCore::DiffusionVolumeFluxOp<TData>(expansionList, components)
    {
    }
    ~FluxVolDiffOp() override = default;

    bool m_requireEnergy = true;

    LibUtilities::Field<TData, FieldState::Phys> *m_b_unit   = nullptr;
    LibUtilities::Field<TData, FieldState::Phys> *m_friction = nullptr;
};
} // namespace PENKNIFE
