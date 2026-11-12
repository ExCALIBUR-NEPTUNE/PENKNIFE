#pragma once

#include <SolverCore/TraceFlux/TraceFluxOp.hpp>

namespace PENKNIFE
{
using namespace Nektar;

template <typename TData>
class FluxTraceAdvDiffOp : public SolverCore::TraceFluxOp<TData>
{
public:
    using Field = LibUtilities::Field<TData, FieldState::Phys>;
    static inline const std::string name = "FluxTraceAdvDiff";
    static std::shared_ptr<FluxTraceAdvDiffOp> Create(
        const MultiRegions::ExpListSharedPtr &expansionList,
        const std::vector<std::string> &components,
        const std::string &execStr = "", const std::string &model = "")
    {
        std::string closure;
        expansionList->GetSession()->LoadSolverInfo("Closure", closure,
                                                    "Braginskii");
        if (!model.empty())
            closure = model;
        ASSERTL0(closure == "Braginskii" || closure == "SingleDiffusive" ||
                     closure == "DoubleDiffusive",
                 "Only the Braginskii trace closure is implemented.");
        return std::dynamic_pointer_cast<FluxTraceAdvDiffOp>(
            SolverCore::TraceFluxOp<TData>::Create(
                expansionList, components, name + "LaxFriedrichs" + closure,
                execStr));
    }

    // Non-owning element-local trace fields, extracted with PhysTraceExtract.
    // MagneticTrace has 3 components; CurrentTrace has one in number-current
    // units, matching QuasineutralKernel. Keep them fixed through an exchange /
    // Apply / ApplyParallel cycle. Current contents are sampled every cycle.
    virtual void SetMagneticTrace(Field &b)              = 0;
    virtual void SetParallelCurrentTrace(Field &current) = 0;
    virtual void SetDiffusionScale(TData scale)          = 0;
    virtual void SetAdvectionScale(TData scale)          = 0;
    virtual void SetExtraAdvectionTrace(Field &velocities, Field &omegaFlux,
                                        unsigned omega)  = 0;
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

    FluxTraceAdvDiffOp(const MultiRegions::ExpListSharedPtr &expansionList,
                       const std::vector<std::string> &components)
        : SolverCore::TraceFluxOp<TData>(expansionList, components)
    {
    }
};
} // namespace PENKNIFE
