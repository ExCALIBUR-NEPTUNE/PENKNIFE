#pragma once

#include "../FluxTraceAdvDiv/FluxTraceAdvDiffOp.hpp"
#include "../FluxVolAdv/FluxVolAdvOp.hpp"
#include "PlasmaSystem.hpp"
#include <MultiRegions/ElmtOps/PhysDeriv/PhysDerivOp.hpp>
#include <SolverCore/AdvDiffusion/AdvWeakDGDiffusionIP/AdvWeakDGDiffusionIPOp.hpp>
namespace PENKNIFE
{
class ReducedBraginskii : public PlasmaSystem
{
public:
    friend class MemoryManager<ReducedBraginskii>;
    static std::string class_name;
    static SC::EquationSystemSharedPtr create(
        const LU::SessionReaderSharedPtr &session,
        const SD::MeshGraphSharedPtr &graph)
    {
        auto p =
            MemoryManager<ReducedBraginskii>::AllocateSharedPtr(session, graph);
        p->InitObject();
        return p;
    }
    ~ReducedBraginskii() override = default;

protected:
    ReducedBraginskii(const LU::SessionReaderSharedPtr &session,
                      const SD::MeshGraphSharedPtr &graph)
        : PlasmaSystem(session, graph)
    {
    }
    void v_InitObject(bool declareFields = true) override;
    void v_InitialiseFields() override;
    void v_InitialiseOperators() override;
    SolverCore::TraceFluxOp<double>::ViscousTraceForm ReadViscousTraceForm();
    virtual void DoOdeRhs(PhysField &in, PhysField &out, const double &time,
                          const double &factor);
    virtual void PrepareAdvection(PhysField &in, double time);
    virtual void CompleteVelocities(PhysField &)
    {
    }
    void AddPressureWork(PhysField &in, PhysField &out, double factor);

    std::shared_ptr<SC::AdvWeakDGDiffusionIPOp<double>> m_advDiffusionOp;
    std::shared_ptr<FluxTraceAdvDiffOp<double>> m_traceFlux;
    std::shared_ptr<FluxVolAdvOp<double>> m_advectionFlux;
    std::shared_ptr<MR::PhysTraceExtractOp<double>> m_traceExtract;
    std::shared_ptr<MR::PhysDerivOp<double>> m_velocityDerivative;
    PhysField m_currentTrace, m_velocity, m_velocityGradient;
};
} // namespace PENKNIFE
