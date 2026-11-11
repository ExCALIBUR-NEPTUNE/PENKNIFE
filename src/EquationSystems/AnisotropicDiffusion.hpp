#pragma once
#include "../FluxTraceAdvDiv/FluxTraceAdvDiffOp.hpp"
#include "../FluxVolDiff/FluxVolDiffOp.hpp"
#include "PlasmaSystem.hpp"
#include <MultiRegions/BndCondOps/NeuBndCond/NeuBndCondOp.hpp>
#include <MultiRegions/ElmtOps/IProductWRTDerivBase/IProductWRTDerivBaseOp.hpp>
#include <MultiRegions/ElmtOps/PhysDeriv/PhysDerivOp.hpp>
#include <SolverCore/Diffusion/DiffusionIP/DiffusionIPOp.hpp>
namespace PENKNIFE
{
class AnisotropicDiffusion : public PlasmaSystem
{
public:
    friend class MemoryManager<AnisotropicDiffusion>;
    static std::string class_name;
    static SC::EquationSystemSharedPtr create(
        const LU::SessionReaderSharedPtr &session,
        const SD::MeshGraphSharedPtr &graph)
    {
        auto p = MemoryManager<AnisotropicDiffusion>::AllocateSharedPtr(session,
                                                                        graph);
        p->InitObject();
        return p;
    }
    ~AnisotropicDiffusion() override = default;

protected:
    AnisotropicDiffusion(const LU::SessionReaderSharedPtr &session,
                         const SD::MeshGraphSharedPtr &graph);
    void v_InitObject(bool declareFields = true) override;
    void v_InitialiseFields() override;
    void v_InitialiseOperators() override;
    std::string DiffusionModel() const;
    void DoOdeRhs(PhysField &in, PhysField &out, const double &time,
                  const double &factor);

    PhysField m_magneticTrace, m_currentTrace, m_gradient, m_flux;
    CoeffField m_rhs, m_massSolution, m_stateCoeffs;
    std::shared_ptr<MR::RobBndCondOp<double>> m_robin;
    std::shared_ptr<MR::IProductWRTBaseOp<double>> m_localInnerProduct;
    std::shared_ptr<MR::MultiplyByElmtInvMassOp<double>> m_localInverseMass;
    std::shared_ptr<SC::DiffusionIPOp<double>> m_diffusionOp;
    std::shared_ptr<FluxVolDiffOp<double>> m_volumeFlux;
    std::shared_ptr<FluxTraceAdvDiffOp<double>> m_traceFlux;
    std::shared_ptr<MR::PhysTraceExtractOp<double>> m_traceExtract;
    std::shared_ptr<MR::PhysDerivOp<double>> m_derivative;
    std::shared_ptr<MR::IProductWRTDerivBaseOp<FieldState::Coeff, double>>
        m_weakDerivative;
    std::shared_ptr<MR::NeuBndCondOp<double>> m_neumann;
    std::shared_ptr<SC::LinearSolverOp<double>> m_massSolver;
    std::string m_diffusionModel;
};
} // namespace PENKNIFE
