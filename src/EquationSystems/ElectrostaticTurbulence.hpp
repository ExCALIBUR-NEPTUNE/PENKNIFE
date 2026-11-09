#pragma once
#include "ReducedBraginskii.hpp"
#include <MultiRegions/ContField.h>
#include <SolverCore/GlobalLinSysOps/LinearSystems/PoissonSolve/PoissonSolveOp.hpp>
namespace PENKNIFE
{
class ElectrostaticTurbulence : public ReducedBraginskii
{
public:
    friend class MemoryManager<ElectrostaticTurbulence>;
    static std::string class_name;
    static SC::EquationSystemSharedPtr create(
        const LU::SessionReaderSharedPtr &session,
        const SD::MeshGraphSharedPtr &graph)
    {
        auto p = MemoryManager<ElectrostaticTurbulence>::AllocateSharedPtr(
            session, graph);
        p->InitObject();
        return p;
    }
    ~ElectrostaticTurbulence() override = default;

protected:
    ElectrostaticTurbulence(const LU::SessionReaderSharedPtr &session,
                            const SD::MeshGraphSharedPtr &graph);
    void LoadParams() override;
    void v_InitialiseFields() override;
    void v_InitialiseOperators() override;
    void PrepareAdvection(PhysField &in, double time) override;
    void UpdateElectricField(PhysField &in, double time) override;
    SC::CoefficientFieldViews v_GetOutputFieldViews() override;
    void CompleteVelocities(PhysField &in) override;

    unsigned omega_idx;
    PhysField m_gradient;
    std::shared_ptr<MR::PhysDerivOp<double>> m_derivative;
    MR::ExpListSharedPtr m_phiExpansion;
    PhysField phi, m_phiRhs, m_phiGradient, m_perpCoeff;
    CoeffField phi_coeff, m_phiOutputRhs;
    std::shared_ptr<MR::IProductWRTBaseOp<double>> m_phiOutputInnerProduct;
    std::shared_ptr<MR::MultiplyByElmtInvMassOp<double>> m_phiOutputInverseMass;
    PhysField m_extraVelocity, m_extraVelocityTrace, m_omegaFlux,
        m_omegaFluxTrace;
    PhysField m_tensor, m_tensorGradient;
    std::shared_ptr<SC::PoissonSolveOp<double>> m_phiSolve;
    std::shared_ptr<MR::BwdTransOp<double>> m_phiBackward;
    std::shared_ptr<MR::PhysDerivOp<double>> m_phiDerivative,
        m_tensorDerivative;
};
} // namespace PENKNIFE
