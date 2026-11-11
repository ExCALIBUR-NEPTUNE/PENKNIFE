#include "ElectrostaticTurbulence.hpp"
namespace PENKNIFE
{
std::string ElectrostaticTurbulence::class_name =
    SC::GetEquationSystemFactory().RegisterCreatorFunction(
        "ElectrostaticTurbulence", ElectrostaticTurbulence::create,
        "Electrostatic turbulence with anisotropic diffusion");
ElectrostaticTurbulence::ElectrostaticTurbulence(
    const LU::SessionReaderSharedPtr &session,
    const SD::MeshGraphSharedPtr &graph)
    : ReducedBraginskii(session, graph)
{
    omega_idx = m_variables.size();
    m_variables.push_back("w");
    m_fieldTemplates.push_back("w");
    m_nVariables = m_variables.size();
}
void ElectrostaticTurbulence::LoadParams()
{
    PlasmaSystem::LoadParams();
    std::string boussinesq;
    m_session->LoadSolverInfo("Boussinesq Approximation", boussinesq, "On");
    ASSERTL0(boussinesq == "On",
             "Only Boussinesq Approximation=On is currently implemented.");
}
void ElectrostaticTurbulence::v_InitialiseFields()
{
    ReducedBraginskii::v_InitialiseFields();
    m_phiExpansion = MemoryManager<MR::ContField>::AllocateSharedPtr(
        m_session, m_graph, "phi", false, true);
    m_phiExpansion->SetDataWarehouse();
    auto volume =
        MR::GetBlockAttributes<double, FieldState::Phys>(m_expansionLists[0]);
    auto trace = MR::GetLocTraceBlockAttributes<double, FieldState::Phys>(
        m_expansionLists[0]);
    auto phiAttr =
        MR::GetBlockAttributes<double, FieldState::Phys>(m_phiExpansion);
    phi       = PhysField("phi", phiAttr, 1, 1);
    phi_coeff = CoeffField(
        "phi coefficients",
        MR::GetBlockAttributes<double, FieldState::Coeff>(m_phiExpansion), 1,
        1);
    m_phiOutputRhs = CoeffField(
        "Potential output RHS",
        MR::GetBlockAttributes<double, FieldState::Coeff>(m_phiExpansion), 1,
        1);
    m_phiRhs      = PhysField("Potential forcing", phiAttr, 1, 1);
    m_phiGradient = PhysField("Potential gradient", phiAttr, m_coordim, 1);
    m_perpCoeff   = PhysField("Perpendicular tensor", phiAttr,
                              m_coordim * (m_coordim + 1) / 2, 1);
    m_gradient =
        PhysField("State gradient", volume, m_variables.size() * m_coordim, 1);
    m_extraVelocity      = PhysField("Per-field drift velocities", volume,
                                     m_variables.size() * m_coordim, 1);
    m_extraVelocityTrace = PhysField("Drift velocity traces", trace,
                                     m_variables.size() * m_coordim, 1);
    m_omegaFlux          = PhysField("Vorticity flux", volume, m_coordim, 1);
    m_omegaFluxTrace = PhysField("Vorticity flux traces", trace, m_coordim, 1);
    m_tensor         = PhysField("Polarization velocity tensor", volume,
                                 m_coordim * m_coordim, 1);
    m_tensorGradient = PhysField("Polarization tensor gradient", volume,
                                 m_coordim * m_coordim * m_coordim, 1);
    m_math.zero(phi);
    m_math.zero(phi_coeff);
}
void ElectrostaticTurbulence::v_InitialiseOperators()
{
    ReducedBraginskii::v_InitialiseOperators();
    m_phiSolve  = SC::PoissonSolveOp<double>::Create(m_phiExpansion, {"phi"});
    auto solver = SC::LinearSolverOp<double>::Create(m_phiExpansion, {"phi"});
    m_phiSolve->SetLinearSolver(solver);
    m_phiSolve->SetVarDiffCoeff(m_perpCoeff);
    m_phiSolve->SetPrecon(
        SC::PreconOp<double>::Create(m_phiExpansion, {"phi"}, "Diagonal"));
    m_phiBackward = MR::BwdTransOp<double>::Create(m_phiExpansion, {"phi"});
    m_derivative =
        MR::PhysDerivOp<double>::Create(m_expansionLists[0], m_fieldTemplates);
    m_phiDerivative = MR::PhysDerivOp<double>::Create(m_phiExpansion, {"phi"});
    m_phiOutputInnerProduct =
        MR::IProductWRTBaseOp<double>::Create(m_phiExpansion, {"phi"});
    m_phiOutputInverseMass =
        MR::MultiplyByElmtInvMassOp<double>::Create(m_phiExpansion, {"phi"});
    m_tensorDerivative = MR::PhysDerivOp<double>::Create(
        m_expansionLists[0], m_tensor.GetComponentNames());
    m_advectionFlux->SetExtraFlux(m_extraVelocity, m_omegaFlux, omega_idx);
    m_traceFlux->SetExtraAdvectionTrace(m_extraVelocityTrace, m_omegaFluxTrace,
                                        omega_idx);
}
void ElectrostaticTurbulence::UpdateElectricField(PhysField &in, double time)
{
    m_pointwise->PotentialForcing(in, b_unit, mag_B, m_phiRhs, m_perpCoeff,
                                  omega_idx);
    m_phiSolve->UpdateBndCoeffs(time);
    m_phiSolve->UpdatePrecon();
    m_math.zero(phi_coeff);
    // PoissonSolve uses the physical sign div(D grad(phi)) = rhs.
    m_phiSolve->Apply(m_phiRhs, phi_coeff);
    m_phiBackward->Apply(phi_coeff, phi);
    m_pointwise->CorrectPotential(in, phi);
    m_phiDerivative->Apply(phi, m_phiGradient);
    m_pointwise->ElectricField(m_phiGradient, E);
}
void ElectrostaticTurbulence::PrepareAdvection(PhysField &in, double time)
{
    ReducedBraginskii::PrepareAdvection(in, time);
    UpdateElectricField(in, time);
    m_derivative->Apply(in, m_gradient);
    m_pointwise->DriftVelocities(in, m_gradient, B, mag_B, E, b_unit,
                                 m_extraVelocity, m_velocity, omega_idx);
    // Differentiation is linear: sum mass-weighted ion tensors first.
    // This reuses the tensor workspace and needs one derivative for all ions.
    m_pointwise->PolarizationTensor(in, m_velocity, B, mag_B, m_tensor,
                                    omega_c);
    m_tensorDerivative->Apply(m_tensor, m_tensorGradient);
    m_pointwise->PolarizationFlux(in, m_tensorGradient, m_omegaFlux);
    m_traceExtract->Apply(m_extraVelocity, m_extraVelocityTrace);
    m_traceExtract->Apply(m_omegaFlux, m_omegaFluxTrace);
}
void ElectrostaticTurbulence::CompleteVelocities(PhysField &in)
{
    // Pressure work uses ExB plus parallel flow, before diamagnetic heat flux.
    m_pointwise->CompleteVelocities(in, E, B, mag_B, m_velocity);
}
SC::CoefficientFieldViews ElectrostaticTurbulence::v_GetOutputFieldViews()
{
    auto views = PlasmaSystem::v_GetOutputFieldViews();
    m_phiOutputInnerProduct->Apply(phi, m_phiOutputRhs);
    m_phiOutputInverseMass->Apply(m_phiOutputRhs, phi_coeff);
    views.push_back({phi_coeff, m_phiExpansion, m_coordim, {{"phi", 0}}});
    return views;
}
} // namespace PENKNIFE
