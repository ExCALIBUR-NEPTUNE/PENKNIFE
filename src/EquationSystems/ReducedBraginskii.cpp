#include "ReducedBraginskii.hpp"
namespace PENKNIFE
{
std::string ReducedBraginskii::class_name =
    SC::GetEquationSystemFactory().RegisterCreatorFunction(
        "ReducedBraginskii", ReducedBraginskii::create,
        "Reduced Braginskii equations");
void ReducedBraginskii::v_InitObject(bool declareFields)
{
    PlasmaSystem::v_InitObject(declareFields);
    m_timeOp->DefineExplicitRhs(&ReducedBraginskii::DoOdeRhs, this);
}
void ReducedBraginskii::v_InitialiseFields()
{
    PlasmaSystem::v_InitialiseFields();
    auto attr =
        MR::GetBlockAttributes<double, FieldState::Phys>(m_expansionLists[0]);
    auto trace = MR::GetLocTraceBlockAttributes<double, FieldState::Phys>(
        m_expansionLists[0]);
    m_currentTrace = PhysField("Parallel current trace", trace, 1, 1);

    m_velocity =
        PhysField("Species velocities", attr, (n_species + 1) * m_coordim, 1);
    m_velocityGradient = PhysField("Velocity gradients", attr,
                                   (n_species + 1) * m_coordim * m_coordim, 1);
}
void ReducedBraginskii::v_InitialiseOperators()
{
    PlasmaSystem::v_InitialiseOperators();
    ASSERTL0(m_projectionType == MR::eDiscontinuous,
             "Reduced Braginskii requires DG projection.");
    m_advDiffusionOp = SC::AdvWeakDGDiffusionIPOp<double>::Create(
        m_expansionLists[0], m_fieldTemplates);
    m_traceFlux = FluxTraceAdvDiffOp<double>::Create(m_expansionLists[0],
                                                     m_fieldTemplates);
    m_traceFlux->SetIndices(n_ids, v_ids, e_ids, m_masses, m_charges);
    m_traceFlux->SetElectronEnergyIndex(ee_idx);
    m_traceFlux->SetMagneticTrace(b_unit_trace);
    m_traceFlux->SetParallelCurrentTrace(m_currentTrace);
    std::string form;
    m_session->LoadSolverInfo("ViscousTraceFlux", form, "FluxOfAverages");
    ASSERTL0(form == "FluxOfAverages" || form == "AverageOfNormalFluxes",
             "Unknown ViscousTraceFlux.");
    if (form == "AverageOfNormalFluxes")
        m_traceFlux->SetViscousTraceForm(
            SC::TraceFluxOp<double>::ViscousTraceForm::AverageOfNormalFluxes);
    m_advDiffusionOp->SetAdvDiffTraceFlux(m_traceFlux);
    m_advectionFlux =
        FluxVolAdvOp<double>::Create(m_expansionLists[0], m_fieldTemplates);
    m_advectionFlux->SetMagneticField(b_unit);
    m_advectionFlux->SetElectronFields(ne, Me, me, ee_idx);
    m_advectionFlux->SetIndices(n_ids, v_ids, e_ids, m_masses, m_charges);
    m_advDiffusionOp->SetAdvVolFlux(m_advectionFlux);
    auto diffusion =
        FluxVolDiffOp<double>::Create(m_expansionLists[0], m_fieldTemplates);
    diffusion->SetMagneticField(b_unit);
    diffusion->SetElectronEnergyId(ee_idx);
    diffusion->SetIndices(n_ids, v_ids, e_ids, m_masses, m_charges);
    m_advDiffusionOp->SetDiffVolFlux(diffusion);
    m_traceExtract = MR::PhysTraceExtractOp<double>::Create(m_expansionLists[0],
                                                            m_fieldTemplates);

    m_velocityDerivative = MR::PhysDerivOp<double>::Create(
        m_expansionLists[0], m_velocity.GetComponentNames());
}
void ReducedBraginskii::PrepareAdvection(PhysField &in, double time)
{
    if (transient_field)
        UpdateMagneticField(time);
    j_par.SetInterleaveWidth(in);
    m_pointwise->Quasineutrality(in, j_par, ne, Me);
    m_traceExtract->Apply(j_par, m_currentTrace);
}
void ReducedBraginskii::DoOdeRhs(PhysField &in, PhysField &out,
                                 const double &time, const double &factor)
{
    PrepareAdvection(in, time);
    m_traceFlux->UpdateBndPhys(time);
    m_advDiffusionOp->SetScale(factor);
    m_advDiffusionOp->Apply(in, out);
    AddPressureWork(in, out, factor);
    AddSources(in, out, time, factor);
}
void ReducedBraginskii::AddPressureWork(PhysField &in, PhysField &out,
                                        double factor)
{
    m_pointwise->ParallelVelocities(in, b_unit, ne, Me, m_velocity);
    CompleteVelocities(in);
    m_velocityDerivative->Apply(m_velocity, m_velocityGradient);
    m_pointwise->PressureWork(in, out, m_velocityGradient, E, b_unit, factor);
}
} // namespace PENKNIFE
