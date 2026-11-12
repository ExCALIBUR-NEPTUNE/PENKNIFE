#include "AnisotropicDiffusion.hpp"
#include <MultiRegions/ElmtOps/Mass/MassOp.hpp>
#include <SolverCore/GlobalLinSysOps/LinearSystems/LHSOp.hpp>

namespace PENKNIFE
{
std::string AnisotropicDiffusion::class_name =
    SC::GetEquationSystemFactory().RegisterCreatorFunction(
        "AnisotropicDiffusion", AnisotropicDiffusion::create,
        "Anisotropic density and heat diffusion");

AnisotropicDiffusion::AnisotropicDiffusion(
    const LU::SessionReaderSharedPtr &session,
    const SD::MeshGraphSharedPtr &graph)
    : PlasmaSystem(session, graph)
{
    m_session->LoadSolverInfo("Closure", m_diffusionModel);
    ASSERTL0(m_diffusionModel == "SingleDiffusive" ||
                 m_diffusionModel == "DoubleDiffusive",
             "AnisotropicDiffusion requires Closure=SingleDiffusive or "
             "Closure=DoubleDiffusive.");
    if (m_diffusionModel == "SingleDiffusive")
    {
        m_variables.pop_back();
        m_fieldTemplates.pop_back();
        ee_idx       = InvalidId;
        m_nVariables = m_variables.size();
    }
}
std::string AnisotropicDiffusion::DiffusionModel() const
{
    return m_diffusionModel;
}
void AnisotropicDiffusion::v_InitObject(bool declareFields)
{
    PlasmaSystem::v_InitObject(declareFields);
    m_timeOp->DefineExplicitRhs(&AnisotropicDiffusion::DoOdeRhs, this);
}
void AnisotropicDiffusion::v_InitialiseFields()
{
    PlasmaSystem::v_InitialiseFields();
    auto volume =
        MR::GetBlockAttributes<double, FieldState::Phys>(m_expansionLists[0]);
    if (m_projectionType == MR::eDiscontinuous)
    {
        auto trace = MR::GetLocTraceBlockAttributes<double, FieldState::Phys>(
            m_expansionLists[0]);
        m_magneticTrace = PhysField("Magnetic trace", trace, 3, 1);
        m_currentTrace  = PhysField("Current trace", trace, 1, 1);
    }
    else
    {
        m_gradient = PhysField("State gradient", volume,
                               m_variables.size() * m_coordim, 1);
        m_flux     = PhysField("Diffusion flux", volume,
                               m_variables.size() * m_coordim, 1);
        auto coeff = MR::GetBlockAttributes<double, FieldState::Coeff>(
            m_expansionLists[0]);
        m_rhs = CoeffField("Explicit mass RHS", coeff, m_variables, 1);
        m_massSolution =
            CoeffField("Explicit mass solution", coeff, m_variables, 1);
    }
}
void AnisotropicDiffusion::v_InitialiseOperators()
{
    PlasmaSystem::v_InitialiseOperators();
    m_volumeFlux = FluxVolDiffOp<double>::Create(
        m_expansionLists[0], m_fieldTemplates, "", DiffusionModel());
    m_volumeFlux->SetMagneticField(B);
    m_volumeFlux->SetElectronEnergyId(ee_idx == InvalidId ? 0 : ee_idx);
    m_volumeFlux->SetSpecies(n_ids, v_ids, e_ids, m_masses, m_charges);
    m_traceExtract = MR::PhysTraceExtractOp<double>::Create(m_expansionLists[0],
                                                            m_fieldTemplates);
    if (m_projectionType == MR::eDiscontinuous)
    {
        m_diffusionOp = SC::DiffusionIPOp<double>::Create(m_expansionLists[0],
                                                          m_fieldTemplates);
        m_traceFlux   = FluxTraceAdvDiffOp<double>::Create(
            m_expansionLists[0], m_fieldTemplates, "", DiffusionModel());
        m_traceFlux->SetSpecies(n_ids, v_ids, e_ids, m_masses, m_charges);
        m_traceFlux->SetElectronEnergyId(ee_idx == InvalidId ? 0 : ee_idx);
        m_traceFlux->SetMagneticTrace(m_magneticTrace);
        m_traceFlux->SetParallelCurrentTrace(m_currentTrace);
        m_traceFlux->SetAdvectionScale(0);
        m_diffusionOp->SetVolumeFluxOp(m_volumeFlux);
        m_diffusionOp->SetTraceFluxOp(m_traceFlux);
    }
    else
    {
        m_derivative = MR::PhysDerivOp<double>::Create(m_expansionLists[0],
                                                       m_fieldTemplates);
        m_weakDerivative =
            MR::IProductWRTDerivBaseOp<FieldState::Coeff, double>::Create(
                m_expansionLists[0], m_fieldTemplates);
        m_weakDerivative->SetScale(1);
        bool robin = false;
        for (const auto &expansion : m_expansionLists)
            for (const auto &condition : expansion->GetBndConditions())
                robin = robin ||
                        condition->GetBoundaryConditionType() == SD::eRobin;
        if (robin)
        {
            m_robin = MR::RobBndCondOp<double>::Create(m_expansionLists[0],
                                                       m_fieldTemplates);
            m_localInnerProduct = MR::IProductWRTBaseOp<double>::Create(
                m_expansionLists[0], m_fieldTemplates);
            m_localInverseMass = MR::MultiplyByElmtInvMassOp<double>::Create(
                m_expansionLists[0], m_fieldTemplates);
            m_stateCoeffs =
                CoeffField("Robin state coefficients",
                           MR::GetBlockAttributes<double, FieldState::Coeff>(
                               m_expansionLists[0]),
                           m_variables, 1);
        }
        m_neumann = MR::NeuBndCondOp<double>::Create(m_expansionLists[0],
                                                     m_fieldTemplates);
        // This solve only inverts the CG mass matrix to evaluate an explicit
        // derivative. Robin stiffness does not belong in the mass matrix.
        class MassLHS : public SC::LHSOp<double>
        {
            class NoRobin : public MR::RobBndCondOp<double>
            {
            public:
                NoRobin(const MR::ExpListSharedPtr &exp,
                        const std::vector<std::string> &names)
                    : MR::RobBndCondOp<double>(exp, names)
                {
                }
                void v_Apply(CoeffField &, CoeffField &) override
                {
                }
            };

        public:
            MassLHS(const MR::ExpListSharedPtr &exp,
                    const std::vector<std::string> &names)
                : SC::LHSOp<double>(exp, names)
            {
                m_robBCOp = std::make_shared<NoRobin>(exp, names);
            }
        };
        auto lhs =
            std::make_shared<MassLHS>(m_expansionLists[0], m_fieldTemplates);
        lhs->SetElmtOp(
            MR::MassOp<double>::Create(m_expansionLists[0], m_fieldTemplates));
        m_massSolver = SC::LinearSolverOp<double>::Create(
            m_expansionLists[0], m_fieldTemplates, "ConjGrad");
        m_massSolver->SetLHS(lhs);
        m_massSolver->SetPrecon(SC::PreconOp<double>::Create(
            m_expansionLists[0], m_fieldTemplates, "Diagonal"));
        m_massSolver->UpdatePrecon();
    }
}
void AnisotropicDiffusion::DoOdeRhs(PhysField &in, PhysField &out,
                                    const double &time, const double &factor)
{
    if (transient_field)
        UpdateMagneticField(time);
    B.SetInterleaveWidth(in);
    if (m_diffusionOp)
    {
        m_traceExtract->Apply(B, m_magneticTrace);
        m_traceExtract->Apply(j_par, m_currentTrace);
        m_diffusionOp->UpdateBndPhys(time);
        m_diffusionOp->SetScale(factor);
        m_diffusionOp->Apply(in, out);
    }
    else
    {
        m_derivative->Apply(in, m_gradient);
        B.SetInterleaveWidth(in);
        m_volumeFlux->SetAppend(false);
        m_volumeFlux->Apply(in, m_gradient, m_flux);
        m_weakDerivative->Apply(m_flux, m_rhs);
        if (m_robin)
        {
            m_localInnerProduct->Apply(in, m_massSolution);
            m_localInverseMass->Apply(m_massSolution, m_stateCoeffs);
            m_robin->Apply(m_stateCoeffs, m_rhs);
        }
        m_math.mul(-1., m_rhs, m_rhs);
        m_neumann->UpdateBndCoeffs(time);
        m_neumann->Apply(m_rhs);
        m_math.zero(m_massSolution);
        m_massSolver->Apply(m_rhs, m_massSolution);
        m_bwdTransOp->Apply(m_massSolution, out);
        m_math.mul(factor, out, out);
    }
    AddSources(in, out, time, factor);
}
} // namespace PENKNIFE
