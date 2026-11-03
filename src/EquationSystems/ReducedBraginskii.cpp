#include "ReducedBraginskii.hpp"
#include "../FluxVolAdvDiff/FluxVolAdvDiffOp.hpp"
#include <LibUtilities/BasicUtils/StringUtils.hpp>

namespace PENKNIFE
{

/// Name of class
static std::string class_name;
std::string ReducedBraginskii::class_name =
    SU::GetEquationSystemFactory().RegisterCreatorFunction(
        "ReducedBraginskii", ReducedBraginskii::create,
        "Solves reduced Braginskii equations");

ReducedBraginskii::ReducedBraginskii(const LU::SessionReaderSharedPtr &session,
                                     const SD::MeshGraphSharedPtr &graph)
    : PlasmaSystem(session, graph)
{
}

/**
 * @brief Initialise the class.
 */
void ReducedBraginskii::v_InitObject(bool DeclareFields)
{
    PlasmaSystem::v_InitObject(DeclareFields);
    m_timeOp->DefineImplicit(&ImplicitHelper::ImplicitTimeInt, m_implHelper);
    m_timeOp->DefineExplicitRhs(&ReducedBraginskii::DoOdeRhs, this);

    this->ee_idx      = m_indfields.size() - this->n_indep_fields;
    m_varConv->ee_idx = this->ee_idx;

    std::string closureName;
    m_session->LoadSolverInfo("Closure", closureName, "Braginskii");
    m_closure = GetClosureFactory().CreateInstance(
        closureName, as<PlasmaSystem>(), m_spacedim);
    m_closure->ee_idx = this->ee_idx;

    std::string diffName;
    m_session->LoadSolverInfo("DiffusionType", diffName, "LDG");
    m_diffusion =
        SolverUtils::GetDiffusionFactory().CreateInstance(diffName, diffName);
    m_diffusion->SetFluxVector(&ReducedBraginskii::GetFluxVectorDiff, this);

    // workaround for bug in DiffusionLDG
    m_difffields = Array<OneD, MR::ExpListSharedPtr>(m_indfields.size());
    friction     = Array<OneD, Array<OneD, NekDouble>>(m_difffields.size());

    for (int f = 0; f < m_difffields.size(); ++f)
    {
        m_difffields[f] = m_indfields[f];
        friction[f]     = Array<OneD, NekDouble>(n_pts, 0.0);
    }

    m_diffusion->InitObject(m_session, m_difffields);

    InitAdvection();

    m_ode.DefineOdeRhs(&ReducedBraginskii::DoOdeRhs, this);

    if (this->particles_enabled)
    {
        std::vector<Sym<REAL>> src_syms;
        std::vector<int> src_components;
        std::vector<Sym<REAL>> out_syms;

        int cnt = 0;
        for (const auto s : this->GetIons())
        {
            const auto &name_s = this->GetSpeciesNames().at(s);
            this->src_fields.emplace_back(
                MemoryManager<MR::DisContField>::AllocateSharedPtr(
                    *std::dynamic_pointer_cast<MR::DisContField>(m_fields[0])));
            src_syms.push_back(Sym<REAL>(name_s + "_SOURCE_DENSITY"));
            src_components.push_back(0);
            ni_src_idx.push_back(cnt++);
            out_syms.push_back(Sym<REAL>(name_s + "_SOURCE_DENSITY"));

            for (int d = 0; d < this->m_spacedim; ++d)
            {
                this->src_fields.emplace_back(
                    MemoryManager<MR::DisContField>::AllocateSharedPtr(
                        *std::dynamic_pointer_cast<MR::DisContField>(
                            m_fields[0])));
                src_syms.push_back(Sym<REAL>(name_s + "_SOURCE_MOMENTUM"));
                src_components.push_back(d);
            }
            vi_src_idx.push_back(cnt);
            cnt += m_spacedim;
            out_syms.push_back(Sym<REAL>(name_s + "_SOURCE_MOMENTUM"));

            this->src_fields.emplace_back(
                MemoryManager<MR::DisContField>::AllocateSharedPtr(
                    *std::dynamic_pointer_cast<MR::DisContField>(m_fields[0])));

            src_syms.push_back(Sym<REAL>(name_s + "_SOURCE_ENERGY"));
            src_components.push_back(0);
            ei_src_idx.push_back(cnt++);
            out_syms.push_back(Sym<REAL>(name_s + "_SOURCE_ENERGY"));
        }
        this->src_fields.emplace_back(
            MemoryManager<MR::DisContField>::AllocateSharedPtr(
                *std::dynamic_pointer_cast<MR::DisContField>(m_fields[0])));
        src_syms.push_back(Sym<REAL>("ELECTRON_SOURCE_ENERGY"));
        src_components.push_back(0);

        this->particle_sys->finish_setup(this->src_fields, src_syms,
                                         src_components);
        this->particle_sys->output_setup(out_syms);
    }
}

/*
 *  Create and initialise Fields.
 */
void ReducedBraginskii::v_InitialiseFields()
{
    // Initialise solution "m_fields" via EquationSystem routine
    PlasmaSystem::v_InitialiseFields();

    // Get interleave width
    auto interleaveWidth =
        Operator<double>::GetDefaultInterleaveWidth(m_session);

    // Create blocks.
    auto bAtr_phys = MultiRegions::GetBlockAttributes<double, FieldState::Phys>(
        m_expansionLists[0], interleaveWidth);
    m_advectionVel = LibUtilities::Field<double, FieldState::Phys>(
        "advectionVel", bAtr_phys, m_coordim, m_npointsZ);

    auto bAtr_phys_trace =
        MultiRegions::GetLocTraceBlockAttributes<double, FieldState::Phys>(
            m_expansionLists[0]);
    m_traceAdvVel = LibUtilities::Field<double, FieldState::Phys>(
        "traceAdvVel", bAtr_phys_trace, m_coordim, m_npointsZ);
}

void ReducedBraginskii::SetUpBoundaryConditions()
{
    // Bohm Boundary Condition Type
    auto bndCondBohmOp =
        BndCondBohmOp<double>::Create(m_expansionLists[0], m_variables);
    m_advDiffusionOp->AddBndCondUpdateOp(bndCondBohmOp);
}

void ReducedBraginskii::v_InitialiseOperators()
{
    PlasmaSystem::v_InitialiseOperators();
    std::string execName =
        MultiRegions::Operator<double>::GetOpExecSpace(m_session);

    // Create advection operators
    m_advDiffusionOp = AdvWeakDGDiffusionIPOp<double>::Create(
        m_expansionLists[0], m_variables);

    std::string riemannMethod = m_session->GetSolverInfo("UpwindType");

    std::string AdvDiffMethod =
        "FluxTraceAdvDiff" + riemannMethod; // add prefix for local version
    auto advDiffusionTraceFluxOp = SolverCore::TraceFluxOp<double>::Create(
        m_expansionLists[0], m_variables, AdvDiffMethod);
    advDiffusionTraceFluxOp->SetViscousTraceForm(ReadViscousTraceForm());
    advDiffusionTraceFluxOp->SetTraceAdvVel(m_traceAdvVel);
    advDiffusionTraceFluxOp->SetTraceB(this->b_unit_trace);
    advDiffusionTraceFluxOp->SetIndices(n_ids, v_ids, e_ids, m_masses,
                                        m_charges);
    m_advDiffusionOp->SetAdvDiffTraceFlux(advDiffusionTraceFluxOp);

    auto advVolFluxOp =
        FluxVolAdvOp<double>::Create(m_expansionLists[0], m_variables);
    advVolFluxOp->SetAdvVel(m_advectionVel);
    advVolFluxOp->SetIndices(n_ids, v_ids, e_ids, m_masses, m_charges);
    m_advDiffusionOp->SetAdvVolFlux(advVolFluxOp);

    m_physTraceExtractOp = PhysTraceExtractOp<double>::Create(
        m_expansionLists[0], m_advectionVel.GetComponentNames());

    auto diffVolFluxOp =
        FluxVolDiffOp<double>::Create(m_expansionLists[0], m_variables);
    diffVolFluxOp->SetMagneticField(this->b_unit);
    diffVolFluxOp->SetElectronEnergyId(ee_idx);
    diffVolFluxOp->SetIndices(n_ids, v_ids, e_ids, m_masses, m_charges);
    m_advDiffusionOp->SetDiffVolFlux(diffVolFluxOp);

    // The average-of-normal-fluxes form needs the viscous flux on its own as
    // well as the total; the CFE operator forms both in one pass.
    if (advDiffusionTraceFluxOp->GetViscousTraceForm() ==
        SolverCore::TraceFluxOp<
            double>::ViscousTraceForm::AverageOfNormalFluxes)
    {
        auto advDiffVolFlux = FluxVolAdvDiffOp<double>::Create(
            m_expansionLists[0], m_variables);
        advDiffVolFlux->SetMagneticField(this->b_unit);
        advDiffVolFlux->SetElectronEnergyId(ee_idx);
        advDiffVolFlux->SetElectronFields(ne, Me, me);
        advDiffVolFlux->SetIndices(n_ids, v_ids, e_ids, m_masses, m_charges);
        m_advDiffusionOp->SetAdvDiffVolFlux(advDiffVolFlux);
    }

    SetUpBoundaryConditions();
}

/**
 * @brief The viscous trace form SOLVERINFO ViscousTraceFlux selects.
 *
 * FluxOfAverages (default) evaluates the viscous trace flux from the
 * averaged state and gradient; AverageOfNormalFluxes averages the two
 * elements' own normal viscous fluxes. Both add the same interior-penalty
 * term, and boundary traces keep the first form in either case.
 */
SolverCore::TraceFluxOp<double>::ViscousTraceForm ReducedBraginskii::
    ReadViscousTraceForm()
{
    using Form = SolverCore::TraceFluxOp<double>::ViscousTraceForm;
    std::string form;
    m_session->LoadSolverInfo("ViscousTraceFlux", form, "FluxOfAverages");
    if (StringUtils::IEquals(form, "AverageOfNormalFluxes"))
    {
        return Form::AverageOfNormalFluxes;
    }
    ASSERTL0(StringUtils::IEquals(form, "FluxOfAverages"),
             "Unknown ViscousTraceFlux '" + form +
                 "'; expected FluxOfAverages or AverageOfNormalFluxes.");
    return Form::FluxOfAverages;
}

/**
 * @brief Initialise the advection object.
 */
void ReducedBraginskii::InitAdvection()
{
    for (unsigned int s = 0; s < this->GetNumSpecies(); ++s)
    {
        if (this->HasMomentum(s))
        {
            int ni_idx = this->GetDensityId(s);
            int vi_idx = this->GetMomentumId(s);
            this->advected_fields.push_back(ni_idx);
            this->advected_fields.push_back(vi_idx);

            if (this->HasEnergy(s))
            {
                int ei_idx = this->GetEnergyId(s);
                this->advected_fields.push_back(ei_idx);
            }
        }
    }

    this->advected_fields.push_back(ee_idx);
    m_advfields = Array<OneD, MR::ExpListSharedPtr>(advected_fields.size());
    for (int a : this->advected_fields)
    {
        m_advfields[a] = m_indfields[a];
    }

    this->adv_vel = Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(
        this->advected_fields.size());

    for (int i = 0; i < this->adv_vel.size(); ++i)
    {
        this->adv_vel[i] = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
        for (int d = 0; d < m_spacedim; ++d)
        {
            this->adv_vel[i][d] = Array<OneD, NekDouble>(this->n_pts, 0.0);
        }
    }

    if (m_indfields[0]->GetTrace())
    {
        auto nTrace = GetTraceNpoints();

        this->trace_vel_norm =
            Array<OneD, Array<OneD, NekDouble>>(adv_vel.size());
        this->trace_b_norm = Array<OneD, NekDouble>(nTrace, 0.0);
        this->adv_vel_trace =
            Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(adv_vel.size());
        for (int i = 0; i < this->trace_vel_norm.size(); ++i)
        {
            this->adv_vel_trace[i] =
                Array<OneD, Array<OneD, NekDouble>>(m_spacedim);

            this->trace_vel_norm[i] = Array<OneD, NekDouble>(nTrace, 0.0);
            for (int d = 0; d < m_spacedim; ++d)
            {
                this->adv_vel_trace[i][d] = Array<OneD, NekDouble>(nTrace, 0.0);
            }
        }
    }

    // Create Riemann solver and set normal velocity
    // callback functions
    m_session->LoadSolverInfo("UpwindType", this->riemann_solver_type,
                              "MultiFieldUpwind");
    this->riemann_solver = SU::GetRiemannSolverFactory().CreateInstance(
        this->riemann_solver_type, m_session);
    auto t = std::dynamic_pointer_cast<PlasmaSolver>(this->riemann_solver);
    this->riemann_solver->SetVector("Vn", &ReducedBraginskii::GetAdvVelNorm,
                                    this);

    // Setup advection object
    m_session->LoadSolverInfo("AdvectionType", this->adv_type, "WeakDG");
    m_advection = SU::GetAdvectionFactory().CreateInstance(this->adv_type,
                                                           this->adv_type);
    m_advection->SetFluxVector(&ReducedBraginskii::GetFluxVector, this);
    m_advection->SetRiemannSolver(this->riemann_solver);
    m_advection->InitObject(m_session, m_indfields);
}

bool ReducedBraginskii::v_PreIntegrate(int step)
{
    if (this->particles_enabled)
    {
        Vmath::Vdiv(this->n_pts, m_indfields[ee_idx]->GetPhys(), 1,
                    m_fields[0]->GetPhys(), 1, Te->UpdatePhys(), 1);
        Vmath::Smul(this->n_pts, 2.0 / 3.0, Te->GetPhys(), 1, Te->UpdatePhys(),
                    1);
    }

    return PlasmaSystem::v_PreIntegrate(step);
}

bool ReducedBraginskii::v_PostIntegrate(int step)
{
    m_fields[0]->FwdTrans(m_fields[0]->GetPhys(), m_fields[0]->UpdateCoeffs());
    m_fields[1]->FwdTrans(m_fields[1]->GetPhys(), m_fields[1]->UpdateCoeffs());
    m_fields[2]->FwdTrans(m_fields[2]->GetPhys(), m_fields[2]->UpdateCoeffs());

    // Writes a step of the particle trajectory.

    return PlasmaSystem::v_PostIntegrate(step);
}

/**
 * @brief Populate rhs array ( @p outarray )
 *
 * @param inarray physical values of all fields
 * @param[out] outarray output array (RHSs of time integration equations)
 */
void ReducedBraginskii::DoOdeRhs(LU::Field<double, FieldState::Phys> &in,
                                 LU::Field<double, FieldState::Phys> &out,
                                 [[maybe_unused]] const double &time,
                                 const double &dt)
{
    m_quasineutralityOp->Apply(in, j_par, ne, Me);

    // Calculate E
    ComputeE();
    // Calculate ExB, parallel and diamagnetic velocities
    CalcVelocities(in);

    m_advDiffusionOp->SetScale(dt);
    m_advDiffusionOp->Apply(in, out);

    if (this->particles_enabled)
    {
        DoParticles(in, out);
    }

    // Add forcing terms
    for (unsigned int i = 0; i < m_forcing.size(); ++i)
    {
        m_forcing[i]->SetAppend(true);
        m_forcing[i]->Apply(in, out, time, dt);
    }
}

/**
 * @brief Compute the advection terms for the right-hand side
 */
void ReducedBraginskii::DoParticles(LU::Field<double, FieldState::Phys> &in,
                                    LU::Field<double, FieldState::Phys> &out)
{
    m_math.add(in, src_fields, out);
}

/**
 * @brief Compute the gradient of phi for evaluation at the particle positions.
 */
void ReducedBraginskii::ComputeE()
{
    Vmath::Neg(this->n_pts, this->E[0]->UpdatePhys(), 1);
    Vmath::Neg(this->n_pts, this->E[1]->UpdatePhys(), 1);
    Vmath::Neg(this->n_pts, this->E[2]->UpdatePhys(), 1);

    this->E[0]->FwdTrans(this->E[0]->GetPhys(), this->E[0]->UpdateCoeffs());
    this->E[1]->FwdTrans(this->E[1]->GetPhys(), this->E[1]->UpdateCoeffs());
    this->E[2]->FwdTrans(this->E[2]->GetPhys(), this->E[2]->UpdateCoeffs());
}

void ReducedBraginskii::CalcVelocities(LU::Field<double, FieldState::Phys> &in)
{
    for (int f = 0; f < this->adv_vel.size(); ++f)
    {
        for (int d = 0; d < this->adv_vel[f].size(); ++d)
        {
            Vmath::Zero(this->n_pts, this->adv_vel[f][d], 1);
        }
    }
    const Array<OneD, NekDouble> &ne = m_fields[0]->GetPhys();
    // Zero Electron velocity
    Array<OneD, NekDouble> &j_i = m_fields[1]->UpdatePhys();
    Vmath::Zero(this->n_pts, j_i, 1);
    for (const auto s : GetIons())
    {
        const double mass_s   = GetMasses().at(s);
        const double charge_s = GetCharges().at(s);
        int ni_idx            = GetDensityId(s);
        int vi_idx            = GetMomentumId(s);
        int ei_idx            = GetEnergyId(s);

        for (int p = 0; p < this->n_pts; ++p)
        {
            j_i[p] += charge_s * inarray[vi_idx][p] / mass_s;
            double v_i_par = inarray[vi_idx][p] / (mass_s * inarray[ni_idx][p]);
            for (int d = 0; d < m_spacedim; ++d)
            {
                this->adv_vel[ni_idx][d][p] = v_i_par * this->b_unit[d][p];
                this->adv_vel[vi_idx][d][p] = v_i_par * this->b_unit[d][p];
                this->adv_vel[ei_idx][d][p] = v_i_par * this->b_unit[d][p];
            }
        }
    }
    Array<OneD, NekDouble> j_par(this->n_pts, 0.0);
    // TODO calculate conductivity
    double sigma = 0;
    for (int p = 0; p < this->n_pts; ++p)
    {
        for (int d = 0; d < m_spacedim; ++d)
        {
            j_par[p] += sigma * this->E[d]->GetPhys()[p] * this->b_unit[d][p];
        }
    }
    for (int p = 0; p < this->n_pts; ++p)
    {
        for (int d = 0; d < m_spacedim; ++d)
        {
            this->adv_vel[ee_idx][d][p] =
                this->b_unit[d][p] * (j_i[p] - j_par[p]) / ne[p];
        }
    }

    for (const auto s : GetNeutrals())
    {
        const double mass_s = GetMasses().at(s);
        int nn_idx          = GetDensityId(s);
        int vn_idx          = GetMomentumId(s);
        int en_idx          = GetEnergyId(s);

        for (int p = 0; p < this->n_pts; ++p)
        {
            double v_i_par = inarray[vn_idx][p] / (mass_s * inarray[nn_idx][p]);
            for (int d = 0; d < m_spacedim; ++d)
            {
                this->adv_vel[nn_idx][d][p] = v_i_par * this->b_unit[d][p];
                this->adv_vel[vn_idx][d][p] = v_i_par * this->b_unit[d][p];
                this->adv_vel[en_idx][d][p] = v_i_par * this->b_unit[d][p];
            }
        }
    }

    m_physTraceExtractOp->Apply(m_advectionVel, m_traceAdvVel);
}

void ReducedBraginskii::CalcNeutralSources_nvp(
    const double m, int ee_idx, int ni_idx, int vi_idx, int ei_idx, int nn_idx,
    int vn_idx, int en_idx, const Array<OneD, Array<OneD, NekDouble>> &inarray,
    const Array<OneD, NekDouble> &ne,
    Array<OneD, Array<OneD, NekDouble>> &outarray, Array<OneD, NekDouble> &See)
{
    const Array<OneD, NekDouble> &ee = inarray[ee_idx];
    const Array<OneD, NekDouble> &nn = inarray[nn_idx];
    const Array<OneD, NekDouble> &vn = inarray[vn_idx];
    const Array<OneD, NekDouble> &en = inarray[en_idx];
    const Array<OneD, NekDouble> &ni = inarray[ni_idx];
    const Array<OneD, NekDouble> &vi = inarray[vi_idx];
    const Array<OneD, NekDouble> &ei = inarray[ei_idx];

    for (int p = 0; p < this->n_pts; ++p)
    {
        double exponent = 13.6 / ee[p];
        double krec     = 0.7e-19 * std::sqrt(exponent);
        double kIZ      = (2e-13 / (6 + 1.0 / exponent)) *
                          std::sqrt(1.0 / exponent) * std::exp(-exponent);
        double kCX      = 3.2e-15 * std::sqrt(ei[p] / 0.026);

        double SN = -kIZ * ne[p] * nn[p] + krec * ne[p] * ni[p];

        double SGN = -kIZ * ne[p] * vn[p] + krec * ne[p] * vi[p] +
                     kCX * nn[p] * vi[p] - kCX * ni[p] * vn[p];

        double vdiffsq = (vi[p] / ni[p] - vn[p] / nn[p]) *
                         (vi[p] / ni[p] - vn[p] / nn[p]) / (m * m);
        double SP = -kIZ * ne[p] * en[p] + krec * ne[p] * ei[p] +
                    kCX * nn[p] * ei[p] - kCX * ni[p] * en[p] +
                    (m * ni[p] / 3) * (krec * ne[p] + kCX * nn[p]) * vdiffsq;

        double SPI = ni[p] * nn[p] * (kIZ + kCX) *
                     ((en[p] / nn[p] - ei[p] / ni[p]) + (m / 3) * vdiffsq);

        double Wiz  = 1.0; // TODO get from AMJUEL
        double Wrec = 1.0;

        See[p] -= ne[p] * ee[p] * (kIZ * nn[p] - krec * ni[p]) +
                  ne[p] * (2.0 / 3.0) * (Wiz * nn[p] + Wrec * ni[p]);

        outarray[nn_idx][p] += SN;
        outarray[ni_idx][p] -= SN;
        outarray[vn_idx][p] += m * SGN;
        outarray[vi_idx][p] -= m * SGN + SN * vi[p] / ni[p];
        outarray[en_idx][p] += SP;
        outarray[ei_idx][p] += SPI;
    }
}

void ReducedBraginskii::CalcNeutralSources_nv(
    const double m, int ee_idx, int ni_idx, int vi_idx, int nn_idx, int vn_idx,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    const Array<OneD, NekDouble> &ne,
    Array<OneD, Array<OneD, NekDouble>> &outarray, Array<OneD, NekDouble> &See)
{
    const Array<OneD, NekDouble> &ee = inarray[ee_idx];
    const Array<OneD, NekDouble> &nn = inarray[nn_idx];
    const Array<OneD, NekDouble> &vn = inarray[vn_idx];
    const Array<OneD, NekDouble> &ni = inarray[ni_idx];
    const Array<OneD, NekDouble> &vi = inarray[vi_idx];

    for (int p = 0; p < this->n_pts; ++p)
    {
        double exponent = 13.6 / ee[p];
        double krec     = 0.7e-19 * std::sqrt(exponent);
        double kIZ      = (2e-13 / (6 + 1.0 / exponent)) *
                          std::sqrt(1.0 / exponent) * std::exp(-exponent);
        double kCX      = 3.2e-15 * std::sqrt(2.0 / 0.026);
        // TODO use background temp instead of 2.0

        double SN = -kIZ * ne[p] * nn[p] + krec * ne[p] * ni[p];

        double SGN = -kIZ * ne[p] * vn[p] + krec * ne[p] * vi[p] +
                     kCX * nn[p] * vi[p] - kCX * ni[p] * vn[p];

        double Wiz  = 1.0; // TODO get from AMJUEL
        double Wrec = 1.0;

        See[p] -= ne[p] * ee[p] * (kIZ * nn[p] - krec * ni[p]) +
                  ne[p] * (2.0 / 3.0) * (Wiz * nn[p] + Wrec * ni[p]);

        outarray[nn_idx][p] += SN;
        outarray[ni_idx][p] -= SN;
        outarray[vn_idx][p] += m * SGN;
        outarray[vi_idx][p] -= m * SGN + SN * vi[p] / ni[p];
    }
}

void ReducedBraginskii::CalcNeutralSources_nv(
    const double m, int ee_idx, int ni_idx, int vi_idx, int ei_idx, int nn_idx,
    int vn_idx, const Array<OneD, Array<OneD, NekDouble>> &inarray,
    const Array<OneD, NekDouble> &ne,
    Array<OneD, Array<OneD, NekDouble>> &outarray, Array<OneD, NekDouble> &See)
{
    const Array<OneD, NekDouble> &ee = inarray[ee_idx];
    const Array<OneD, NekDouble> &nn = inarray[nn_idx];
    const Array<OneD, NekDouble> &vn = inarray[vn_idx];
    const Array<OneD, NekDouble> &ni = inarray[ni_idx];
    const Array<OneD, NekDouble> &vi = inarray[vi_idx];
    const Array<OneD, NekDouble> &ei = inarray[ei_idx];

    for (int p = 0; p < this->n_pts; ++p)
    {
        double exponent = 13.6 / ee[p];
        double krec     = 0.7e-19 * std::sqrt(exponent);
        double kIZ      = (2e-13 / (6 + 1.0 / exponent)) *
                          std::sqrt(1.0 / exponent) * std::exp(-exponent);
        double kCX      = 3.2e-15 * std::sqrt(ei[p] / 0.026);

        double SN = -kIZ * ne[p] * nn[p] + krec * ne[p] * ni[p];

        double SGN = -kIZ * ne[p] * vn[p] + krec * ne[p] * vi[p] +
                     kCX * nn[p] * vi[p] - kCX * ni[p] * vn[p];

        double vdiffsq = (vi[p] / ni[p] - vn[p] / nn[p]) *
                         (vi[p] / ni[p] - vn[p] / nn[p]) / (m * m);

        double SPI = ni[p] * nn[p] * (kIZ + kCX) * (m / 3) * vdiffsq;

        double Wiz  = 1.0; // TODO get from AMJUEL
        double Wrec = 1.0;

        See[p] -= ne[p] * ee[p] * (kIZ * nn[p] - krec * ni[p]) +
                  ne[p] * (2.0 / 3.0) * (Wiz * nn[p] + Wrec * ni[p]);

        outarray[nn_idx][p] += SN;
        outarray[ni_idx][p] -= SN;
        outarray[vn_idx][p] += m * SGN;
        outarray[vi_idx][p] -= m * SGN + SN * vi[p] / ni[p];
        outarray[ei_idx][p] += SPI;
    }
}

void ReducedBraginskii::CalcNeutralSources_np(
    const double m, int ee_idx, int ni_idx, int ei_idx, int nn_idx, int en_idx,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    const Array<OneD, NekDouble> &ne,
    Array<OneD, Array<OneD, NekDouble>> &outarray, Array<OneD, NekDouble> &See)
{
    const Array<OneD, NekDouble> &ee = inarray[ee_idx];
    const Array<OneD, NekDouble> &nn = inarray[nn_idx];
    const Array<OneD, NekDouble> &en = inarray[en_idx];
    const Array<OneD, NekDouble> &ni = inarray[ni_idx];
    const Array<OneD, NekDouble> &ei = inarray[ei_idx];

    for (int p = 0; p < this->n_pts; ++p)
    {
        double exponent = 13.6 / ee[p];
        double krec     = 0.7e-19 * std::sqrt(exponent);
        double kIZ      = (2e-13 / (6 + 1.0 / exponent)) *
                          std::sqrt(1.0 / exponent) * std::exp(-exponent);
        double kCX      = 3.2e-15 * std::sqrt(ei[p] / 0.026);

        double SN = -kIZ * ne[p] * nn[p] + krec * ne[p] * ni[p];

        double SP = -kIZ * ne[p] * en[p] + krec * ne[p] * ei[p] +
                    kCX * nn[p] * ei[p] - kCX * ni[p] * en[p];

        double SPI =
            ni[p] * nn[p] * (kIZ + kCX) * (en[p] / nn[p] - ei[p] / ni[p]);

        double Wiz  = 1.0; // TODO get from AMJUEL
        double Wrec = 1.0;

        See[p] -= ne[p] * ee[p] * (kIZ * nn[p] - krec * ni[p]) +
                  ne[p] * (2.0 / 3.0) * (Wiz * nn[p] + Wrec * ni[p]);

        outarray[nn_idx][p] += SN;
        outarray[ni_idx][p] -= SN;
        outarray[en_idx][p] += SP;
        outarray[ei_idx][p] += SPI;
    }
}

void ReducedBraginskii::CalcNeutralSources_n(
    const double m, int ee_idx, int ni_idx, int nn_idx,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    const Array<OneD, NekDouble> &ne,
    Array<OneD, Array<OneD, NekDouble>> &outarray, Array<OneD, NekDouble> &See)
{
    const Array<OneD, NekDouble> &ee = inarray[ee_idx];
    const Array<OneD, NekDouble> &nn = inarray[nn_idx];
    const Array<OneD, NekDouble> &ni = inarray[ni_idx];

    for (int p = 0; p < this->n_pts; ++p)
    {
        double exponent = 13.6 / ee[p];
        double krec     = 0.7e-19 * std::sqrt(exponent);
        double kIZ      = (2e-13 / (6 + 1.0 / exponent)) *
                          std::sqrt(1.0 / exponent) * std::exp(-exponent);

        double SN = -kIZ * ne[p] * nn[p] + krec * ne[p] * ni[p];

        double Wiz  = 1.0; // TODO get from AMJUEL
        double Wrec = 1.0;

        See[p] -= ne[p] * ee[p] * (kIZ * nn[p] - krec * ni[p]) +
                  ne[p] * (2.0 / 3.0) * (Wiz * nn[p] + Wrec * ni[p]);

        outarray[nn_idx][p] += SN;
        outarray[ni_idx][p] -= SN;
    }
}

void ReducedBraginskii::AddNeutralSources(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    auto ne = m_fields[0]->GetPhys();
    Array<OneD, NekDouble> tmp(this->n_pts, 0.0);

    for (const auto s : this->GetNeutrals())
    {
        const double mass_s = this->GetMasses().at(s);
        int nn_idx          = this->GetDensityId(s);
        int ni_idx          = this->GetDensityId(this->GetNeutralIonId(s));

        if (this->HasEnergy(s))
        {
            int en_idx = this->GetEnergyId(s);
            int ei_idx = this->GetEnergyId(this->GetNeutralIonId(s));

            if (this->HasMomentum(s))
            {
                // Neutral n,v,p; Ion n,v,p
                int vn_idx = this->GetMomentumId(s);
                int vi_idx = this->GetMomentumId(this->GetNeutralIonId(s));
                CalcNeutralSources_nvp(mass_s, ee_idx, ni_idx, vi_idx, ei_idx,
                                       nn_idx, vn_idx, en_idx, inarray, ne,
                                       outarray, tmp);
            }
            else
            {
                // Neutral n,p; Ion n,v,p; Ion n,p;
                CalcNeutralSources_np(mass_s, ee_idx, ni_idx, ei_idx, nn_idx,
                                      en_idx, inarray, ne, outarray, tmp);
            }
        }

        else if (this->HasMomentum(s))
        {
            int vn_idx = this->GetMomentumId(s);
            int vi_idx = this->GetMomentumId(this->GetNeutralIonId(s));
            if (this->HasEnergy(this->GetNeutralIonId(s)))
            {
                // Neutral n,v; Ion n,v,p
                int ei_idx = this->GetEnergyId(this->GetNeutralIonId(s));
                CalcNeutralSources_nv(mass_s, ee_idx, ni_idx, vi_idx, ei_idx,
                                      nn_idx, vn_idx, inarray, ne, outarray,
                                      tmp);
            }
            else
            {
                // Neutral n,v; Ion n,v
                CalcNeutralSources_nv(mass_s, ee_idx, ni_idx, vi_idx, nn_idx,
                                      vn_idx, inarray, ne, outarray, tmp);
            }
        }
        else
        {
            // Neutral n; Ion n; Ion n,v; Ion n,v,p
            CalcNeutralSources_n(mass_s, ee_idx, ni_idx, nn_idx, inarray, ne,
                                 outarray, tmp);
        }
    }
    Vmath::Vadd(this->n_pts, outarray[ee_idx], 1, tmp, 1, outarray[ee_idx], 1);
}

/**
 * @brief After reading ICs, calculate phi and grad(phi)
 */
void ReducedBraginskii::v_SetInitialConditions(NekDouble init_time,
                                               bool dump_ICs, const int domain)
{
    PlasmaSystem::v_SetInitialConditions(init_time, dump_ICs, domain);
    Checkpoint_Output(0);
}

void ReducedBraginskii::v_ExtraFldOutput(
    std::vector<Array<OneD, NekDouble>> &fieldcoeffs,
    std::vector<std::string> &variables)
{
    PlasmaSystem::v_ExtraFldOutput(fieldcoeffs, variables);

    if (this->particles_enabled)
    {
        int cnt = 0;
        for (const auto k : this->GetIons())
        {
            const auto &name_k = this->GetSpeciesNames().at(k);
            variables.emplace_back(name_k + "_SOURCE_DENSITY");
            Array<OneD, NekDouble> SrcFwd1(this->n_coeffs);
            m_fields[0]->FwdTransLocalElmt(this->src_fields[cnt++]->GetPhys(),
                                           SrcFwd1);
            fieldcoeffs.emplace_back(SrcFwd1);

            for (int d = 0; d < this->m_spacedim; ++d)
            {
                variables.emplace_back(name_k + "_SOURCE_MOMENTUM" +
                                       std::to_string(d));
                Array<OneD, NekDouble> SrcFwd(this->n_coeffs);
                m_fields[0]->FwdTransLocalElmt(
                    this->src_fields[cnt++]->GetPhys(), SrcFwd);
                fieldcoeffs.emplace_back(SrcFwd);
            }

            variables.emplace_back(name_k + "_SOURCE_ENERGY");
            Array<OneD, NekDouble> SrcFwd2(this->n_coeffs);
            m_fields[0]->FwdTransLocalElmt(this->src_fields[cnt++]->GetPhys(),
                                           SrcFwd2);
            fieldcoeffs.emplace_back(SrcFwd2);
        }
    }
}
} // namespace PENKNIFE