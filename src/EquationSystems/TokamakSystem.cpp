#include "TokamakSystem.hpp"

#include <LibUtilities/BasicUtils/Vmath.hpp>

#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <boost/core/ignore_unused.hpp>

namespace NESO::Solvers::tokamak
{

/// Name of class
static std::string class_name;
std::string TokamakSystem::class_name =
    SU::GetEquationSystemFactory().RegisterCreatorFunction(
        "Tokamak", TokamakSystem::create,
        "Tokamak equation system. Runs in either a 2D or 3D "
        "domain.");
/**
 * @brief Creates an instance of this class.
 */
static SU::EquationSystemSharedPtr create(
    const LU::SessionReaderSharedPtr &session,
    const SD::MeshGraphSharedPtr &graph)
{
    SU::EquationSystemSharedPtr p =
        MemoryManager<TokamakSystem>::AllocateSharedPtr(session, graph);
    p->InitObject();
    return p;
}

TokamakSystem::TokamakSystem(const LU::SessionReaderSharedPtr &session,
                             const SD::MeshGraphSharedPtr &graph)
    : TimeEvoEqnSysBase<SU::UnsteadySystem, ParticleSystem>(session, graph)
{
}

std::shared_ptr<ParticleSystem> TokamakSystem::GetParticleSystem()
{
    return this->particle_sys;
}

/**
 * @brief Load all required session parameters into member variables.
 */
void TokamakSystem::load_params()
{
    TimeEvoEqnSysBase<SU::UnsteadySystem, ParticleSystem>::load_params();

    m_session->LoadParameter("mesh_length", this->mesh_length, 1.);
    m_session->LoadParameter("Nnorm", this->Nnorm, 1e18);
    m_session->LoadParameter("Tnorm", this->Tnorm, 100.);
    m_session->LoadParameter("Bnorm", this->Bnorm, 1.);
    this->me = 1. / 1836;

    this->cs = std::sqrt(constants::qeomp *
                         this->Tnorm); // Reference sound speed [m/s]
    this->omega_c =
        constants::qeomp * this->Bnorm;     // Ion cyclotron frequency [1/s]
    this->rho_s = this->cs / this->omega_c; // Length scale [m]

    std::string transient_field_str;
    m_session->LoadSolverInfo("MagneticFieldEvolution", transient_field_str,
                              "Static");
    this->transient_field = (transient_field_str == "Transient");

    // Particle-related parameters
    m_session->LoadParameter("particle_output_freq", particle_output_freq, 0);
    m_session->LoadParameter("num_particle_steps_per_fluid_step",
                             this->num_part_substeps, 1);
    this->part_timestep = m_timestep / this->num_part_substeps;
}

/**
 * @brief Perform projection into correct polynomial space.
 *
 * @details This routine projects the @p in_arr input and ensures the @p
 * out_arr output lives in the correct space. Since we are hard-coding DG,
 * this corresponds to a simple copy from in to out, since no elemental
 * connectivity is required and the output of the RHS function is
 * polynomial.
 *
 * @param in_arr Unprojected values
 * @param[out] out_arr Projected values
 * @param time Current simulation time
 *
 */
void TokamakSystem::DoOdeProjection(
    const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
    Array<OneD, Array<OneD, NekDouble>> &out_arr, const NekDouble time)
{
    int i;
    int num_vars = in_arr.size();
    int npoints  = GetNpoints();
    SetBoundaryConditions(time);

    switch (m_projectionType)
    {
        case MultiRegions::eDiscontinuous:
        {
            // Just copy over array
            if (in_arr != out_arr)
            {
                int npoints = GetNpoints();

                for (i = 0; i < num_vars; ++i)
                {
                    Vmath::Vcopy(npoints, in_arr[i], 1, out_arr[i], 1);
                }
            }
            break;
        }
        case MultiRegions::eGalerkin:
        {
            Array<OneD, NekDouble> coeffs(m_fields[0]->GetNcoeffs());

            for (i = 0; i < num_vars; ++i)
            {
                m_indfields[i]->FwdTrans(in_arr[i], coeffs);
                m_indfields[i]->BwdTrans(coeffs, out_arr[i]);
            }
            break;
        }
        default:
        {
            ASSERTL0(false, "Unknown projection scheme");
            break;
        }
    }
}

void TokamakSystem::v_ExtraFldOutput(
    std::vector<Array<OneD, NekDouble>> &fieldcoeffs,
    std::vector<std::string> &variables)
{
    bool extraFields;
    m_session->MatchSolverInfo("OutputSpeciesFields", "True", extraFields,
                               true);
    if (extraFields)
    {
        const int nPhys   = m_fields[0]->GetNpoints();
        const int nCoeffs = m_fields[0]->GetNcoeffs();

        for (const auto &[k, v] : this->GetSpecies())
        {
            for (const auto &[f, fi] : v.fields)
            {
                variables.push_back(m_session->GetVariable(f) + "_" + v.name);
                Array<OneD, NekDouble> Fwd(nCoeffs);
                this->m_indfields[fi]->FwdTransLocalElmt(
                    this->m_indfields[fi]->GetPhys(), Fwd);
                fieldcoeffs.push_back(Fwd);
            }
        }
        if (Te)
        {
            variables.push_back("Te");
            Array<OneD, NekDouble> Fwd(nCoeffs);
            this->Te->FwdTransLocalElmt(this->Te->GetPhys(), Fwd);
            fieldcoeffs.push_back(Fwd);
        }
    }

    m_session->MatchSolverInfo("OutputEMFields", "True", extraFields, true);
    if (extraFields)
    {
        const int nPhys   = m_fields[0]->GetNpoints();
        const int nCoeffs = m_fields[0]->GetNcoeffs();

        variables.push_back("Bx");
        Array<OneD, NekDouble> BxFwd(nCoeffs);
        m_fields[0]->FwdTransLocalElmt(this->B[0]->GetPhys(), BxFwd);
        fieldcoeffs.push_back(BxFwd);

        variables.push_back("By");
        Array<OneD, NekDouble> ByFwd(nCoeffs);
        m_fields[0]->FwdTransLocalElmt(this->B[1]->GetPhys(), ByFwd);
        fieldcoeffs.push_back(ByFwd);

        variables.push_back("Bz");
        Array<OneD, NekDouble> BzFwd(nCoeffs);
        m_fields[0]->FwdTransLocalElmt(this->B[2]->GetPhys(), BzFwd);
        fieldcoeffs.push_back(BzFwd);

        variables.push_back("Ex");
        Array<OneD, NekDouble> ExFwd(nCoeffs);
        m_fields[0]->FwdTransLocalElmt(this->E[0]->GetPhys(), ExFwd);
        fieldcoeffs.push_back(ExFwd);

        variables.push_back("Ey");
        Array<OneD, NekDouble> EyFwd(nCoeffs);
        m_fields[0]->FwdTransLocalElmt(this->E[1]->GetPhys(), EyFwd);
        fieldcoeffs.push_back(EyFwd);

        variables.push_back("Ez");
        Array<OneD, NekDouble> EzFwd(nCoeffs);
        m_fields[0]->FwdTransLocalElmt(this->E[2]->GetPhys(), EzFwd);
        fieldcoeffs.push_back(EzFwd);
    }
    m_session->MatchSolverInfo("OutputPartitions", "True", extraFields, false);
    if (extraFields)
    {
        const int nPhys   = m_fields[0]->GetNpoints();
        const int nCoeffs = m_fields[0]->GetNcoeffs();
        variables.push_back("Rank");
        Array<OneD, NekDouble> Rank(nPhys,
                                    this->m_session->GetComm()->GetRank());
        Array<OneD, NekDouble> RankFwd(nCoeffs);
        m_fields[0]->FwdTransLocalElmt(Rank, RankFwd);
        fieldcoeffs.push_back(RankFwd);
    }
}

/**
 * @brief Post-construction class initialisation.
 *
 * @param create_field if true, create a new field object and add it to
 * m_fields. Optional, defaults to true.
 */
void TokamakSystem::v_InitObject(bool create_field)
{
    TimeEvoEqnSysBase::v_InitObject(create_field);

    m_domains          = m_graph->GetDomain();
    m_dom_to_offset[0] = 0;
    int size           = 0;

    for (size_t d = 0; d < m_domains.size(); ++d)
    {
        std::vector<int> global;
        std::vector<int> local;
        std::vector<int> phys_offset;
        std::vector<int> coeff_offset;

        for (auto &compIter : m_domains[d])
        {
            for (auto &x : compIter.second->m_geomVec)
            {

                global.push_back(x->GetGlobalID());

                local.push_back(
                    m_fields[0]->GetElmtToExpId().at(global.back()));

                phys_offset.push_back(
                    m_fields[0]->GetPhys_Offset(local.back()));

                size += m_fields[0]->GetExp(local.back())->GetTotPoints();
            }
        }
        m_dom_to_offset[d + 1] = size;
    }

    // Store FieldSharedPtr casts of fields in a map, indexed by name

    this->E = Array<OneD, MR::DisContFieldSharedPtr>(3);
    this->B = Array<OneD, MR::DisContFieldSharedPtr>(3);
    for (int d = 0; d < 3; ++d)
    {
        this->E[d] = MemoryManager<MR::DisContField>::AllocateSharedPtr(
            *std::dynamic_pointer_cast<MR::DisContField>(m_fields[0]));
        this->B[d] = MemoryManager<MR::DisContField>::AllocateSharedPtr(
            *std::dynamic_pointer_cast<MR::DisContField>(m_fields[0]));
        this->B[d]->GetTrace();
    }
    this->mag_field = std::make_shared<MagneticField>(
        m_session, as<TokamakSystem>(), this->B, m_spacedim);
    this->b_unit = this->mag_field->b_unit;
    this->mag_B  = this->mag_field->mag_B;
    this->mag_field->Update(0);

    this->ve = Array<OneD, MR::DisContFieldSharedPtr>(3);

    int nfields = 0;
    std::map<std::string, int> ions_by_name;
    int s_idx = 0;
    int f_idx = 0;
    for (auto &[k, v] : this->neso_config->get_species())
    {
        double mass, charge;
        this->neso_config->load_species_parameter(k, "Mass", mass, 1.0);
        this->neso_config->load_species_parameter(k, "Charge", charge, 0.0);

        Species spec{mass, k};
        for (const auto &str : std::get<2>(v))
        {
            spec.fields[field_to_index[str]] = f_idx++;
        }
        if (std::get<3>(v).has_value())
        {
            int ion = ions_by_name[*std::get<3>(v)];
            Neutral neutral{spec, ion};
            m_neutrals[s_idx] = neutral;
        }
        else
        {
            Ion ion{spec, charge};
            m_ions[s_idx]   = ion;
            ions_by_name[k] = s_idx;
        }
        m_species[s_idx++] = spec;
    }

    this->n_species = m_species.size();

    m_indfields = Array<OneD, MR::ExpListSharedPtr>(n_indep_fields + f_idx);

    for (const auto &[k, v] : this->GetSpecies())
    {
        std::string name = v.name;
        for (const auto &[f, fi] : v.fields)
        {
            if (m_projectionType == MR::eGalerkin ||
                m_projectionType == MR::eMixed_CG_Discontinuous)
            {
                m_indfields[fi] =
                    MemoryManager<MR::ContField>::AllocateSharedPtr(
                        *std::dynamic_pointer_cast<MR::ContField>(m_fields[0]),
                        m_graph, m_session->GetVariable(f), true,
                        m_checkIfSystemSingular[f]);
            }
            else
            {
                m_indfields[fi] =
                    MemoryManager<MR::DisContField>::AllocateSharedPtr(
                        *std::dynamic_pointer_cast<MR::DisContField>(
                            m_fields[0]),
                        m_graph, m_session->GetVariable(f));
            }
        }
    }

    for (int i = 0; i < n_indep_fields; ++i)
    {
        m_indfields[f_idx + i] = m_fields[m_fields.size() - n_indep_fields + i];
    }
    // Ensure DG is setup
    for (int i = 1; i < m_indfields.size(); ++i)
    {
        m_indfields[i]->GetTrace();
    }

    // Bind projection function for time integration object
    m_ode.DefineProjection(&TokamakSystem::DoOdeProjection, this);
    if (m_projectionType == MR::eDiscontinuous)
    {
        // Turn off forward-transform of initial conditions.
        m_homoInitialFwd = false;
        m_implHelper     = std::make_shared<ImplicitHelper>(
            m_session, m_indfields, m_ode, m_indfields.size());
        m_implHelper->InitialiseNonlinSysSolver();
        m_ode.DefineImplicitSolve(&ImplicitHelper::ImplicitTimeInt,
                                  m_implHelper);
    }

    m_varConv = MemoryManager<VariableConverter>::AllocateSharedPtr(
        as<TokamakSystem>(), m_spacedim);

    // Forcing terms
    m_forcing = SU::Forcing::Load(m_session, shared_from_this(), m_indfields,
                                  m_indfields.size());

    m_bndConds = MemoryManager<TokamakBoundaryConditions>::AllocateSharedPtr();

    m_bndConds->Initialize(m_session, as<TokamakSystem>(), m_indfields, B, E,
                           m_spacedim);
}

/**
 * @brief Initialises the time integration scheme (as specified in the
 * session file), and perform the time integration.
 */
void TokamakSystem::v_DoSolve()
{
    ASSERTL0(m_intScheme != nullptr, "No time integration scheme.");

    int i          = 1;
    int nvariables = 0;
    int nfields    = m_indfields.size();
    if (m_intVariables.empty())
    {
        for (i = 0; i < nfields; ++i)
        {
            m_intVariables.push_back(i);
        }
        nvariables = nfields;
    }
    else
    {
        nvariables = m_intVariables.size();
    }

    // Set up wrapper to fields data storage.
    Array<OneD, Array<OneD, NekDouble>> fields(nvariables);

    // Order storage to list time-integrated fields first.
    // @TODO: Refactor to take coeffs (FwdTrans) if boolean flag (in constructor
    // function) says to.
    for (i = 0; i < nvariables; ++i)
    {
        fields[i] = m_indfields[m_intVariables[i]]->UpdatePhys();
        m_indfields[m_intVariables[i]]->SetPhysState(false);
    }

    // @TODO: Virtual function that allows to transform the field space, embed
    // the MultiplyMassMatrix in here.
    // @TODO: Specify what the fields variables are physical or coefficient,
    // boolean in UnsteadySystem class...

    v_ALEPreMultiplyMass(fields);

    // Initialise time integration scheme.
    m_intScheme->InitializeScheme(m_timestep, fields, m_time, m_ode);

    LibUtilities::Timer timer;
    bool doCheckTime      = false;
    int step              = m_initialStep;
    int stepCounter       = 0;
    NekDouble intTime     = 0.0;
    NekDouble cpuTime     = 0.0;
    NekDouble cpuPrevious = 0.0;
    NekDouble elapsed     = 0.0;

    m_lastCheckTime = 0.0;

    Array<OneD, int> abortFlags(2, 0);
    std::string abortFile = "abort";
    if (m_session->DefinesSolverInfo("CheckAbortFile"))
    {
        abortFile = m_session->GetSolverInfo("CheckAbortFile");
    }

    NekDouble tmp_cflSafetyFactor = m_cflSafetyFactor;
    while ((step < m_steps || m_time < m_fintime - NekConstants::kNekZeroTol) &&
           abortFlags[1] == 0)
    {
        if (m_CFLGrowth > 1.0 && m_cflSafetyFactor < m_CFLEnd)
        {
            tmp_cflSafetyFactor =
                std::min(m_CFLEnd, m_CFLGrowth * tmp_cflSafetyFactor);
        }

        // Frozen preconditioner checks.
        if (!m_comm->IsParallelInTime())
        {
            if (v_UpdateTimeStepCheck())
            {
                m_cflSafetyFactor = tmp_cflSafetyFactor;

                if (m_cflSafetyFactor)
                {
                    m_timestep = GetTimeStep(fields);
                }

                // Ensure that the final timestep finishes at the final
                // time, or at a prescribed IO_CheckTime.
                if (m_time + m_timestep > m_fintime && m_fintime > 0.0)
                {
                    m_timestep = m_fintime - m_time;
                }
                else if (m_checktime &&
                         m_time + m_timestep - m_lastCheckTime >= m_checktime)
                {
                    m_lastCheckTime += m_checktime;
                    m_timestep  = m_lastCheckTime - m_time;
                    doCheckTime = true;
                }
            }
        }

        if (m_TimeIncrementFactor > 1.0)
        {
            NekDouble timeincrementFactor = m_TimeIncrementFactor;
            m_timestep *= timeincrementFactor;

            if (m_time + m_timestep > m_fintime && m_fintime > 0.0)
            {
                m_timestep = m_fintime - m_time;
            }
        }

        // Perform any solver-specific pre-integration steps.
        timer.Start();
        if (v_PreIntegrate(
                step)) // Could be possible to put a preintegrate step in the
                       // ALEHelper class, put in the Unsteady Advection class
        {
            break;
        }

        ASSERTL0(m_timestep > 0, "m_timestep < 0");

        fields = m_intScheme->TimeIntegrate(stepCounter, m_timestep);
        timer.Stop();

        m_time += m_timestep;
        elapsed = timer.TimePerTest(1);
        intTime += elapsed;
        cpuTime += elapsed;

        // Write out status information.
        v_PrintStatusInformation(step, cpuTime);
        if (m_infosteps &&
            m_session->GetComm()->GetSpaceComm()->GetRank() == 0 &&
            !((step + 1) % m_infosteps))
        {
            cpuPrevious = cpuTime;
            cpuTime     = 0.0;
        }

        // Transform data into coefficient space
        for (i = 0; i < nvariables; ++i)
        {
            // copy fields into ExpList::m_phys and assign the new
            // array to fields
            m_indfields[m_intVariables[i]]->SetPhys(fields[i]);
            fields[i] = m_indfields[m_intVariables[i]]->UpdatePhys();
            if (v_RequireFwdTrans())
            {
                if (m_comm->IsParallelInTime())
                {
                    m_indfields[m_intVariables[i]]->FwdTrans(
                        m_indfields[m_intVariables[i]]->GetPhys(),
                        m_indfields[m_intVariables[i]]->UpdateCoeffs());
                }
                else
                {
                    m_indfields[m_intVariables[i]]->FwdTransLocalElmt(
                        m_indfields[m_intVariables[i]]->GetPhys(),
                        m_indfields[m_intVariables[i]]->UpdateCoeffs());
                }
            }
            m_indfields[m_intVariables[i]]->SetPhysState(false);
        }

        // Perform any solver-specific post-integration steps.
        if (v_PostIntegrate(step))
        {
            break;
        }

        // Test for abort conditions (nan, or abort file).
        if (m_abortSteps && !((step + 1) % m_abortSteps))
        {
            abortFlags[0] = 0;
            for (i = 0; i < nvariables; ++i)
            {
                if (Vmath::Nnan(
                        m_indfields[m_intVariables[i]]->GetPhys().size(),
                        m_indfields[m_intVariables[i]]->GetPhys(), 1) > 0)
                {
                    abortFlags[0] = 1;
                }
            }

            // Rank zero looks for abort file and deletes it
            // if it exists. Communicates the abort.
            if (m_session->GetComm()->GetSpaceComm()->GetRank() == 0)
            {
                if (fs::exists(abortFile))
                {
                    fs::remove(abortFile);
                    abortFlags[1] = 1;
                }
            }

            m_session->GetComm()->GetSpaceComm()->AllReduce(
                abortFlags, LibUtilities::ReduceMax);

            ASSERTL0(!abortFlags[0], "NaN found during time integration.");
        }

        // Write out checkpoint files.
        if ((m_checksteps && !((step + 1) % m_checksteps)) || doCheckTime)
        {

            Checkpoint_Output(m_nchk);
            m_nchk++;

            doCheckTime = false;
        }

        // Step advance.
        ++step;
        ++stepCounter;
    }

    // Print out summary statistics.
    v_PrintSummaryStatistics(intTime);

    // If homogeneous, transform back into physical space if necessary.

    for (i = 0; i < nvariables; ++i)
    {
        m_indfields[m_intVariables[i]]->SetPhysState(true);
    }
    if (this->particle_sys)
    {
        this->particle_sys->free();
    }
}

/**
 * @brief Override v_PostIntegrate to do particle output and compute
 * diagnostics, if enabled
 * @param step Time step number
 */
bool TokamakSystem::v_PostIntegrate(int step)
{
    this->solver_callback_handler.call_post_integrate(this);

    // Writes a step of the particle trajectory.

    return TimeEvoEqnSysBase<SU::UnsteadySystem,
                             ParticleSystem>::v_PostIntegrate(step);
}

/**
 * @brief Override v_PreIntegrate to do particle system integration
 *
 * @param step Time step number
 */
bool TokamakSystem::v_PreIntegrate(int step)
{
    this->solver_callback_handler.call_pre_integrate(this);

    if (this->transient_field)
    {
        this->mag_field->Update(m_time);
    }

    if (this->Te)
    {
        Vmath::Vdiv(n_pts,
                    m_indfields[m_indfields.size() - n_indep_fields]->GetPhys(),
                    1, m_fields[0]->GetPhys(), 1, Te->UpdatePhys(), 1);
        Te->FwdTransLocalElmt(this->Te->GetPhys(), Te->UpdateCoeffs());
    }

    if (this->particles_enabled)
    {
        // Integrate the particle system to the requested time.
        this->particle_sys->evaluate_fields(this->E, this->B, this->ne,
                                            this->Te, this->ve);
        if (particle_output_freq > 0 && (step % particle_output_freq) == 0)
        {
            this->particle_sys->write(step);
        }
        for (auto &fld : this->src_fields)
        {
            Vmath::Zero(fld->GetNpoints(), fld->UpdatePhys(), 1);
        }
        this->particle_sys->integrate(m_time + m_timestep, this->part_timestep);
    }

    return UnsteadySystem::v_PreIntegrate(step);
}

NESOSessionFunctionSharedPtr TokamakSystem::get_species_function(
    const std::string &s, std::string name, const MR::ExpListSharedPtr &field,
    bool cache)
{
    MR::ExpListSharedPtr vField = field;
    if (!field)
    {
        vField = m_fields[0];
    }

    if (cache)
    {
        if ((m_nesoSessionFunctions[s].find(name) ==
             m_nesoSessionFunctions[s].end()) ||
            (m_nesoSessionFunctions[s][name]->GetSession() !=
             this->neso_config) ||
            (m_nesoSessionFunctions[s][name]->GetExpansion() != vField))
        {
            m_nesoSessionFunctions[s][name] =
                MemoryManager<NESOSessionFunction>::AllocateSharedPtr(
                    s, this->neso_config, vField, name, cache);
        }

        return m_nesoSessionFunctions[s][name];
    }
    else
    {
        return std::make_shared<NESOSessionFunction>(s, this->neso_config,
                                                     vField, name, cache);
    }
}

/**
 * @brief After reading ICs, calculate phi and grad(phi)
 */
void TokamakSystem::v_SetInitialConditions(NekDouble init_time, bool dump_ICs,
                                           [[maybe_unused]] const int domain)
{
    if (m_session->DefinesFunction("InitialConditions"))
    {
        for (int i = 0; i < m_fields.size(); ++i)
        {
            Array<OneD, NekDouble> phys = m_fields[i]->UpdatePhys();
            Array<OneD, NekDouble> tmpphys(phys.size(), 0.0);
            for (size_t d = 0; d < m_domains.size(); ++d)
            {
                GetFunction("InitialConditions")
                    ->Evaluate(m_session->GetVariables()[i], tmpphys, m_time,
                               d);
                Vmath::Vcopy(m_dom_to_offset[d + 1] - m_dom_to_offset[d],
                             tmpphys.data() + m_dom_to_offset[d], 1,
                             phys.data() + m_dom_to_offset[d], 1);
            }

            if (m_fields[i]->GetWaveSpace())
            {
                m_fields[i]->HomogeneousFwdTrans(m_fields[i]->GetTotPoints(),
                                                 m_fields[i]->GetPhys(),
                                                 m_fields[i]->UpdatePhys());
            }
            m_fields[i]->FwdTransLocalElmt(m_fields[i]->GetPhys(),
                                           m_fields[i]->UpdateCoeffs());
        }

        // Enforce C0 Continuity of initial condiiton
        if ((m_projectionType == MultiRegions::eGalerkin) ||
            (m_projectionType == MultiRegions::eMixed_CG_Discontinuous))
        {
            for (int i = 0; i < m_fields.size(); ++i)
            {
                m_fields[i]->LocalToGlobal();
                m_fields[i]->GlobalToLocal();
                m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),
                                      m_fields[i]->UpdatePhys());
            }
        }

        if (m_session->GetComm()->GetRank() == 0)
        {
            std::cout << "Initial Conditions:" << std::endl;
            for (int i = 0; i < m_fields.size(); ++i)
            {
                std::string varName = m_session->GetVariable(i);
                std::cout << "  - Field " << varName << std::endl;
                for (size_t d = 0; d < m_domains.size(); ++d)
                {
                    std::cout << "\tDomain " << d << ": "
                              << GetFunction("InitialConditions")
                                     ->Describe(varName, d)
                              << std::endl;
                }
            }
        }
    }
    else
    {
        int nq = m_fields[0]->GetNpoints();
        for (int i = 0; i < m_fields.size(); i++)
        {
            Vmath::Zero(nq, m_fields[i]->UpdatePhys(), 1);
            m_fields[i]->SetPhysState(true);
            Vmath::Zero(m_fields[i]->GetNcoeffs(), m_fields[i]->UpdateCoeffs(),
                        1);
            if (m_session->GetComm()->GetRank() == 0)
            {
                std::cout << "Initial Conditions:" << std::endl;
                std::cout << "  - Field " << m_session->GetVariable(i)
                          << ": 0 (default)" << std::endl;
            }
        }
    }
    int s = 0;
    for (const auto &[k, v] : this->GetSpecies())
    {
        if (this->neso_config->defines_species_function(v.name,
                                                        "InitialConditions"))
        {
            auto fn = this->get_species_function(v.name, "InitialConditions");
            for (const auto &[f, fi] : v.fields)
            {
                Array<OneD, NekDouble> phys = m_indfields[fi]->UpdatePhys();
                Array<OneD, NekDouble> tmpphys(phys.size(), 0.0);
                for (size_t d = 0; d < m_domains.size(); ++d)
                {
                    fn->Evaluate(m_session->GetVariables()[f], tmpphys, m_time,
                                 d);
                    Vmath::Vcopy(m_dom_to_offset[d + 1] - m_dom_to_offset[d],
                                 tmpphys.data() + m_dom_to_offset[d], 1,
                                 phys.data() + m_dom_to_offset[d], 1);
                }
                if (m_indfields[fi]->GetWaveSpace())
                {
                    m_indfields[fi]->HomogeneousFwdTrans(
                        m_indfields[fi]->GetTotPoints(),
                        m_indfields[fi]->GetPhys(),
                        m_indfields[fi]->UpdatePhys());
                }
                m_indfields[fi]->FwdTransLocalElmt(
                    m_indfields[fi]->GetPhys(),
                    m_indfields[fi]->UpdateCoeffs());

                if ((m_projectionType == MultiRegions::eGalerkin) ||
                    (m_projectionType == MultiRegions::eMixed_CG_Discontinuous))
                {

                    m_indfields[fi]->LocalToGlobal();
                    m_indfields[fi]->GlobalToLocal();
                    m_indfields[fi]->BwdTrans(m_indfields[fi]->GetCoeffs(),
                                              m_indfields[fi]->UpdatePhys());
                }
                if (m_session->GetComm()->GetRank() == 0)
                {
                    std::string varName = m_session->GetVariable(f);
                    std::cout << "  - Species " << v.name << " Field "
                              << varName << std::endl;
                    for (size_t d = 0; d < m_domains.size(); ++d)
                    {
                        std::cout << "\tDomain " << d << ": "
                                  << this->get_species_function(
                                             v.name, "InitialConditions")
                                         ->Describe(varName, d)
                                  << std::endl;
                    }
                }
            }
        }
        else
        {
            int nq = m_indfields[0]->GetNpoints();
            for (const auto &[f, fi] : v.fields)
            {
                Vmath::Zero(nq, m_indfields[fi]->UpdatePhys(), 1);
                m_indfields[fi]->SetPhysState(true);
                Vmath::Zero(m_indfields[fi]->GetNcoeffs(),
                            m_indfields[fi]->UpdateCoeffs(), 1);
                if (m_session->GetComm()->GetRank() == 0)
                {
                    std::cout << "Initial Conditions:" << std::endl;
                    std::cout << "  - Field " << m_session->GetVariable(f)
                              << ": 0 (default)" << std::endl;
                }
            }
        }
        s++;
    }
    if (m_session->GetComm()->GetRank() == 0)
    {
        std::cout << "============================================="
                     "================"
                     "=========="
                  << std::endl
                  << std::endl;
    }
    if (dump_ICs && m_checksteps && m_nchk == 0 && !m_comm->IsParallelInTime())
    {
        Checkpoint_Output(m_nchk);
    }
    ++m_nchk;
}

void TokamakSystem::SetBoundaryConditions(NekDouble time)
{
    EquationSystem::SetBoundaryConditions(time);
    std::string varName;
    for (const auto &[k, v] : this->GetSpecies())
    {
        for (const auto &[f, fi] : v.fields)
        {
            varName = m_session->GetVariable(f);
            m_indfields[fi]->EvaluateBoundaryConditions(time, varName);
        }
    }
}

} // namespace NESO::Solvers::tokamak
