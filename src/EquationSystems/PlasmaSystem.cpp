#include <algorithm>
#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/BasicUtils/Vmath.hpp>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>

#include "PlasmaSystem.hpp"

namespace PENKNIFE
{

PlasmaSystem::PlasmaSystem(const LU::SessionReaderSharedPtr &session,
                           const SD::MeshGraphSharedPtr &graph)
    : SC::UnsteadySystem(session, graph),
      field_to_index(session->GetVariables())
{
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
                    m_expansionLists[0]->GetElmtToExpId().at(global.back()));

                phys_offset.push_back(
                    m_expansionLists[0]->GetPhys_Offset(local.back()));

                size +=
                    m_expansionLists[0]->GetExp(local.back())->GetTotPoints();
            }
        }
        m_dom_to_offset[d + 1] = size;
    }

    this->particles_enabled = false;
    this->neso_config       = std::make_shared<NESOReader>(session);

    this->neso_config->read_species();
    this->neso_config->read_vantage();

    // Session n/v/e variables describe templates. The evolved fields are
    // constructed once from the species declarations, followed by ee.
    m_variables.clear();
    for (const auto &[name, spec] : this->neso_config->get_species())
    {
        const unsigned int s = m_names.size();
        m_names.push_back(name);
        n_ids.push_back(InvalidId);
        v_ids.push_back(InvalidId);
        e_ids.push_back(InvalidId);
        m_neutral_ions.push_back(InvalidId);

        double mass, charge;
        this->neso_config->load_species_parameter(name, "Mass", mass, 1.0);
        this->neso_config->load_species_parameter(name, "Charge", charge, 0.0);
        m_masses.push_back(mass);
        m_charges.push_back(charge);

        for (const auto &field : std::get<2>(spec))
        {
            std::vector<unsigned int> *ids = nullptr;
            if (field == "n")
                ids = &n_ids;
            else if (field == "v")
                ids = &v_ids;
            else if (field == "e")
                ids = &e_ids;
            else
                ASSERTL0(false, "Unknown field '" + field + "' for " + name +
                                    "; expected n, v or e.");

            ASSERTL0((*ids)[s] == InvalidId,
                     "Duplicate field '" + field + "' for species " + name);
            (*ids)[s] = m_variables.size();
            m_variables.push_back(name + "_" + field);
        }

        if (std::get<3>(spec).has_value())
            m_neutrals.push_back(s);
        else
            m_ions.push_back(s);
    }
    n_species = m_names.size();

    // Resolve neutral associations after all names and ion IDs are known.
    // A neutral may sort before its associated ion in the configuration map.
    unsigned int s = 0;
    for (const auto &[name, spec] : this->neso_config->get_species())
    {
        if (const auto &ion_name = std::get<3>(spec); ion_name.has_value())
        {
            const auto ion =
                std::find(m_names.begin(), m_names.end(), *ion_name);
            ASSERTL0(ion != m_names.end(),
                     "Unknown ion '" + *ion_name + "' associated with " + name);
            const unsigned int ion_id = std::distance(m_names.begin(), ion);
            ASSERTL0(std::find(m_ions.begin(), m_ions.end(), ion_id) !=
                         m_ions.end(),
                     "Species '" + *ion_name + "' associated with " + name +
                         " is not an ion.");
            m_neutral_ions[s] = ion_id;
        }
        ++s;
    }

    ee_idx = m_variables.size();
    m_variables.push_back("ee");
    m_nVariables = m_variables.size();

    if (session->DefinesSolverInfo("PARTTYPE"))
    {
        NESOASSERT(this->neso_config->get_particle_species().size(),
                   "ParticleSystem specified in <INFO> but no kinetic species "
                   "found in <SPECIES>");

        std::string part_sys_name = session->GetSolverInfo("PARTTYPE");
        NESOASSERT(GetParticleSystemFactory().ModuleExists(part_sys_name),
                   "ParticleSystem '" + part_sys_name +
                       "' is not defined.\n"
                       "Ensure particle system name is correct and module is "
                       "compiled.\n");

        // The PartSysBase ptr returned from the factory is cast back to the
        // solver-specific PARTSYS type to allow the eqn_sys to use
        // solver-specific polymorphism
        this->particle_sys = GetParticleSystemFactory().CreateInstance(
            part_sys_name, neso_config, graph);
        this->particles_enabled = true;
        this->particle_sys->init_object();
    }
    else if (this->neso_config->get_particle_species().size())
    {
        NESOASSERT(
            false,
            "Kinetic species "
            "found in <SPECIES> but no ParticleSystem specified in <INFO>");
    }
}

void PlasmaSystem::v_GenerateSummary(SummaryList &s)
{
    std::stringstream ss;
    ss << R"(
╔═══════════════════════════════════════════════════════════════════════════╗
║                                                                           ║
║  ██████╗  ███████╗ ███╗   ██╗ ██╗  ██╗ ███╗   ██╗ ██╗ ███████╗ ███████╗   ║
║  ██╔══██╗ ██╔════╝ ████╗  ██║ ██║ ██╔╝ ████╗  ██║ ██║ ██╔════╝ ██╔════╝   ║
║  ██████╔╝ █████╗   ██╔██╗ ██║ █████╔╝  ██╔██╗ ██║ ██║ █████╗   █████╗     ║
║  ██╔═══╝  ██╔══╝   ██║╚██╗██║ ██╔═██╗  ██║╚██╗██║ ██║ ██╔══╝   ██╔══╝     ║
║  ██║      ███████╗ ██║ ╚████║ ██║  ██╗ ██║ ╚████║ ██║ ██║      ███████╗   ║
║  ╚═╝      ╚══════╝ ╚═╝  ╚═══╝ ╚═╝  ╚═╝ ╚═╝  ╚═══╝ ╚═╝ ╚═╝      ╚══════╝   ║
║ Plasma Edge Numerics using Kinetic Neutrals Integrated in Finite Elements ║
╚═══════════════════════════════════════════════════════════════════════════╝
)";
    AddSummaryItem(s, "", ss.str());
    SC::UnsteadySystem::v_GenerateSummary(s);

    AddSummaryItem(s, "Mesh length scale", this->mesh_length);
    AddSummaryItem(s, "Density normalisation", this->Nnorm);
    AddSummaryItem(s, "Temperature normalisation", this->Tnorm);
    AddSummaryItem(s, "Magnetic Field normalisation", this->Bnorm);
}

std::shared_ptr<ParticleSystem> PlasmaSystem::GetParticleSystem()
{
    return this->particle_sys;
}

/**
 * @brief Load all required session parameters into member variables.
 */
void PlasmaSystem::LoadParams()
{
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
}

void PlasmaSystem::v_InitialiseFields()
{
    EquationSystem::v_InitialiseFields();
    auto interleaveWidth =
        Operator<double>::GetDefaultInterleaveWidth(m_session);
    auto bAtr_phys = MultiRegions::GetBlockAttributes<double, FieldState::Phys>(
        m_expansionLists[0], interleaveWidth);

    this->B = LibUtilities::Field<double, FieldState::Phys>(
        "Magnetic Field", bAtr_phys, 3, m_npointsZ);
    this->b_unit = LibUtilities::Field<double, FieldState::Phys>(
        "Unit Magnetic Field", bAtr_phys, 3, m_npointsZ);

    this->mag_B = LibUtilities::Field<double, FieldState::Phys>(
        "Squared Magnetic Field Magnitude", bAtr_phys, 1, m_npointsZ);
    this->E = LibUtilities::Field<double, FieldState::Phys>(
        "Electric Field", bAtr_phys, 3, m_npointsZ);

    this->ne = LU::Field<double, FieldState::Phys>(
        "Electron Density", bAtr_phys, 1, m_npointsZ);
    m_math.zero(this->ne);
    this->j_par = LU::Field<double, FieldState::Phys>(
        "Parallel Current", bAtr_phys, 1, m_npointsZ);
    m_math.zero(this->j_par);
    this->Me = LU::Field<double, FieldState::Phys>(
        "Parallel Electron Momentum", bAtr_phys, 1, m_npointsZ);
    m_math.zero(this->Me);

    auto bAtr_phys_trace =
        MultiRegions::GetLocTraceBlockAttributes<double, FieldState::Phys>(
            m_expansionLists[0]);
    this->b_unit_trace = LibUtilities::Field<double, FieldState::Phys>(
        "Unit B on trace", bAtr_phys_trace, m_coordim, m_npointsZ);
}

void PlasmaSystem::v_InitialiseOperators()
{
    EquationSystem::v_InitialiseOperators();
    std::string execName =
        MultiRegions::Operator<double>::GetOpExecSpace(m_session);

    m_quasineutralityOp = QuasineutralityOp<double>::Create(
        m_expansionLists[0], m_variables, execName);
    std::vector<unsigned int> ionDensityIds;
    std::vector<unsigned int> ionMomentumIds;
    std::vector<double> ionCharges, ionMasses;
    ionDensityIds.reserve(m_ions.size());
    ionMomentumIds.reserve(m_ions.size());
    ionCharges.reserve(m_ions.size());
    ionMasses.reserve(m_ions.size());
    for (const auto s : m_ions)
    {
        ionDensityIds.push_back(GetDensityId(s));
        ionMomentumIds.push_back(v_ids[s]);
        ionCharges.push_back(m_charges[s]);
        ionMasses.push_back(m_masses[s]);
    }
    m_quasineutralityOp->SetSpecies(ionDensityIds, ionMomentumIds, ionCharges,
                                  ionMasses, me);

    m_BTraceExtractOp = PhysTraceExtractOp<double>::Create(
        m_expansionLists[0], this->b_unit.GetComponentNames());
    this->mag_field = MagneticFieldOp<double>::Create(
        m_expansionLists[0], this->B.GetComponentNames());
    this->mag_field->SetFields(this->B, this->b_unit, this->mag_B);
}

/**
 * @brief Construct the flux vector for the anisotropic diffusion problem.
 * @param fieldcoeffs field coefficients to be appended to
 * @param variables variable names to be appended to
 */
void PlasmaSystem::v_WriteFld(const std::string &outname)
{
    EquationSystem::v_WriteFld(outname);
}

/**
 * @brief Post-construction class initialisation.
 *
 * @param create_field if true, create a new field object and add it to
 * m_fields. Optional, defaults to true.
 */
void PlasmaSystem::v_InitObject(bool declareExpansionLists)
{
    SC::UnsteadySystem::v_InitObject(declareExpansionLists);
    // Load parameters
    LoadParams();

    ASSERTL0(m_session->DefinesSolverInfo("UPWINDTYPE"),
             "No UPWINDTYPE defined in session.");

    /// Create Field for solution m_fields and others
    InitialiseFields();

    // Create and initialise all operators
    InitialiseOperators();
    // Configure Time-integration
    InitialiseTimeOp();

    // this->b_unit = this->mag_field->b_unit;
    // this->mag_B  = this->mag_field->mag_B;
    this->mag_field->Apply(m_time);
    m_BTraceExtractOp->Apply(this->b_unit, this->b_unit_trace);

    // Forcing terms
    m_forcing = SC::Forcing::Load(m_session, m_expansionLists[0], m_variables);
}

void PlasmaSystem::v_DoInitialise(bool dump_initial_conditions)
{
    if (this->m_session->GetComm()->TreatAsRankZero() &&
        this->particles_enabled)
    {
        particle_sys->add_params_report();
    }
    UnsteadySystem::v_DoInitialise(dump_initial_conditions);
}

/**
 * @brief Override v_PostIntegrate to do particle output and compute
 * diagnostics, if enabled
 * @param step Time step number
 */
bool PlasmaSystem::v_PostIntegrate()
{
    return UnsteadySystem::v_PostIntegrate();
}

/**
 * @brief Override v_PreIntegrate to do particle system integration
 *
 * @param step Time step number
 */
bool PlasmaSystem::v_PreIntegrate()
{
    if (this->transient_field)
    {
        this->mag_field->Apply(m_time);
        m_BTraceExtractOp->Apply(this->b_unit, this->b_unit_trace);
    }

    if (this->particles_enabled)
    {
        m_math.zero(src_fields);
        this->particle_sys->integrate(m_time + m_timestep, m_timestep,
                                      m_timeOp->GetStep());
    }

    return UnsteadySystem::v_PreIntegrate();
}

/**
 * @brief Returns a session function particular to a species
 */
// NESOSessionFunctionSharedPtr PlasmaSystem::get_species_function(
//     const std::string &s, std::string name, const MR::ExpListSharedPtr
//     &field, bool cache)
// {
//     MR::ExpListSharedPtr vField = field;
//     if (!field)
//     {
//         vField = m_fields[0];
//     }

//     if (cache)
//     {
//         if ((m_nesoSessionFunctions[s].find(name) ==
//              m_nesoSessionFunctions[s].end()) ||
//             (m_nesoSessionFunctions[s][name]->GetSession() !=
//              this->neso_config) ||
//             (m_nesoSessionFunctions[s][name]->GetExpansion() != vField))
//         {
//             m_nesoSessionFunctions[s][name] =
//                 MemoryManager<NESOSessionFunction>::AllocateSharedPtr(
//                     s, this->neso_config, vField, name, cache);
//         }

//         return m_nesoSessionFunctions[s][name];
//     }
//     else
//     {
//         return std::make_shared<NESOSessionFunction>(s, this->neso_config,
//                                                      vField, name, cache);
//     }
// }

/**
 * @brief After reading ICs, calculate phi and grad(phi)
 */
// void PlasmaSystem::v_SetInitialConditions(NekDouble init_time)
// {
//     if (m_session->DefinesFunction("InitialConditions"))
//     {
//         for (int i = 0; i < m_fields.size(); ++i)
//         {
//             Array<OneD, NekDouble> phys = m_fields[i]->UpdatePhys();
//             Array<OneD, NekDouble> tmpphys(phys.size(), 0.0);
//             for (size_t d = 0; d < m_domains.size(); ++d)
//             {
//                 GetFunction("InitialConditions")
//                     ->Evaluate(m_session->GetVariables()[i], tmpphys, m_time,
//                                d);
//                 Vmath::Vcopy(m_dom_to_offset[d + 1] - m_dom_to_offset[d],
//                              tmpphys.data() + m_dom_to_offset[d], 1,
//                              phys.data() + m_dom_to_offset[d], 1);
//             }

//             if (m_fields[i]->GetWaveSpace())
//             {
//                 m_fields[i]->HomogeneousFwdTrans(m_fields[i]->GetTotPoints(),
//                                                  m_fields[i]->GetPhys(),
//                                                  m_fields[i]->UpdatePhys());
//             }
//             m_fields[i]->FwdTransLocalElmt(m_fields[i]->GetPhys(),
//                                            m_fields[i]->UpdateCoeffs());
//         }

//         // Enforce C0 Continuity of initial condiiton
//         if ((m_projectionType == MultiRegions::eGalerkin) ||
//             (m_projectionType == MultiRegions::eMixed_CG_Discontinuous))
//         {
//             for (int i = 0; i < m_fields.size(); ++i)
//             {
//                 m_fields[i]->AvgAssemble();
//                 m_fields[i]->GlobalToLocal();
//                 m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),
//                                       m_fields[i]->UpdatePhys());
//             }
//         }

//         if (m_session->GetComm()->GetRank() == 0)
//         {
//             std::cout << "Initial Conditions:" << std::endl;
//             for (int i = 0; i < m_fields.size(); ++i)
//             {
//                 std::string varName = m_session->GetVariable(i);
//                 std::cout << "  - Field " << varName << std::endl;
//                 for (size_t d = 0; d < m_domains.size(); ++d)
//                 {
//                     std::cout << "\tDomain " << d << ": "
//                               << GetFunction("InitialConditions")
//                                      ->Describe(varName, d)
//                               << std::endl;
//                 }
//             }
//         }
//     }
//     else
//     {
//         for (int i = 0; i < m_fields.size(); i++)
//         {
//             Vmath::Zero(this->n_pts, m_fields[i]->UpdatePhys(), 1);
//             m_fields[i]->SetPhysState(true);
//             Vmath::Zero(this->n_coeffs, m_fields[i]->UpdateCoeffs(), 1);
//             if (m_session->GetComm()->GetRank() == 0)
//             {
//                 std::cout << "Initial Conditions:" << std::endl;
//                 std::cout << "  - Field " << m_session->GetVariable(i)
//                           << ": 0 (default)" << std::endl;
//             }
//         }
//     }
//     int s = 0;
//     for (const auto &[k, v] : this->GetSpecies())
//     {
//         if (this->neso_config->defines_species_function(v.name,
//                                                         "InitialConditions"))
//         {
//             auto fn = this->get_species_function(v.name,
//             "InitialConditions"); for (const auto &[f, fi] : v.fields)
//             {
//                 Array<OneD, NekDouble> phys = m_indfields[fi]->UpdatePhys();
//                 Array<OneD, NekDouble> tmpphys(phys.size(), 0.0);
//                 for (size_t d = 0; d < m_domains.size(); ++d)
//                 {
//                     fn->Evaluate(m_session->GetVariables()[f], tmpphys,
//                     m_time,
//                                  d);
//                     Vmath::Vcopy(m_dom_to_offset[d + 1] - m_dom_to_offset[d],
//                                  tmpphys.data() + m_dom_to_offset[d], 1,
//                                  phys.data() + m_dom_to_offset[d], 1);
//                 }
//                 if (m_indfields[fi]->GetWaveSpace())
//                 {
//                     m_indfields[fi]->HomogeneousFwdTrans(
//                         m_indfields[fi]->GetTotPoints(),
//                         m_indfields[fi]->GetPhys(),
//                         m_indfields[fi]->UpdatePhys());
//                 }
//                 m_indfields[fi]->FwdTransLocalElmt(
//                     m_indfields[fi]->GetPhys(),
//                     m_indfields[fi]->UpdateCoeffs());

//                 if ((m_projectionType == MultiRegions::eGalerkin) ||
//                     (m_projectionType ==
//                     MultiRegions::eMixed_CG_Discontinuous))
//                 {

//                     m_indfields[fi]->AvgAssemble();
//                     m_indfields[fi]->GlobalToLocal();
//                     m_indfields[fi]->BwdTrans(m_indfields[fi]->GetCoeffs(),
//                                               m_indfields[fi]->UpdatePhys());
//                 }
//                 if (m_session->GetComm()->GetRank() == 0)
//                 {
//                     std::string varName = m_session->GetVariable(f);
//                     std::cout << "  - Species " << v.name << " Field "
//                               << varName << std::endl;
//                     for (size_t d = 0; d < m_domains.size(); ++d)
//                     {
//                         std::cout << "\tDomain " << d << ": "
//                                   << this->get_species_function(
//                                              v.name, "InitialConditions")
//                                          ->Describe(varName, d)
//                                   << std::endl;
//                     }
//                 }
//             }
//         }
//         else
//         {
//             for (const auto &[f, fi] : v.fields)
//             {
//                 Vmath::Zero(this->n_pts, m_indfields[fi]->UpdatePhys(), 1);
//                 m_indfields[fi]->SetPhysState(true);
//                 Vmath::Zero(this->n_coeffs, m_indfields[fi]->UpdateCoeffs(),
//                 1); if (m_session->GetComm()->GetRank() == 0)
//                 {
//                     std::cout << "Initial Conditions:" << std::endl;
//                     std::cout << "  - Field " << m_session->GetVariable(f)
//                               << ": 0 (default)" << std::endl;
//                 }
//             }
//         }
//         s++;
//     }
//     Vmath::Zero(this->n_pts, m_fields[0]->UpdatePhys(), 1);
//     for (const auto &[s, v] : GetIons())
//     {
//         int ni_idx = v.fields.at(field_to_index.at("n"));

//         Vmath::Svtvp(this->n_pts, v.charge, m_indfields[ni_idx]->GetPhys(),
//         1,
//                      m_fields[0]->UpdatePhys(), 1, m_fields[0]->UpdatePhys(),
//                      1);
//     }
//     m_fields[0]->FwdTrans(m_fields[0]->GetPhys(),
//     m_fields[0]->UpdateCoeffs());

//     if (m_session->GetComm()->GetRank() == 0)
//     {
//         std::cout << "============================================="
//                      "================"
//                      "=========="
//                   << std::endl
//                   << std::endl;
//     }
//     if (dump_ICs && m_checksteps && m_nchk == 0 &&
//     !m_comm->IsParallelInTime())
//     {
//         Checkpoint_Output(m_nchk);
//     }
//     ++m_nchk;
// }

} // namespace PENKNIFE
