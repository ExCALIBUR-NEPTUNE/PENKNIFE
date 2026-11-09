#include <FieldUtils/Interpolator.h>
#include <SolverUtils/Core/SessionFunction.h>
#include <algorithm>
#include <boost/core/ignore_unused.hpp>
#include <nektar_interface/solver_base/neso_session_function.hpp>
#include <set>

#include <LibUtilities/BasicUtils/Vmath.hpp>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>

#include "PlasmaSystem.hpp"

namespace PENKNIFE
{

PlasmaSystem::PlasmaSystem(const LU::SessionReaderSharedPtr &session,
                           const SD::MeshGraphSharedPtr &graph)
    : SC::UnsteadySystem(session, graph)
{
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
            m_fieldTemplates.push_back(field);
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
    m_fieldTemplates.push_back("e");
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

void PlasmaSystem::v_GenerateSummary(SC::SummaryList &s)
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
    SC::AddSummaryItem(s, "", ss.str());
    SC::UnsteadySystem::v_GenerateSummary(s);

    SC::AddSummaryItem(s, "Mesh length scale", this->mesh_length);
    SC::AddSummaryItem(s, "Density normalisation", this->Nnorm);
    SC::AddSummaryItem(s, "Temperature normalisation", this->Tnorm);
    SC::AddSummaryItem(s, "Magnetic Field normalisation", this->Bnorm);
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
    m_session->LoadParameter("T_bg", m_backgroundTemperature, 0.0);

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

    this->ne = LU::Field<double, FieldState::Phys>("Electron Density",
                                                   bAtr_phys, 1, m_npointsZ);
    m_math.zero(this->ne);
    this->j_par = LU::Field<double, FieldState::Phys>("Parallel Current",
                                                      bAtr_phys, 1, m_npointsZ);
    m_math.zero(this->j_par);
    this->Me = LU::Field<double, FieldState::Phys>("Parallel Electron Momentum",
                                                   bAtr_phys, 1, m_npointsZ);
    m_math.zero(this->Me);
    m_math.zero(this->E);

    auto bAtr_phys_trace =
        MultiRegions::GetLocTraceBlockAttributes<double, FieldState::Phys>(
            m_expansionLists[0]);
    this->b_unit_trace = LibUtilities::Field<double, FieldState::Phys>(
        "Unit B on trace", bAtr_phys_trace, 3, m_npointsZ);
}

void PlasmaSystem::v_InitialiseOperators()
{
    const auto variables = m_variables;
    m_variables          = m_fieldTemplates;
    EquationSystem::v_InitialiseOperators();
    m_variables = variables;
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
    // Initialise Nektar's session-sized metadata before mapping the evolved
    // species components onto the existing n/v/e expansion templates. This
    // also keeps its restart scan within the original session variable list.
    const auto variables = m_variables;
    const auto count     = m_nVariables;
    const auto templates = m_session->GetVariables();
    m_variables          = templates;
    m_nVariables         = templates.size();
    SC::UnsteadySystem::v_InitObject(declareExpansionLists);
    m_variables  = variables;
    m_nVariables = count;
    if (declareExpansionLists)
    {
        const auto expansions = m_expansionLists;
        const auto singular   = m_checkIfSystemSingular;
        m_expansionLists.clear();
        m_checkIfSystemSingular.clear();
        for (const auto &field : m_fieldTemplates)
        {
            const auto it =
                std::find(templates.begin(), templates.end(), field);
            ASSERTL0(it != templates.end(),
                     "Missing session field template: " + field);
            const auto index = std::distance(templates.begin(), it);
            m_expansionLists.push_back(expansions[index]);
            m_checkIfSystemSingular.push_back(singular[index]);
        }
    }
    // Load parameters
    LoadParams();

    /// Create Field for solution m_fields and others
    InitialiseFields();

    // Create and initialise all operators
    InitialiseOperators();
    // Configure Time-integration
    InitialiseTimeOp();

    UpdateMagneticField(m_time);
    SetupParticles();

    // Forcing terms
    m_forcing =
        SC::Forcing::Load(m_session, m_expansionLists[0], m_fieldTemplates);
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
    if (particles_enabled)
        particle_sys->diag_project();
    return UnsteadySystem::v_PostIntegrate();
}

/**
 * @brief Override v_PreIntegrate to do particle system integration
 *
 * @param step Time step number
 */
bool PlasmaSystem::v_PreIntegrate()
{
    if (transient_field)
        UpdateMagneticField(m_time);
    j_par.SetInterleaveWidth(m_fields);
    m_quasineutralityOp->Apply(m_fields, j_par, ne, Me);
    UpdateElectricField(m_fields, m_time);

    if (this->particles_enabled)
    {
        UpdateParticleFields();
        m_math.zero(m_particleSources);
        m_math.zero(m_particleSourceCoeffs);
        this->particle_sys->integrate(m_time + m_timestep, m_timestep,
                                      m_timeOp->GetStep());
    }

    return UnsteadySystem::v_PreIntegrate();
}

void PlasmaSystem::UpdateMagneticField(double time)
{
    mag_field->Apply(time);
    m_BTraceExtractOp->Apply(b_unit, b_unit_trace);
}

void PlasmaSystem::v_InitialiseTimeOp()
{
    const auto &scheme = m_session->GetTimeIntScheme();
    if (scheme.method == "RungeKutta" || scheme.method == "AdamsBashforth" ||
        scheme.method == "ExplicitSDC" || scheme.method == "ExplicitGEM")
        UnsteadySystem::v_InitialiseTimeOp();
    else
    {
        m_timeOp = SC::TimeOp<double>::Create(m_expansionLists[0], m_variables,
                                              "RungeKutta", 4, "");
        m_timeOp->DefineProjection(&SC::EquationSystem::DoProjection, this);
        WARNINGL0(false, "PENKNIFE currently uses explicit RK4 in place of the "
                         "requested implicit or IMEX scheme.");
    }
}

void PlasmaSystem::SetupParticles()
{
    if (!particles_enabled)
        return;
    ASSERTL0(m_npointsZ == 1,
             "NESO Field coupling requires non-homogeneous fields.");
    auto attr =
        MR::GetBlockAttributes<double, FieldState::Phys>(m_expansionLists[0]);
    std::vector<Sym<REAL>> evalSyms, sourceSyms, outputSyms;
    std::vector<int> evalComps, sourceComps;
    std::vector<std::string> evalNames, sourceNames, projectionNames;
    auto evaluate = [&](const std::string &name, int count)
    {
        for (int d = 0; d < count; ++d)
        {
            evalSyms.emplace_back(name);
            evalComps.push_back(d);
            evalNames.push_back(name + std::to_string(d));
        }
    };
    evaluate("ELECTRON_DENSITY", 1);
    evaluate("ELECTRON_TEMPERATURE", 1);
    evaluate("ELECTRON_FLOW_SPEED", 3);
    evaluate("ELECTRIC_FIELD", 3);
    evaluate("MAGNETIC_FIELD", 3);
    for (const auto &name : m_names)
    {
        evaluate(name + "_DENSITY", 1);
        evaluate(name + "_TEMPERATURE", 1);
        evaluate(name + "_FLOW_SPEED", 3);
    }
    auto source = [&](const std::string &name, unsigned count,
                      const std::string &variable)
    {
        outputSyms.emplace_back(name);
        const unsigned first = sourceNames.size();
        for (unsigned d = 0; d < count; ++d)
        {
            sourceSyms.emplace_back(name);
            sourceComps.push_back(d);
            sourceNames.push_back(name + std::to_string(d));
            projectionNames.push_back(variable);
        }
        return first;
    };
    m_sourceDensity.assign(n_species, InvalidId);
    m_sourceMomentum.assign(n_species, InvalidId);
    m_sourceEnergy.assign(n_species, InvalidId);
    for (unsigned s = 0; s < n_species; ++s)
    {
        if (HasDensity(s))
            m_sourceDensity[s] = source(m_names[s] + "_SOURCE_DENSITY", 1,
                                        m_fieldTemplates[n_ids[s]]);
        if (HasMomentum(s))
            m_sourceMomentum[s] = source(m_names[s] + "_SOURCE_MOMENTUM", 3,
                                         m_fieldTemplates[v_ids[s]]);
        if (HasEnergy(s))
            m_sourceEnergy[s] = source(m_names[s] + "_SOURCE_ENERGY", 1,
                                       m_fieldTemplates[e_ids[s]]);
    }
    if (ee_idx != InvalidId)
        m_sourceElectron =
            source("ELECTRON_SOURCE_ENERGY", 1, m_fieldTemplates[ee_idx]);
    m_particleFields =
        PhysField("Particle evaluation fields", attr, evalNames, 1);
    m_particleSources = PhysField("Particle sources", attr, projectionNames, 1);

    m_particleSourceCoeffs = CoeffField(
        "Particle source coefficients",
        MR::GetBlockAttributes<double, FieldState::Coeff>(m_expansionLists[0]),
        projectionNames, 1);
    m_particleSourceNames = sourceNames;
    m_math.zero(m_particleSources);
    m_math.zero(m_particleSourceCoeffs);
    particle_sys->setup_evaluate_fields(m_expansionLists[0], m_particleFields,
                                        evalSyms, evalComps);
    particle_sys->finish_setup(m_expansionLists[0], m_particleSources,
                               m_particleSourceCoeffs, sourceSyms, sourceComps);
    particle_sys->diag_setup();
    particle_sys->output_setup(outputSyms);
}

void PlasmaSystem::UpdateParticleFields()
{
    m_particleFields.SetInterleaveWidth(m_fields);
    for (auto *field : {&B, &E, &b_unit, &ne, &Me})
        field->SetInterleaveWidth(m_fields);
    m_math.zero(m_particleFields);
    for (unsigned blk = 0; blk < m_fields.GetBlocks().size(); ++blk)
    {
        auto &block         = m_fields.GetBlocks()[blk];
        const size_t stride = block.CompSize();
        auto u              = block.GetPtr<NektarSpaces::HostSpace, ReadOnly>();
        auto n =
            ne.GetBlocks()[blk].GetPtr<NektarSpaces::HostSpace, ReadOnly>();
        auto momentum =
            Me.GetBlocks()[blk].GetPtr<NektarSpaces::HostSpace, ReadOnly>();
        auto b =
            b_unit.GetBlocks()[blk].GetPtr<NektarSpaces::HostSpace, ReadOnly>();
        auto magnetic =
            B.GetBlocks()[blk].GetPtr<NektarSpaces::HostSpace, ReadOnly>();
        auto electric =
            E.GetBlocks()[blk].GetPtr<NektarSpaces::HostSpace, ReadOnly>();
        auto values = m_particleFields.GetBlocks()[blk]
                          .GetPtr<NektarSpaces::HostSpace, ReadWrite>();
        for (size_t p = 0; p < stride; ++p)
        {
            values[p]          = n[p];
            values[stride + p] = n[p] > 0 && ee_idx != InvalidId
                                     ? (2. / 3.) * u[ee_idx * stride + p] / n[p]
                                     : m_backgroundTemperature;
            for (unsigned d = 0; d < 3; ++d)
            {
                values[(2 + d) * stride + p] =
                    n[p] > 0 ? momentum[p] * b[d * stride + p] / (me * n[p])
                             : 0;
                values[(5 + d) * stride + p] = electric[d * stride + p];
                values[(8 + d) * stride + p] = magnetic[d * stride + p];
            }
            for (unsigned s = 0; s < n_species; ++s)
            {
                const size_t c  = 11 + 5 * s;
                const double ns = HasDensity(s) ? u[n_ids[s] * stride + p] : 0;
                values[c * stride + p] = ns;
                values[(c + 1) * stride + p] =
                    ns > 0 && HasEnergy(s)
                        ? (2. / 3.) * u[e_ids[s] * stride + p] / ns
                        : m_backgroundTemperature;
                for (unsigned d = 0; d < 3; ++d)
                    values[(c + 2 + d) * stride + p] =
                        ns > 0 && HasMomentum(s)
                            ? u[v_ids[s] * stride + p] * b[d * stride + p] /
                                  (m_masses[s] * ns)
                            : 0;
            }
        }
    }
}

void PlasmaSystem::AddSources(PhysField &in, PhysField &out, double time,
                              double factor)
{
    if (particles_enabled)
    {
        m_particleSources.SetInterleaveWidth(out);
        b_unit.SetInterleaveWidth(out);
        for (unsigned blk = 0; blk < out.GetBlocks().size(); ++blk)
        {
            const auto stride = out.GetBlocks()[blk].CompSize();
            auto source       = m_particleSources.GetBlocks()[blk]
                              .GetPtr<NektarSpaces::HostSpace, ReadOnly>();
            auto b = b_unit.GetBlocks()[blk]
                         .GetPtr<NektarSpaces::HostSpace, ReadOnly>();
            auto rhs = out.GetBlocks()[blk]
                           .GetPtr<NektarSpaces::HostSpace, ReadWrite>();
            for (size_t p = 0; p < stride; ++p)
            {
                for (unsigned s = 0; s < n_species; ++s)
                {
                    if (HasDensity(s))
                        rhs[n_ids[s] * stride + p] +=
                            factor * source[m_sourceDensity[s] * stride + p] /
                            m_timestep;
                    if (HasEnergy(s))
                        rhs[e_ids[s] * stride + p] +=
                            factor * source[m_sourceEnergy[s] * stride + p] /
                            m_timestep;
                    if (HasMomentum(s))
                        for (unsigned d = 0; d < 3; ++d)
                            rhs[v_ids[s] * stride + p] +=
                                factor * b[d * stride + p] *
                                source[(m_sourceMomentum[s] + d) * stride + p] /
                                m_timestep;
                }
                if (m_sourceElectron != InvalidId)
                    rhs[ee_idx * stride + p] +=
                        factor * source[m_sourceElectron * stride + p] /
                        m_timestep;
            }
        }
    }
    for (auto &forcing : m_forcing)
    {
        forcing->SetAppend(true);
        forcing->Apply(in, out, time, factor);
    }
}

SC::CoefficientFieldViews PlasmaSystem::v_GetOutputFieldViews()
{
    auto views = EquationSystem::v_GetOutputFieldViews();
    if (particles_enabled)
    {
        views.push_back({m_particleSourceCoeffs, m_expansionLists[0], m_coordim,
                         SC::MakeFieldComponents(m_particleSourceNames)});
        particle_sys->append_diagnostics(views);
    }
    return views;
}

void PlasmaSystem::v_SetInitialConditions(double time)
{
    m_time = time;
    m_math.zero(m_fields);
    m_math.zero(m_fields_coeff);
    const auto expansion = m_expansionLists[0];
    auto collections     = MR::GetCollections(expansion);
    auto attrs = MR::GetBlockAttributes<double, FieldState::Phys>(expansion);
    PhysField values("Initial component", attrs, 1, 1);
    m_math.zero(values);
    auto expression = MR::ExpressionOp<double>::Create(expansion, {"initial"});
    for (const auto &[domain, composites] : m_graph->GetDomain())
    {
        std::set<int> geometries;
        for (const auto &[id, composite] : composites)
            for (const auto &geom : composite->m_geomVec)
                geometries.insert(geom->GetGlobalID());
        for (unsigned f = 0; f < m_variables.size(); ++f)
        {
            std::string species, field = m_fieldTemplates[f];
            for (unsigned s = 0; s < n_species; ++s)
                if (n_ids[s] == f || v_ids[s] == f || e_ids[s] == f)
                {
                    species = m_names[s];
                    break;
                }
            const bool specific =
                !species.empty() && neso_config->defines_species_function(
                                        species, "InitialConditions");
            const bool generic =
                m_session->DefinesFunction("InitialConditions", m_variables[f],
                                           domain) ||
                m_session->DefinesFunction("InitialConditions", field, domain);
            if (!specific && !generic)
                continue;
            const std::string variable =
                m_session->DefinesFunction("InitialConditions", m_variables[f],
                                           domain)
                    ? m_variables[f]
                    : field;
            auto type = specific
                            ? neso_config->get_species_function_type(
                                  species, "InitialConditions", field, domain)
                            : m_session->GetFunctionType("InitialConditions",
                                                         variable, domain);
            if (type == LU::eFunctionTypeExpression)
            {
                auto equation =
                    specific ? neso_config->get_species_function(
                                   species, "InitialConditions", field, domain)
                             : m_session->GetFunction("InitialConditions",
                                                      variable, domain);
                expression->SetExpressions({equation});
                expression->SetTime(time);
                expression->Apply(values, values);
            }
            else
            {
                // File import is setup only; copy by element into the native
                // Field blocks, rather than assuming expansion/collection
                // order.
                Array<OneD, double> imported(expansion->GetTotPoints(), 0.);
                if (specific)
                {
                    NESO::NESOSessionFunction function(
                        species, neso_config, expansion, "InitialConditions");
                    function.Evaluate(field, imported, time, domain);
                }
                else
                {
                    SolverUtils::SessionFunction function(m_session, expansion,
                                                          "InitialConditions");
                    function.Evaluate(variable, imported, time, domain);
                }
                values.Initialize<NektarSpaces::HostSpace>(0.);
                for (unsigned blk = 0; blk < collections.size(); ++blk)
                {
                    auto &block = values.GetBlocks()[blk];
                    auto dst =
                        block.GetPtr<NektarSpaces::HostSpace, ReadWrite>();
                    const size_t width   = block.GetInterleaveWidth(),
                                 np      = block.GetNumData();
                    const auto &elements = collections[blk].GetExpVector();
                    for (unsigned e = 0; e < elements.size(); ++e)
                    {
                        const auto local = expansion->GetElmtToExpId().at(
                            elements[e]->GetGeom()->GetGlobalID());
                        const size_t offset = expansion->GetPhys_Offset(local);
                        for (size_t p = 0; p < np; ++p)
                            dst[(e / width) * np * width + p * width +
                                e % width] = imported[offset + p];
                    }
                }
            }
            values.SetInterleaveWidth(m_fields);
            for (unsigned blk = 0; blk < collections.size(); ++blk)
            {
                auto &block = m_fields.GetBlocks()[blk];
                auto dst = block.GetPtr<NektarSpaces::HostSpace, ReadWrite>();
                auto src = values.GetBlocks()[blk]
                               .GetPtr<NektarSpaces::HostSpace, ReadOnly>();
                const size_t width   = block.GetInterleaveWidth(),
                             np      = block.GetNumData();
                const auto &elements = collections[blk].GetExpVector();
                for (unsigned e = 0; e < elements.size(); ++e)
                    if (geometries.count(elements[e]->GetGeom()->GetGlobalID()))
                        for (size_t p = 0; p < np; ++p)
                        {
                            const size_t point = (e / width) * np * width +
                                                 p * width + e % width;
                            dst[f * block.CompSize() + point] = src[point];
                        }
            }
        }
    }
    // Use the redesigned projection for CG continuity and stage boundary data.
    DoProjection(m_fields, m_fields, time);
    m_quasineutralityOp->Apply(m_fields, j_par, ne, Me);
    UpdateElectricField(m_fields, time);
    if (particles_enabled)
        UpdateParticleFields();
}

} // namespace PENKNIFE
