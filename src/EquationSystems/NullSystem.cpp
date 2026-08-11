#include "NullSystem.hpp"

namespace PENKNIFE
{

/// Name of class
static std::string class_name;
std::string NullSystem::class_name =
    SU::GetEquationSystemFactory().RegisterCreatorFunction(
        "NullSystem", NullSystem::create,
        "Only advances particles with no fluid solve");
/**
 * @brief Creates an instance of this class.
 */
static SU::EquationSystemSharedPtr create(
    const LU::SessionReaderSharedPtr &session,
    const SD::MeshGraphSharedPtr &graph)
{
    SU::EquationSystemSharedPtr p =
        MemoryManager<NullSystem>::AllocateSharedPtr(session, graph);
    p->InitObject();
    return p;
}

NullSystem::NullSystem(const LU::SessionReaderSharedPtr &session,
                       const SD::MeshGraphSharedPtr &graph)
    : PlasmaSystem(session, graph)
{
}

/**
 * @brief load necessary parameters
 */
void NullSystem::load_params()
{
    PlasmaSystem::load_params();
    auto variables = m_session->GetVariables();
    if (std::find(variables.begin(), variables.end(), "e") != variables.end())
    {
        this->n_indep_fields = 1;
    }
    else
    {
        this->n_indep_fields = 0;
    }
}

/**
 * @brief Initialise the class.
 */
void NullSystem::v_InitObject(bool DeclareFields)
{
    PlasmaSystem::v_InitObject(DeclareFields);

    if (this->n_indep_fields)
        this->ee_idx = m_indfields.size() - this->n_indep_fields;

    int npoints = m_indfields[0]->GetNpoints();

    m_session->MatchSolverInfo("SpectralVanishingViscosity", "True",
                               m_useSpecVanVisc, false);
    m_session->LoadParameter("epsilon", m_epsilon, 1.0);
    if (m_useSpecVanVisc)
    {
        m_session->LoadParameter("SVVCutoffRatio", m_sVVCutoffRatio, 0.75);
        m_session->LoadParameter("SVVDiffCoeff", m_sVVDiffCoeff, 0.1);
    }
    // Setup diffusion object
    m_ode.DefineOdeRhs(&NullSystem::DoOdeRhs, this);

    if (this->particles_enabled)
    {
        this->ne = std::dynamic_pointer_cast<MR::DisContField>(m_fields[0]);
        this->particle_sys->setup_evaluate_fields(this->E, this->B, this->ne,
                                                  this->Te, this->ve);
        std::vector<Sym<REAL>> src_syms;
        std::vector<int> src_components;
        std::vector<Sym<REAL>> out_syms;

        int cnt = 0;
        for (const auto &[s, v] : this->GetIons())
        {

            this->src_fields.emplace_back(
                MemoryManager<MR::DisContField>::AllocateSharedPtr(
                    *std::dynamic_pointer_cast<MR::DisContField>(m_fields[0])));
            src_syms.push_back(Sym<REAL>(v.name + "_SOURCE_DENSITY"));
            src_components.push_back(0);
            ni_src_idx.push_back(cnt++);
            out_syms.push_back(Sym<REAL>(v.name + "_SOURCE_DENSITY"));

            for (int d = 0; d < this->m_spacedim; ++d)
            {
                this->src_fields.emplace_back(
                    MemoryManager<MR::DisContField>::AllocateSharedPtr(
                        *std::dynamic_pointer_cast<MR::DisContField>(
                            m_fields[0])));
                src_syms.push_back(Sym<REAL>(v.name + "_SOURCE_MOMENTUM"));
                src_components.push_back(d);
            }
            vi_src_idx.push_back(cnt);
            cnt += m_spacedim;
            out_syms.push_back(Sym<REAL>(v.name + "_SOURCE_MOMENTUM"));

            this->src_fields.emplace_back(
                MemoryManager<MR::DisContField>::AllocateSharedPtr(
                    *std::dynamic_pointer_cast<MR::DisContField>(m_fields[0])));

            src_syms.push_back(Sym<REAL>(v.name + "_SOURCE_ENERGY"));
            src_components.push_back(0);
            ei_src_idx.push_back(cnt++);
            out_syms.push_back(Sym<REAL>(v.name + "_SOURCE_ENERGY"));
        }

        this->src_fields.emplace_back(
            MemoryManager<MR::DisContField>::AllocateSharedPtr(
                *std::dynamic_pointer_cast<MR::DisContField>(m_fields[0])));
        src_syms.push_back(Sym<REAL>("ELECTRON_SOURCE_DENSITY"));
        src_components.push_back(0);
        out_syms.push_back(Sym<REAL>("ELECTRON_SOURCE_ENERGY"));

        this->particle_sys->finish_setup(this->src_fields, src_syms,
                                         src_components);

        std::vector<int> diag_components = {0};
        std::vector<Sym<REAL>> diag_syms = {Sym<REAL>("WEIGHT")};

        for (const auto &[s, v] : this->particle_sys->get_species())
        {
            this->diag_fields[v.id].emplace_back(
                MemoryManager<MR::DisContField>::AllocateSharedPtr(
                    *std::dynamic_pointer_cast<MR::DisContField>(m_fields[0])));
        }
        this->particle_sys->diag_setup(this->diag_fields, diag_syms,
                                       diag_components);
        this->particle_sys->output_setup(out_syms);
    }
}

/**
 * @brief Populate rhs array ( @p out_arr ) with diffused quantitities
 * @param in_arr physical values of all fields
 * @param[out] out_arr output array (RHSs of time integration equations)
 * @param time simulation time
 */
void NullSystem::DoOdeRhs(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    // RHS should be set to zero.
    for (int i = 0; i < outarray.size(); ++i)
    {
        Vmath::Zero(outarray[i].size(), &outarray[i][0], 1);
    }

    if (this->particles_enabled)
    {
        // DoParticles(inarray, outarray);
    }

    // Add forcing terms
    for (auto &x : m_forcing)
    {
        x->Apply(m_fields, inarray, outarray, time);
    }
}

/**
 * @brief Compute the advection terms for the right-hand side
 */
void NullSystem::DoParticles(const Array<OneD, Array<OneD, NekDouble>> &inarray,
                             Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    int cnt = 0;
    for (const auto &[s, v] : this->GetIons())
    {
        int ni_idx = v.fields.at(field_to_index["n"]);
        //  Add contribution to ion density
        Vmath::Vadd(this->n_pts, outarray[ni_idx], 1,
                    this->src_fields[ni_src_idx[cnt]]->GetPhys(), 1,
                    outarray[ni_idx], 1);
        if (v.fields.find(field_to_index.get_idx("v")) != v.fields.end())
        {
            int vi_idx = v.fields.at(field_to_index["v"]);

            for (int d = 0; d < m_spacedim; ++d)
            {
                Vmath::Vvtvp(this->n_pts, this->b_unit[d], 1,
                             this->src_fields[vi_src_idx[cnt] + d]->GetPhys(),
                             1, outarray[vi_idx], 1, outarray[vi_idx], 1);
            }
        }

        if (v.fields.find(field_to_index.get_idx("e")) != v.fields.end())
        {
            int ei_idx = v.fields.at(field_to_index["e"]);

            // Add contribution to ion energy
            Vmath::Vadd(this->n_pts, outarray[ei_idx], 1,
                        this->src_fields[ei_src_idx[cnt]]->GetPhys(), 1,
                        outarray[ei_idx], 1);

            // Add number density source contribution to ion energy
            Array<OneD, NekDouble> dynamic_energy(this->n_pts);
            m_varConv->GetIonDynamicEnergy(s, v.mass, inarray, dynamic_energy);
            Vmath::Vvtvp(this->n_pts, dynamic_energy, 1,
                         this->src_fields[ni_src_idx[cnt]]->GetPhys(), 1,
                         outarray[ei_idx], 1, outarray[ei_idx], 1);
        }
        cnt++;
    }
    // Add contribution to electron energy
    if (this->n_indep_fields)
        Vmath::Vadd(this->n_pts, outarray[ee_idx], 1,
                    this->src_fields.back()->GetPhys(), 1, outarray[ee_idx], 1);
}

/**
 * @brief Post-integration step
 */
bool NullSystem::v_PostIntegrate(int step)
{
    Vmath::Zero(this->n_pts, m_fields[0]->UpdatePhys(), 1);
    for (const auto &[s, v] : GetIons())
    {
        int ni_idx = v.fields.at(field_to_index.at("n"));

        Vmath::Svtvp(this->n_pts, v.charge, m_indfields[ni_idx]->GetPhys(), 1,
                     m_fields[0]->UpdatePhys(), 1, m_fields[0]->UpdatePhys(),
                     1);
    }
    m_fields[0]->FwdTrans(m_fields[0]->GetPhys(), m_fields[0]->UpdateCoeffs());

    if (this->particles_enabled)
    {
        this->particle_sys->diag_project();
    }
    // Writes a step of the particle trajectory.

    return PlasmaSystem::v_PostIntegrate(step);
}

/**
 * @brief Write extra fields.
 * @param fieldcoeffs field coefficients to be appended to
 * @param variables variable names to be appended to
 */
void NullSystem::v_ExtraFldOutput(
    std::vector<Array<OneD, NekDouble>> &fieldcoeffs,
    std::vector<std::string> &variables)
{
    PlasmaSystem::v_ExtraFldOutput(fieldcoeffs, variables);
    const int nPhys   = m_fields[0]->GetNpoints();
    const int nCoeffs = m_fields[0]->GetNcoeffs();

    if (this->particles_enabled)
    {
        int cnt = 0;

        for (auto &[k, v] : this->GetIons())
        {
            variables.push_back(v.name + "_SOURCE_DENSITY");
            Array<OneD, NekDouble> SrcFwd(nCoeffs);
            m_fields[0]->FwdTransLocalElmt(
                this->src_fields[ni_src_idx[cnt]]->GetPhys(), SrcFwd);
            fieldcoeffs.push_back(SrcFwd);

            for (int d = 0; d < this->m_spacedim; ++d)
            {
                variables.emplace_back(v.name + "_SOURCE_MOMENTUM" +
                                       std::to_string(d));
                Array<OneD, NekDouble> SrcFwd(nCoeffs);
                m_fields[0]->FwdTransLocalElmt(
                    this->src_fields[vi_src_idx[cnt] + d]->GetPhys(), SrcFwd);
                fieldcoeffs.emplace_back(SrcFwd);
            }

            variables.emplace_back(v.name + "_SOURCE_ENERGY");
            Array<OneD, NekDouble> SrcFwd2(nCoeffs);
            m_fields[0]->FwdTransLocalElmt(
                this->src_fields[ei_src_idx[cnt]]->GetPhys(), SrcFwd2);
            fieldcoeffs.emplace_back(SrcFwd2);
        }
        variables.push_back("ELECTRON_SOURCE_ENERGY");
        Array<OneD, NekDouble> ESrcFwd(nCoeffs);
        m_fields[0]->FwdTransLocalElmt(this->src_fields.back()->GetPhys(),
                                       ESrcFwd);
        fieldcoeffs.push_back(ESrcFwd);

        for (auto &[k, v] : this->particle_sys->get_species())
        {
            variables.emplace_back(k + "_DENSITY");
            Array<OneD, NekDouble> DiagFwd(nCoeffs);
            m_fields[0]->FwdTransLocalElmt(
                this->diag_fields[v.id][0]->GetPhys(), DiagFwd);
            fieldcoeffs.push_back(DiagFwd);
        }
    }
}
} // namespace PENKNIFE