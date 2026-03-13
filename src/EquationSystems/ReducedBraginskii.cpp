#include "ReducedBraginskii.hpp"
#include "../RiemannSolvers/PlasmaSolver.hpp"
#include <SolverUtils/Advection/AdvectionNonConservative.h>
#include <SolverUtils/Advection/AdvectionWeakDG.h>

namespace PENKNIFE
{

/// Name of class
static std::string class_name;
std::string ReducedBraginskii::class_name =
    SU::GetEquationSystemFactory().RegisterCreatorFunction(
        "ReducedBraginskii", ReducedBraginskii::create,
        "Solves reduced Braginskii equations");
/**
 * @brief Creates an instance of this class.
 */
static SU::EquationSystemSharedPtr create(
    const LU::SessionReaderSharedPtr &session,
    const SD::MeshGraphSharedPtr &graph)
{
    SU::EquationSystemSharedPtr p =
        MemoryManager<ReducedBraginskii>::AllocateSharedPtr(session, graph);
    p->InitObject();
    return p;
}

ReducedBraginskii::ReducedBraginskii(const LU::SessionReaderSharedPtr &session,
                                     const SD::MeshGraphSharedPtr &graph)
    : PlasmaSystem(session, graph)
{
    this->n_indep_fields = 1; // p_e
}

void ReducedBraginskii::v_InitObject(bool DeclareFields)
{
    PlasmaSystem::v_InitObject(DeclareFields);

    InitBraginskii();

    std::string diffName;
    m_session->LoadSolverInfo("DiffusionType", diffName, "LDG");
    m_diffusion =
        SolverUtils::GetDiffusionFactory().CreateInstance(diffName, diffName);
    m_diffusion->SetFluxVector(&ReducedBraginskii::GetFluxVectorDiff, this);

    // workaround for bug in DiffusionLDG
    m_difffields = Array<OneD, MR::ExpListSharedPtr>(m_indfields.size());
    for (int f = 0; f < m_difffields.size(); ++f)
    {
        m_difffields[f] = m_indfields[f];
    }

    m_diffusion->InitObject(m_session, m_difffields);
    this->ee_idx = m_indfields.size() - this->n_indep_fields;

    // Create storage for velocities
    this->m_kpar   = Array<OneD, NekDouble>(this->n_pts, 0.0);
    this->m_kcross = Array<OneD, NekDouble>(this->n_pts, 0.0);
    this->m_kperp  = Array<OneD, NekDouble>(this->n_pts, 0.0);

    InitAdvection();

    m_ode.DefineOdeRhs(&ReducedBraginskii::DoOdeRhs, this);

    if (this->particles_enabled)
    {
        std::vector<Sym<REAL>> src_syms;
        std::vector<int> src_components;

        int cnt = 0;
        for (const auto &[s, v] : this->GetIons())
        {
            this->src_fields.emplace_back(
                MemoryManager<MR::DisContField>::AllocateSharedPtr(
                    *std::dynamic_pointer_cast<MR::DisContField>(m_fields[0])));
            src_syms.push_back(Sym<REAL>(v.name + "_SOURCE_DENSITY"));
            src_components.push_back(0);
            ni_src_idx.push_back(cnt++);

            if (v.fields.find(field_to_index["v"]) != v.fields.end())
            {
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
            }
            if (v.fields.find(field_to_index["e"]) != v.fields.end())
            {
                this->src_fields.emplace_back(
                    MemoryManager<MR::DisContField>::AllocateSharedPtr(
                        *std::dynamic_pointer_cast<MR::DisContField>(
                            m_fields[0])));

                src_syms.push_back(Sym<REAL>(v.name + "_SOURCE_ENERGY"));
                src_components.push_back(0);
                ei_src_idx.push_back(cnt++);
            }
        }
        this->src_fields.emplace_back(
            MemoryManager<MR::DisContField>::AllocateSharedPtr(
                *std::dynamic_pointer_cast<MR::DisContField>(m_fields[0])));
        src_syms.push_back(Sym<REAL>("ELECTRON_SOURCE_ENERGY"));
        src_components.push_back(0);

        this->particle_sys->finish_setup(this->src_fields, src_syms,
                                         src_components);
    }
}

void ReducedBraginskii::InitAdvection()
{
    for (const auto &[s, v] : this->GetSpecies())
    {
        if (v.fields.find(field_to_index["v"]) != v.fields.end())
        {
            int ni_idx = v.fields.at(field_to_index["n"]);
            int vi_idx = v.fields.at(field_to_index["v"]);
            this->advected_fields.push_back(ni_idx);
            this->advected_fields.push_back(vi_idx);

            if (v.fields.find(field_to_index["e"]) != v.fields.end())
            {
                int ei_idx = v.fields.at(field_to_index["e"]);
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
void ReducedBraginskii::DoOdeRhs(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    // Get field indices
    int nTracePts = GetTraceTotPoints();
    for (int f = 0; f < outarray.size(); ++f)
    {
        Vmath::Zero(this->n_pts, outarray[f], 1);
    }

    int nvariables = inarray.size();

    m_varConv->GetElectronDensity(inarray, m_fields[0]->UpdatePhys());

    // Store forwards/backwards space along trace space
    Array<OneD, Array<OneD, NekDouble>> Fwd(nvariables);
    Array<OneD, Array<OneD, NekDouble>> Bwd(nvariables);

    for (int i = 0; i < nvariables; ++i)
    {
        Fwd[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
        Bwd[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
        m_indfields[i]->GetFwdBwdTracePhys(inarray[i], Fwd[i], Bwd[i]);
    }

    // Calculate E
    ComputeE();
    // Calculate ExB, parallel and diamagnetic velocities
    CalcVelocities(inarray);

    // Perform advection
    DoAdvection(inarray, outarray, time, Fwd, Bwd);

    m_bndConds->Update(inarray, time);

    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Neg(this->n_pts, outarray[i], 1);
    }

    // Perform Diffusion
    DoDiffusion(inarray, outarray, Fwd, Bwd);

    if (this->particles_enabled)
    {
        DoParticles(inarray, outarray);
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
void ReducedBraginskii::DoAdvection(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
    const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    const Array<OneD, Array<OneD, NekDouble>> &pBwd)
{
    int nvariables = this->advected_fields.size();
    int nTracePts  = GetTraceTotPoints();

    Array<OneD, Array<OneD, NekDouble>> advVel(m_spacedim);

    Array<OneD, Array<OneD, NekDouble>> outarrayAdv(nvariables);
    for (int i = 0; i < nvariables; ++i)
    {
        outarrayAdv[i] = Array<OneD, NekDouble>(this->n_pts, 0.0);
    }

    Array<OneD, Array<OneD, NekDouble>> inarrayAdv(nvariables);
    Array<OneD, Array<OneD, NekDouble>> inFwd(nvariables);
    Array<OneD, Array<OneD, NekDouble>> inBwd(nvariables);

    for (int i = 0; i < nvariables; ++i)
    {
        inarrayAdv[i] = inarray[advected_fields[i]];
        inFwd[i]      = pFwd[advected_fields[i]];
        inBwd[i]      = pBwd[advected_fields[i]];
    }

    m_advection->Advect(nvariables, m_advfields, advVel, inarrayAdv,
                        outarrayAdv, time, inFwd, inBwd);

    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Vadd(this->n_pts, outarrayAdv[i], 1,
                    outarray[this->advected_fields[i]], 1,
                    outarray[this->advected_fields[i]], 1);
    }
}

/**
 * @brief Compute the advection terms for the right-hand side
 */
void ReducedBraginskii::DoParticles(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    // Add contribution to electron energy
    Vmath::Vadd(this->n_pts, outarray[ee_idx], 1,
                this->src_fields[0]->GetPhys(), 1, outarray[ee_idx], 1);

    for (const auto &[s, v] : this->GetIons())
    {
        int ni_idx = v.fields.at(field_to_index["n"]);
        //  Add contribution to ion density
        Vmath::Vadd(this->n_pts, outarray[ni_idx], 1,
                    this->src_fields[ni_src_idx[s]]->GetPhys(), 1,
                    outarray[ni_idx], 1);

        if (v.fields.find(field_to_index["v"]) != v.fields.end())
        {
            int vi_idx = v.fields.at(field_to_index["v"]);

            for (int d = 0; d < m_spacedim; ++d)
            {
                Vmath::Vvtvp(this->n_pts, this->b_unit[d], 1,
                             this->src_fields[vi_src_idx[s] + d]->GetPhys(), 1,
                             outarray[vi_idx], 1, outarray[vi_idx], 1);
            }
        }

        if (v.fields.find(field_to_index["e"]) != v.fields.end())
        {
            int ei_idx = v.fields.at(field_to_index["e"]);

            // Add contribution to ion energy
            Vmath::Vadd(this->n_pts, outarray[ei_idx], 1,
                        this->src_fields[ei_src_idx[s]]->GetPhys(), 1,
                        outarray[ei_idx], 1);

            // Add number density source contribution to ion energy
            Array<OneD, NekDouble> dynamic_energy(this->n_pts);
            m_varConv->GetIonDynamicEnergy(s, v.mass, inarray, dynamic_energy);
            Vmath::Vvtvp(this->n_pts, dynamic_energy, 1,
                         this->src_fields[ni_src_idx[s]]->GetPhys(), 1,
                         outarray[ei_idx], 1, outarray[ei_idx], 1);
        }
    }
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

void ReducedBraginskii::CalcVelocities(
    const Array<OneD, Array<OneD, NekDouble>> &inarray)
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
    for (const auto &[s, v] : GetIons())
    {
        int ni_idx = v.fields.at(field_to_index["n"]);
        int vi_idx = v.fields.at(field_to_index["v"]);
        int ei_idx = v.fields.at(field_to_index["e"]);

        for (int p = 0; p < this->n_pts; ++p)
        {
            j_i[p] += v.charge * inarray[vi_idx][p] / v.mass;
            double v_i_par = inarray[vi_idx][p] / (v.mass * inarray[ni_idx][p]);
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

    for (const auto &[s, v] : GetNeutrals())
    {
        int nn_idx = v.fields.at(field_to_index["n"]);
        int vn_idx = v.fields.at(field_to_index["v"]);
        int en_idx = v.fields.at(field_to_index["e"]);

        for (int p = 0; p < this->n_pts; ++p)
        {
            double v_i_par = inarray[vn_idx][p] / (v.mass * inarray[nn_idx][p]);
            for (int d = 0; d < m_spacedim; ++d)
            {
                this->adv_vel[nn_idx][d][p] = v_i_par * this->b_unit[d][p];
                this->adv_vel[vn_idx][d][p] = v_i_par * this->b_unit[d][p];
                this->adv_vel[en_idx][d][p] = v_i_par * this->b_unit[d][p];
            }
        }
    }
}

/**
 *  @brief Compute components of advection velocities normal to trace elements
 * (faces, in 3D).
 *
 * @param[in,out] trace_vel_norm Trace normal velocities for each field
 * @param         adv_vel_trace        Advection velocities for each field
 */
Array<OneD, Array<OneD, NekDouble>> &ReducedBraginskii::GetAdvVelNorm()
{
    // Number of trace (interface) points
    int num_trace_pts = GetTraceNpoints();
    // Auxiliary variable to compute normal velocities

    // Compute advection vel dot trace normals and store
    for (int j = 0; j < this->adv_vel_trace.size(); ++j)
    {
        Array<OneD, Array<OneD, NekDouble>> normals(m_spacedim);

        for (int d = 0; d < m_spacedim; ++d)
        {
            normals[d] = Array<OneD, NekDouble>(num_trace_pts);
        }
        m_indfields[j]->GetTrace()->GetNormals(normals);
        // Ensure output array is zeroed
        Vmath::Zero(num_trace_pts, this->trace_vel_norm[j], 1);
        for (int d = 0; d < this->adv_vel_trace[j].size(); ++d)
        {
            m_indfields[j]->ExtractTracePhys(this->adv_vel[j][d],
                                             this->adv_vel_trace[j][d]);
            for (int p = 0; p < num_trace_pts; ++p)
            {
                this->trace_vel_norm[j][p] +=
                    normals[d][p] * this->adv_vel_trace[j][d][p];
            }
        }
    }
    return this->trace_vel_norm;
}

/**
 *  @brief Construct flux array.
 *
 * @param  field_vals Physical values for each advection field
 * @param[out] flux       Flux array
 */
void ReducedBraginskii::GetFluxVector(
    const Array<OneD, Array<OneD, NekDouble>> &field_vals,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fluxes)
{
    for (int d = 0; d < m_spacedim; ++d)
    {
        for (int p = 0; p < this->n_pts; ++p)
        {
            fluxes[ee_idx][d][p] =
                this->adv_vel[ee_idx][d][p] * field_vals[ee_idx][p];
        }
    }

    for (const auto &[s, v] : this->GetSpecies())
    {
        int ni_idx = v.fields.at(field_to_index["n"]);
        int vi_idx = v.fields.at(field_to_index["v"]);
        int ei_idx = v.fields.at(field_to_index["e"]);

        for (int d = 0; d < m_spacedim; ++d)
        {
            for (int p = 0; p < this->n_pts; ++p)
            {
                fluxes[ni_idx][d][p] =
                    this->adv_vel[ni_idx][d][p] * field_vals[ni_idx][p];
                fluxes[vi_idx][d][p] =
                    this->adv_vel[vi_idx][d][p] * field_vals[vi_idx][p];
                fluxes[ei_idx][d][p] =
                    this->adv_vel[ei_idx][d][p] * field_vals[ei_idx][p];
            }
        }
    }
}

void ReducedBraginskii::DoDiffusion(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray,
    const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    const Array<OneD, Array<OneD, NekDouble>> &pBwd)
{
    int nvariables = inarray.size();
    int nTracePts  = GetTraceTotPoints();

    // this should be preallocated
    Array<OneD, Array<OneD, NekDouble>> outarrayDiff(nvariables);
    for (int i = 0; i < nvariables; ++i)
    {
        outarrayDiff[i] = Array<OneD, NekDouble>(this->n_pts, 0.0);
    }

    Array<OneD, Array<OneD, NekDouble>> inarrayDiff(nvariables);
    Array<OneD, Array<OneD, NekDouble>> inFwd(nvariables);
    Array<OneD, Array<OneD, NekDouble>> inBwd(nvariables);

    for (int i = 0; i < nvariables; ++i)
    {
        inarrayDiff[i] = Array<OneD, NekDouble>(this->n_pts, 0.0);
        inFwd[i]       = Array<OneD, NekDouble>(nTracePts, 0.0);
        inBwd[i]       = Array<OneD, NekDouble>(nTracePts, 0.0);
    }

    // Extract temperature
    m_varConv->GetElectronTemperature(inarray, inarrayDiff[ee_idx]);

    for (const auto &[s, v] : this->GetIons())
    {
        int ni_idx = v.fields.at(field_to_index["n"]);
        int ei_idx = v.fields.at(field_to_index["e"]);

        inarrayDiff[ni_idx] = inarray[ni_idx];
        m_varConv->GetIonTemperature(s, v.mass, inarray, inarrayDiff[ei_idx]);
    }

    // Repeat calculation for trace space
    if (pFwd == NullNekDoubleArrayOfArray || pBwd == NullNekDoubleArrayOfArray)
    {
        inFwd = NullNekDoubleArrayOfArray;
        inBwd = NullNekDoubleArrayOfArray;
    }
    else
    {
        m_varConv->GetElectronTemperature(pFwd, inFwd[ee_idx]);
        m_varConv->GetElectronTemperature(pBwd, inBwd[ee_idx]);
        for (const auto &[s, v] : this->GetIons())
        {
            int ni_idx    = v.fields.at(field_to_index["n"]);
            int ei_idx    = v.fields.at(field_to_index["e"]);
            inFwd[ni_idx] = pFwd[ni_idx];
            inBwd[ni_idx] = pBwd[ni_idx];

            m_varConv->GetIonTemperature(s, v.mass, pFwd, inFwd[ei_idx]);
            m_varConv->GetIonTemperature(s, v.mass, pBwd, inBwd[ei_idx]);
        }
    }

    m_diffusion->Diffuse(nvariables, m_difffields, inarrayDiff, outarrayDiff,
                         inFwd, inBwd);

    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Vadd(this->n_pts, outarrayDiff[i], 1, outarray[i], 1,
                    outarray[i], 1);
    }
}

void ReducedBraginskii::CalcDiffTensor()
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Array<OneD, NekDouble> d(n_pts, 0.0);
            for (int k = 0; k < this->n_pts; k++)
            {
                d[k] = (m_kpar[k] - m_kperp[k]) * b_unit[i][k] * b_unit[j][k];
                if (i == j)
                {
                    d[k] += m_kperp[k];
                }
            }
            m_D[vc[i][j]] = d;
        }
    }
}

inline double CoulombLog_ii(double Nnorm, double ni1, double ni2, double Ti1,
                            double Ti2, double A1, double A2, double Z1,
                            double Z2)
{
    return 29.91 - log(sqrt(Nnorm)) -
           log((Z1 * Z2 * (A1 + A2)) / (A1 * Ti2 + A2 * Ti1) *
               sqrt(ni1 * Z1 * Z1 / Ti1 + ni2 * Z2 * Z2 / Ti2));
}

inline double CoulombLog_ee(double Nnorm, double Tnorm, double ne, double Te)
{
    double logTe = log(Tnorm * Te);
    return 30.4 - 0.5 * log(ne) - 0.5 * log(Nnorm) + (5. / 4) * logTe -
           sqrt(1e-5 + (logTe - 2) * (logTe - 2) / 16.);
}

inline double CoulombLog_ei(double Nnorm, double Tnorm, double ni, double ne,
                            double Ti, double Te, double Ai, double Zi)
{
    if ((Te * Tnorm < 0.1) || (ni * Nnorm < 1e10) || (ne * Nnorm < 1e10))
        return 10;
    else if (Te < Ti * constants::m_e_m_p / Ai)
        return 23 - 0.5 * log(ni) + 1.5 * log(Ti) - log(Zi * Zi * Ai) -
               0.5 * log(Nnorm) + 1.5 * log(Tnorm);
    else if (Te * Tnorm < exp(2) * Zi * Zi)
        return 30.0 - 0.5 * log(ne) - log(Zi) + 1.5 * log(Te) -
               0.5 * log(Nnorm) + 1.5 * log(Tnorm);
    else
        return 31.0 - 0.5 * log(ne) + log(Te) - 0.5 * log(Nnorm) + log(Tnorm);
}

void ReducedBraginskii::InitBraginskii()
{
    this->nu_ee = Array<OneD, NekDouble>(this->n_pts);
    this->nu_e  = Array<OneD, NekDouble>(this->n_pts);
    for (const auto &[s, v] : this->GetIons())
    {
        this->nu_i[s]  = Array<OneD, NekDouble>(this->n_pts);
        this->nu_ei[s] = Array<OneD, NekDouble>(this->n_pts);
        for (const auto &[s2, v2] : this->GetIons())
        {
            this->nu_ii[std::minmax(s, s2)] =
                Array<OneD, NekDouble>(this->n_pts);
        }
    }
}

void ReducedBraginskii::CalcCollisionFrequencies(
    const Array<OneD, Array<OneD, NekDouble>> &in_arr)
{
    auto ne = this->m_fields[0]->GetPhys();
    for (int p = 0; p < this->n_pts; ++p)
    {
        const double v1sq =
            2 * in_arr[ee_idx][p] * constants::e / constants::m_e_si;
        double coulomb_log =
            CoulombLog_ee(Nnorm, Tnorm, ne[p], in_arr[ee_idx][p]);

        // Electon collision frequency
        double nu = pow(constants::e, 4) * ne[p] * coulomb_log * 2 /
                    (3 * pow(M_PI * 2 * v1sq, 1.5) *
                     pow(constants::epsilon_0_si * constants::m_e_si, 2));
        this->nu_ee[p] = nu / omega_c;
        this->nu_e[p]  = this->nu_ee[p];
    }
    for (const auto &[s, v] : this->GetIons())
    {
        double Z   = this->m_ions[s].charge;
        double A   = this->m_ions[s].mass;
        int ni_idx = v.fields.at(field_to_index["n"]);
        int ei_idx = v.fields.at(field_to_index["e"]);
        for (int p = 0; p < this->n_pts; ++p)
        {
            const double vesq =
                2 * in_arr[ei_idx][p] * constants::e / constants::m_e_si;
            const double visq =
                2 * in_arr[ei_idx][p] * constants::e / (A * constants::m_p_si);
            double coulomb_log =
                CoulombLog_ei(Nnorm, Tnorm, in_arr[ni_idx][p], ne[p],
                              in_arr[ei_idx][p], in_arr[ee_idx][p], A, Z);
            // Collision frequency
            double nu = Z * Z * pow(constants::e, 4) * in_arr[ni_idx][p] *
                        coulomb_log * (1. + constants::m_e_m_p) /
                        (3 * pow(M_PI * (vesq + visq), 1.5) *
                         pow(constants::epsilon_0_si * constants::m_e_si, 2));
            nu /= omega_c;
            this->nu_ei[s][p] = nu;
            this->nu_e[p] += nu;
            this->nu_i[s][p] =
                constants::m_e_m_p * ne[p] * nu / in_arr[ni_idx][p];
        }
    }
    for (const auto &[s, v] : this->GetIons())
    {
        double Z   = this->m_ions[s].charge;
        double A   = this->m_ions[s].mass;
        int ni_idx = v.fields.at(field_to_index["n"]);
        int ei_idx = v.fields.at(field_to_index["e"]);
        for (const auto &[s2, v2] : this->GetIons())
        {
            if (s2 > s)
                continue;
            double Z2   = this->m_ions[s2].charge;
            double A2   = this->m_ions[s2].mass;
            int ni_idx2 = v2.fields.at(field_to_index["n"]);
            int ei_idx2 = v2.fields.at(field_to_index["e"]);
            for (int p = 0; p < this->n_pts; ++p)
            {
                double coulomb_log = CoulombLog_ii(
                    Nnorm, in_arr[ni_idx][p], in_arr[ni_idx2][p],
                    in_arr[ei_idx][p], in_arr[ei_idx2][p], A, A2, Z, Z2);

                const double v1sq = 2 * in_arr[ei_idx][p] * constants::e /
                                    (A * constants::m_p_si);
                const double v2sq = 2 * in_arr[ei_idx2][p] * constants::e /
                                    (A2 * constants::m_p_si);
                double nu =
                    Z * Z * Z2 * Z2 * pow(constants::e, 4) *
                    in_arr[ni_idx2][p] * coulomb_log * (1. + A / A2) /
                    (3 * pow(M_PI * (v1sq + v2sq), 1.5) *
                     pow(constants::epsilon_0_si * A * constants::m_p_si, 2));
                nu /= omega_c;
                this->nu_ii[std::minmax(s, s2)][p] = nu;
                this->nu_i[s][p] += nu;
                this->nu_i[s2][p] +=
                    (A / A2) * nu * in_arr[ni_idx][p] / in_arr[ni_idx][p];
            }
        }
    }
}

void ReducedBraginskii::CalcK(const Array<OneD, Array<OneD, NekDouble>> &in_arr,
                              int f)
{
    double Z   = this->m_ions[f].charge;
    double A   = this->m_ions[f].mass;
    int ni_idx = this->m_ions[f].fields[field_to_index["n"]];
    auto ne    = this->m_fields[0]->GetPhys();

    for (int p = 0; p < this->n_pts; ++p)
    {
        m_kpar[p] = this->k_ci * this->k_par * pow(in_arr[ee_idx][p], 2.5) /
                    (Z * Z * in_arr[ni_idx][p]);
        m_kperp[p] = this->k_perp * Z * Z * std::sqrt(A) * in_arr[ni_idx][p] /
                     (sqrt(in_arr[ee_idx][p]) * this->mag_B[p]);
        m_kcross[p] =
            this->k_cross * ne[p] * in_arr[ee_idx][p] / (sqrt(this->mag_B[p]));
    }
    if (m_session->DefinesParameter("k_perp"))
    {
        double k = m_session->GetParameter("k_perp");
        Vmath::Fill(this->n_pts, k, m_kperp, 1);
    }
    if (m_session->DefinesParameter("k_par"))
    {
        double k = m_session->GetParameter("k_par");
        Vmath::Fill(this->n_pts, k, m_kpar, 1);
    }
}

// Ion thermal conductivity
void ReducedBraginskii::CalcKappa(
    const Array<OneD, Array<OneD, NekDouble>> &in_arr, int f)
{
    double Z   = this->m_ions[f].charge;
    double A   = this->m_ions[f].mass;
    int ei_idx = this->m_ions[f].fields[field_to_index["e"]];
    int ni_idx = this->m_ions[f].fields[field_to_index["n"]];

    for (int p = 0; p < this->n_pts; ++p)
    {
        this->m_kpar[p] =
            k_ci * in_arr[ee_idx][p] * in_arr[ni_idx][p] / nu_i[f][p];

        // this->m_kperp[p] = this->kappa_i_perp * sqrt(A) * tmp[p] *
        //                    in_arr[ni_idx][p] /
        //                    (this->mag_B[p] * sqrt(in_arr[ei_idx][p]));
        // this->m_kcross[p] = this->kappa_i_cross * in_arr[ni_idx][p] *
        //                     in_arr[ei_idx][p] / (Z * sqrt(this->mag_B[p]));
    }
    if (m_session->DefinesParameter("kappa_perp"))
    {
        double k = m_session->GetParameter("kappa_perp");
        Vmath::Fill(this->n_pts, k, m_kperp, 1);
    }
    if (m_session->DefinesParameter("kappa_par"))
    {
        double k = m_session->GetParameter("kappa_par");
        Vmath::Fill(this->n_pts, k, m_kpar, 1);
    }
}

// Electron thermal conductivity
void ReducedBraginskii::CalcKappa(
    const Array<OneD, Array<OneD, NekDouble>> &in_arr)
{
    auto ne = this->m_fields[0]->GetPhys();
    for (int p = 0; p < this->n_pts; ++p)
    {
        this->m_kpar[p] = k_ce * in_arr[ee_idx][p] * ne[p] / nu_e[p];

        this->m_kperp[p] = this->kappa_e_perp * ne[p] * ne[p] /
                           (this->mag_B[p] * sqrt(in_arr[ee_idx][p]));
        this->m_kcross[p] = this->kappa_e_cross * ne[p] * in_arr[ee_idx][p] /
                            (sqrt(this->mag_B[p]));
    }
}

/**
 * @brief Construct the flux vector for the anisotropic diffusion problem.
 */
void ReducedBraginskii::GetFluxVectorDiff(
    const Array<OneD, Array<OneD, NekDouble>> &in_arr,
    const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &qfield,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fluxes)
{
    unsigned int nDim = qfield.size();
    unsigned int nFld = qfield[0].size();
    CalcCollisionFrequencies(in_arr);

    for (const auto &[s, v] : this->GetIons())
    {
        int ni_idx = v.fields.at(field_to_index["n"]);
        int ei_idx = v.fields.at(field_to_index["e"]);

        // CalcK(in_arr, s);
        // CalcDiffTensor();

        // for (unsigned int j = 0; j < nDim; ++j)
        // {
        //     for (unsigned int k = 0; k < nDim; ++k)
        //     {
        //         const Array<OneD, NekDouble> &D = m_D[vc[j][k]].GetValue();
        //         for (int p = 0; p < this->n_pts; ++p)
        //         {
        //             fluxes[j][ni_idx][p] += D[p] * qfield[k][ni_idx][p];
        //         }
        //     }
        // }

        // if (nDim == 3)
        // {
        //     for (int p = 0; p < this->n_pts; ++p)
        //     {
        //         fluxes[0][ei_idx][p] = b_unit[1][p] * qfield[2][ei_idx][p] -
        //                                b_unit[2][p] * qfield[1][ei_idx][p];
        //         fluxes[1][ei_idx][p] = b_unit[2][p] * qfield[0][ei_idx][p] -
        //                                b_unit[0][p] * qfield[2][ei_idx][p];
        //         fluxes[2][ei_idx][p] = b_unit[0][p] * qfield[1][ei_idx][p] -
        //                                b_unit[1][p] * qfield[0][ei_idx][p];
        //     }
        // }
        // else
        // {
        //     for (int p = 0; p < this->n_pts; ++p)
        //     {
        //         fluxes[0][ei_idx][p] = -b_unit[2][p] * qfield[1][ei_idx][p];
        //         fluxes[1][ei_idx][p] = b_unit[2][p] * qfield[0][ei_idx][p];
        //     }
        // }

        CalcKappa(in_arr, s);
        CalcDiffTensor();
        for (unsigned int j = 0; j < nDim; ++j)
        {
            // Calc diffusion of n with D tensor and n field
            for (unsigned int k = 0; k < nDim; ++k)
            {
                const Array<OneD, NekDouble> &D = m_D[vc[j][k]].GetValue();
                for (int p = 0; p < this->n_pts; ++p)
                {
                    fluxes[j][ei_idx][p] += D[p] * qfield[k][ei_idx][p];
                }
            }
        }
    }
    CalcKappa(in_arr);
    CalcDiffTensor();
    for (unsigned int j = 0; j < nDim; ++j)
    {
        // Calc diffusion of n with D tensor and n field
        for (unsigned int k = 0; k < nDim; ++k)
        {
            const Array<OneD, NekDouble> &D = m_D[vc[j][k]].GetValue();
            for (int p = 0; p < this->n_pts; ++p)
            {
                fluxes[j][ee_idx][p] += D[p] * qfield[k][ee_idx][p];
            }
        }
    }
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
        double kCX = 3.2e-15 * std::sqrt(ei[p] / 0.026);

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
        double kCX = 3.2e-15 * std::sqrt(2.0 / 0.026);
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
        double kCX = 3.2e-15 * std::sqrt(ei[p] / 0.026);

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
        double kCX = 3.2e-15 * std::sqrt(ei[p] / 0.026);

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

    for (const auto &[s, v] : this->GetNeutrals())
    {
        int nn_idx = v.fields.at(field_to_index["n"]);
        int ni_idx = this->m_ions[v.ion].fields.at(field_to_index["n"]);

        if (v.fields.find(field_to_index["e"]) != v.fields.end())
        {
            int en_idx = v.fields.at(field_to_index["e"]);
            int ei_idx = this->m_ions[v.ion].fields.at(field_to_index["e"]);

            if (v.fields.find(field_to_index["v"]) != v.fields.end())
            {
                // Neutral n,v,p; Ion n,v,p
                int vn_idx = v.fields.at(field_to_index["v"]);
                int vi_idx = this->m_ions[v.ion].fields.at(field_to_index["v"]);
                CalcNeutralSources_nvp(v.mass, ee_idx, ni_idx, vi_idx, ei_idx,
                                       nn_idx, vn_idx, en_idx, inarray, ne,
                                       outarray, tmp);
            }
            else
            {
                // Neutral n,p; Ion n,v,p; Ion n,p;
                CalcNeutralSources_np(v.mass, ee_idx, ni_idx, ei_idx, nn_idx,
                                      en_idx, inarray, ne, outarray, tmp);
            }
        }

        else if (v.fields.find(field_to_index["v"]) != v.fields.end())
        {
            int vn_idx = v.fields.at(field_to_index["v"]);
            int vi_idx = this->m_ions[v.ion].fields.at(field_to_index["v"]);
            if (this->m_ions[v.ion].fields.find(field_to_index["e"]) !=
                this->m_ions[v.ion].fields.end())
            {
                // Neutral n,v; Ion n,v,p
                int ei_idx = this->m_ions[v.ion].fields.at(field_to_index["e"]);
                CalcNeutralSources_nv(v.mass, ee_idx, ni_idx, vi_idx, ei_idx,
                                      nn_idx, vn_idx, inarray, ne, outarray,
                                      tmp);
            }
            else
            {
                // Neutral n,v; Ion n,v
                CalcNeutralSources_nv(v.mass, ee_idx, ni_idx, vi_idx, nn_idx,
                                      vn_idx, inarray, ne, outarray, tmp);
            }
        }
        else
        {
            // Neutral n; Ion n; Ion n,v; Ion n,v,p
            CalcNeutralSources_n(v.mass, ee_idx, ni_idx, nn_idx, inarray, ne,
                                 outarray, tmp);
        }
    }
    Vmath::Vadd(this->n_pts, outarray[ee_idx], 1, tmp, 1, outarray[ee_idx], 1);
}

/**
 * @brief Populate rhs array ( @p outarray )
 *
 * @param inarray physical values of all fields
 * @param[out] outarray output array (RHSs of time integration equations)
 */
void ReducedBraginskii::DoOdeImplicitRhs(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    int nvariables = inarray.size();
    int ncoeffs    = m_fields[0]->GetNcoeffs();

    Array<OneD, Array<OneD, NekDouble>> tmpOut(nvariables);
    for (int i = 0; i < nvariables; ++i)
    {
        tmpOut[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    DoOdeRhsCoeff(inarray, tmpOut, time);

    for (int i = 0; i < nvariables; ++i)
    {
        m_fields[i]->BwdTrans(tmpOut[i], outarray[i]);
    }
}

/**
 * @brief Compute the right-hand side.
 */
void ReducedBraginskii::DoOdeRhsCoeff(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{

    int nvariables = inarray.size();
    int nTracePts  = GetTraceTotPoints();
    int ncoeffs    = GetNcoeffs();

    // Store forwards/backwards space along trace space
    Array<OneD, Array<OneD, NekDouble>> Fwd(nvariables);
    Array<OneD, Array<OneD, NekDouble>> Bwd(nvariables);

    for (int i = 0; i < nvariables; ++i)
    {
        Fwd[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
        Bwd[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
        m_indfields[i]->GetFwdBwdTracePhys(inarray[i], Fwd[i], Bwd[i]);
    }

    // Calculate advection
    DoAdvectionCoeff(inarray, outarray, time, Fwd, Bwd);

    // Negate results
    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Neg(ncoeffs, outarray[i], 1);
    }

    // Add diffusion terms

    DoDiffusionCoeff(inarray, outarray, Fwd, Bwd);

    if (this->particles_enabled)
    {
        DoParticlesCoeff(inarray, outarray);
    }

    // Add forcing terms
    for (auto &x : m_forcing)
    {
        x->ApplyCoeff(m_indfields, inarray, outarray, time);
    }
}

/**
 * @brief Compute the advection terms for the right-hand side
 */
void ReducedBraginskii::DoAdvectionCoeff(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
    const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    const Array<OneD, Array<OneD, NekDouble>> &pBwd)
{
    int nvariables = inarray.size();
    Array<OneD, Array<OneD, NekDouble>> advVel(m_spacedim);
    std::dynamic_pointer_cast<SU::AdvectionWeakDG>(m_advection)
        ->AdvectCoeffs(nvariables, m_indfields, advVel, inarray, outarray, time,
                       pFwd, pBwd);
}

void ReducedBraginskii::DoDiffusionCoeff(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray,
    const Array<OneD, const Array<OneD, NekDouble>> &pFwd,
    const Array<OneD, const Array<OneD, NekDouble>> &pBwd)
{
    size_t nvariables = inarray.size();
    size_t ncoeffs    = GetNcoeffs();
    size_t nTracePts  = GetTraceTotPoints();

    Array<OneD, Array<OneD, NekDouble>> outarrayDiff{nvariables};
    for (int i = 0; i < nvariables; ++i)
    {
        outarrayDiff[i] = Array<OneD, NekDouble>{ncoeffs, 0.0};
    }

    // if (m_is_diffIP)
    // {
    //     m_diffusion->DiffuseCoeffs(nvariables, m_fields, inarray,
    //     outarrayDiff,
    //                                m_bndEvaluateTime, pFwd, pBwd);
    //     for (int i = 0; i < nvariables; ++i)
    //     {
    //         Vmath::Vadd(ncoeffs, outarrayDiff[i], 1, outarray[i], 1,
    //                     outarray[i], 1);
    //     }
    // }
    // else
    // {
    ASSERTL1(false, "LDGNS not yet validated for implicit compressible "
                    "flow solver");
    Array<OneD, Array<OneD, NekDouble>> inarrayDiff(nvariables);
    Array<OneD, Array<OneD, NekDouble>> inFwd(nvariables);
    Array<OneD, Array<OneD, NekDouble>> inBwd(nvariables);

    for (int i = 0; i < nvariables; ++i)
    {
        inarrayDiff[i] = Array<OneD, NekDouble>(this->n_pts);
        inFwd[i]       = Array<OneD, NekDouble>(nTracePts);
        inBwd[i]       = Array<OneD, NekDouble>(nTracePts);
    }

    // Extract temperature
    m_varConv->GetElectronTemperature(inarray, inarrayDiff[ee_idx]);
    for (const auto &[s, v] : this->GetIons())
    {
        int ei_idx = v.fields.at(field_to_index["e"]);
        m_varConv->GetIonTemperature(s, v.mass, inarray, inarrayDiff[ei_idx]);
    }

    // Repeat calculation for trace space
    if (pFwd == NullNekDoubleArrayOfArray || pBwd == NullNekDoubleArrayOfArray)
    {
        inFwd = NullNekDoubleArrayOfArray;
        inBwd = NullNekDoubleArrayOfArray;
    }
    else
    {
        m_varConv->GetElectronTemperature(pFwd, inFwd[ee_idx]);
        for (const auto &[s, v] : this->GetIons())
        {
            int ei_idx = v.fields.at(field_to_index["e"]);

            m_varConv->GetIonTemperature(s, v.mass, pFwd, inFwd[ei_idx]);
            m_varConv->GetIonTemperature(s, v.mass, pBwd, inBwd[ei_idx]);
        }
    }

    // Diffusion term in coeff rhs form
    m_diffusion->DiffuseCoeffs(nvariables, m_indfields, inarrayDiff,
                               outarrayDiff, inFwd, inBwd);

    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Vadd(ncoeffs, outarrayDiff[i], 1, outarray[i], 1, outarray[i],
                    1);
    }
    //}
}

/**
 * @brief Compute the advection terms for the right-hand side
 */
void ReducedBraginskii::DoParticlesCoeff(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    int ncoeff = GetNcoeffs();
    Array<OneD, NekDouble> tmp(ncoeff, 0.0);

    // Add contribution to electron energy
    m_indfields[ee_idx]->FwdTrans(this->src_fields[0]->GetPhys(), tmp);
    Vmath::Vadd(this->n_pts, outarray[ee_idx], 1, tmp, 1, outarray[ee_idx], 1);

    for (const auto &[s, v] : this->GetIons())
    {
        int ni_idx = v.fields.at(field_to_index["n"]);
        int vi_idx = v.fields.at(field_to_index["v"]);
        int ei_idx = v.fields.at(field_to_index["e"]);
        //  Add contribution to ion density
        m_indfields[ni_idx]->FwdTrans(
            this->src_fields[ni_src_idx[s]]->GetPhys(), tmp);
        Vmath::Vadd(this->n_pts, outarray[ni_idx], 1, tmp, 1, outarray[ni_idx],
                    1);

        // Add contribution to ion energy
        m_indfields[ei_idx]->FwdTrans(
            this->src_fields[ei_src_idx[s]]->GetPhys(), tmp);
        Vmath::Vadd(this->n_pts, outarray[ei_idx], 1, tmp, 1, outarray[ei_idx],
                    1);

        // Add number density source contribution to ion energy
        Array<OneD, NekDouble> dynamic_energy(this->n_pts);
        m_varConv->GetIonDynamicEnergy(s, v.mass, inarray, dynamic_energy);
        Vmath::Vmul(this->n_pts, dynamic_energy, 1,
                    this->src_fields[ni_src_idx[s]]->GetPhys(), 1,
                    dynamic_energy, 1);
        m_indfields[ei_idx]->FwdTrans(dynamic_energy, tmp);
        Vmath::Vadd(this->n_pts, outarray[ei_idx], 1, tmp, 1, outarray[ei_idx],
                    1);

        Vmath::Zero(this->n_pts, dynamic_energy, 1);
        for (int d = 0; d < m_spacedim; ++d)
        {
            Vmath::Vvtvp(this->n_pts, this->b_unit[d], 1,
                         this->src_fields[vi_src_idx[s] + d]->GetPhys(), 1,
                         dynamic_energy, 1, dynamic_energy, 1);
        }
        m_indfields[vi_idx]->FwdTrans(dynamic_energy, tmp);
        Vmath::Vadd(this->n_pts, outarray[vi_idx], 1, tmp, 1, outarray[vi_idx],
                    1);
    }
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

void ReducedBraginskii::load_params()
{
    PlasmaSystem::load_params();

    NekDouble lambda;

    m_session->LoadParameter("k_ci", this->k_ci, 3.9);
    m_session->LoadParameter("k_ce", this->k_ce, 3.16);
    m_session->LoadParameter("lambda", lambda);

    double scaling_constant =
        this->omega_c * this->mesh_length * this->mesh_length;

    double tau_const = 6.0 * (sqrt(2.0 * pow(M_PI, 3)) / lambda) *
                       constants::epsilon_0 * constants::epsilon_0 * 1e12 *
                       pow(this->Tnorm, 1.5) / (constants::c * this->Nnorm);
    // multiply by sqrt(m in eV) * (T in eV)^1.5 * (n in Nnorm)^-1 for collision
    // time in s

    double t_const = 1.0 / (constants::c * constants::c);
    // multiply by (m in eV)/(B in T)  for gyrotime in s

    k_par = this->k_ce * 1.5 * this->Tnorm * this->Nnorm * tau_const /
            sqrt(constants::m_e);
    // Convert to solver length and time scale
    k_par /= scaling_constant;
    // multiply k_par by Z^-2 n^-1 T^2.5 in solver

    k_perp = 1.5 * this->Tnorm * this->Nnorm * t_const * t_const / tau_const;
    k_perp /= scaling_constant;
    // multiply k_perp by A^0.5 Z^2 n B^-2 T^-0.5 in solver

    k_cross = 3.75 * this->Tnorm * this->Nnorm * t_const / this->Bnorm;
    k_cross /= scaling_constant;

    kappa_i_par = this->k_ci * 1.5 * this->Tnorm * this->Nnorm * tau_const /
                  sqrt(constants::m_p);
    kappa_i_par /= scaling_constant;

    kappa_i_perp = 2 * 1.5 * this->Tnorm * this->Nnorm * t_const * t_const *
                   sqrt(constants::m_p) / tau_const;
    kappa_i_perp /= scaling_constant;

    kappa_i_cross = 3.75 * this->Tnorm * this->Nnorm * t_const / this->Bnorm;
    kappa_i_cross /= scaling_constant;

    kappa_e_par = this->k_ce * 1.5 * this->Tnorm * this->Nnorm * tau_const /
                  sqrt(constants::m_e);
    kappa_e_par /= scaling_constant;

    kappa_e_perp = (sqrt(2.0) + 3.25) * 1.5 * this->Tnorm * this->Nnorm *
                   t_const * t_const * sqrt(constants::m_e) / tau_const;
    kappa_e_perp /= scaling_constant;

    kappa_e_cross = 3.75 * this->Tnorm * this->Nnorm * t_const / this->Bnorm;
    kappa_e_cross /= scaling_constant;
}

void ReducedBraginskii::v_ExtraFldOutput(
    std::vector<Array<OneD, NekDouble>> &fieldcoeffs,
    std::vector<std::string> &variables)
{
    PlasmaSystem::v_ExtraFldOutput(fieldcoeffs, variables);
    const int nCoeffs = m_fields[0]->GetNcoeffs();

    if (this->particles_enabled)
    {
        int cnt = 0;
        for (auto &[k, v] : this->particle_sys->get_species())
        {
            variables.emplace_back(k + "_SOURCE_DENSITY");
            Array<OneD, NekDouble> SrcFwd1(nCoeffs);
            m_fields[0]->FwdTransLocalElmt(this->src_fields[cnt++]->GetPhys(),
                                           SrcFwd1);
            fieldcoeffs.emplace_back(SrcFwd1);

            for (int d = 0; d < this->m_spacedim; ++d)
            {
                variables.emplace_back(k + "_SOURCE_MOMENTUM" +
                                       std::to_string(d));
                Array<OneD, NekDouble> SrcFwd(nCoeffs);
                m_fields[0]->FwdTransLocalElmt(
                    this->src_fields[cnt++]->GetPhys(), SrcFwd);
                fieldcoeffs.emplace_back(SrcFwd);
            }

            variables.emplace_back(k + "_SOURCE_ENERGY");
            Array<OneD, NekDouble> SrcFwd2(nCoeffs);
            m_fields[0]->FwdTransLocalElmt(this->src_fields[cnt++]->GetPhys(),
                                           SrcFwd2);
            fieldcoeffs.emplace_back(SrcFwd2);
        }
    }
}
} // namespace PENKNIFE