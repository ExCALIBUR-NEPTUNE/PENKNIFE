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

    // Create storage for velocities
    int npts       = GetNpoints();
    this->m_kpar   = Array<OneD, NekDouble>(npts, 0.0);
    this->m_kcross = Array<OneD, NekDouble>(npts, 0.0);
    this->m_kperp  = Array<OneD, NekDouble>(npts, 0.0);

    pe_idx = m_indfields.size() - this->n_indep_fields;

    // Parallel velocities
    this->v_e_par = Array<OneD, NekDouble>(npts, 0.0);
    this->v_i_par = std::vector<Array<OneD, NekDouble>>(n_species);

    // Per-field advection velocities
    this->adv_vel =
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(m_indfields.size());

    for (int i = 0; i < this->adv_vel.size(); ++i)
    {
        this->adv_vel[i] = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
    }

    for (int d = 0; d < m_spacedim; ++d)
    {
        this->adv_vel[pe_idx][d] = Array<OneD, NekDouble>(npts, 0.0);
    }

    for (const auto &[s, v] : this->GetSpecies())
    {
        int ni_idx = v.fields.at(field_to_index["n"]);
        int vi_idx = v.fields.at(field_to_index["v"]);
        int pi_idx = v.fields.at(field_to_index["p"]);

        this->v_i_par[s] = Array<OneD, NekDouble>(npts, 0.0);
        for (int d = 0; d < m_graph->GetSpaceDimension(); ++d)
        {
            this->adv_vel[ni_idx][d] = Array<OneD, NekDouble>(npts, 0.0);
            this->adv_vel[vi_idx][d] = Array<OneD, NekDouble>(npts, 0.0);
            this->adv_vel[pi_idx][d] = Array<OneD, NekDouble>(npts, 0.0);
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

    m_ode.DefineOdeRhs(&ReducedBraginskii::DoOdeRhs, this);

    if (this->particles_enabled)
    {
        std::vector<Sym<REAL>> src_syms;
        std::vector<int> src_components;

        for (const auto &[s, v] : this->GetIons())
        {
            ni_src_idx.push_back(s * (2 + m_spacedim));
            pi_src_idx.push_back(1 + s * (2 + m_spacedim));
            vi_src_idx.push_back(2 + s * (2 + m_spacedim));

            this->src_fields.emplace_back(
                MemoryManager<MR::DisContField>::AllocateSharedPtr(
                    *std::dynamic_pointer_cast<MR::DisContField>(m_fields[0])));
            src_syms.push_back(Sym<REAL>(v.name + "_SOURCE_DENSITY"));
            src_components.push_back(0);

            this->src_fields.emplace_back(
                MemoryManager<MR::DisContField>::AllocateSharedPtr(
                    *std::dynamic_pointer_cast<MR::DisContField>(m_fields[0])));

            src_syms.push_back(Sym<REAL>(v.name + "_SOURCE_ENERGY"));
            src_components.push_back(0);

            for (int d = 0; d < this->m_spacedim; ++d)
            {
                this->src_fields.emplace_back(
                    MemoryManager<MR::DisContField>::AllocateSharedPtr(
                        *std::dynamic_pointer_cast<MR::DisContField>(
                            m_fields[0])));
                src_syms.push_back(Sym<REAL>(v.name + "_SOURCE_MOMENTUM"));
                src_components.push_back(d);
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
    int npts      = GetNpoints();
    int nTracePts = GetTraceTotPoints();
    for (int f = 0; f < outarray.size(); ++f)
    {
        Vmath::Zero(npts, outarray[f], 1);
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
    CalcVelocities(inarray, m_fields[0]->GetPhys(), this->adv_vel);

    // Perform advection
    DoAdvection(inarray, outarray, time, Fwd, Bwd);

    m_bndConds->Update(inarray, time);

    // DoExtra(inarray, outarray);

    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Neg(npts, outarray[i], 1);
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
    int nvariables = inarray.size();
    Array<OneD, Array<OneD, NekDouble>> advVel(m_spacedim);

    m_advection->Advect(nvariables, m_indfields, advVel, inarray, outarray,
                        time, pFwd, pBwd);
}

/**
 * @brief Compute the advection terms for the right-hand side
 */
void ReducedBraginskii::DoParticles(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    int npts = GetNpoints();

    // Add contribution to electron energy
    Vmath::Vadd(npts, outarray[pe_idx], 1, this->src_fields[0]->GetPhys(), 1,
                outarray[pe_idx], 1);

    for (const auto &[s, v] : this->GetIons())
    {
        int ni_idx = v.fields.at(field_to_index["n"]);
        int vi_idx = v.fields.at(field_to_index["v"]);
        int pi_idx = v.fields.at(field_to_index["p"]);
        //  Add contribution to ion density
        Vmath::Vadd(npts, outarray[ni_idx], 1,
                    this->src_fields[ni_src_idx[s]]->GetPhys(), 1,
                    outarray[ni_idx], 1);

        // Add contribution to ion energy
        Vmath::Vadd(npts, outarray[pi_idx], 1,
                    this->src_fields[pi_src_idx[s]]->GetPhys(), 1,
                    outarray[pi_idx], 1);
        // Add number density source contribution to ion energy
        Array<OneD, NekDouble> dynamic_energy(npts);
        m_varConv->GetIonDynamicEnergy(s, v.mass, inarray, dynamic_energy);
        Vmath::Vvtvp(npts, dynamic_energy, 1,
                     this->src_fields[ni_src_idx[s]]->GetPhys(), 1,
                     outarray[pi_idx], 1, outarray[pi_idx], 1);

        for (int d = 0; d < m_spacedim; ++d)
        {
            Vmath::Vvtvp(npts, this->b_unit[d], 1,
                         this->src_fields[vi_src_idx[s] + d]->GetPhys(), 1,
                         outarray[vi_idx], 1, outarray[vi_idx], 1);
        }
    }
}

/**
 * @brief Compute the gradient of phi for evaluation at the particle positions.
 */
void ReducedBraginskii::ComputeE()
{
    int npts = GetNpoints();

    Vmath::Neg(npts, this->E[0]->UpdatePhys(), 1);
    Vmath::Neg(npts, this->E[1]->UpdatePhys(), 1);
    Vmath::Neg(npts, this->E[2]->UpdatePhys(), 1);

    this->E[0]->FwdTrans(this->E[0]->GetPhys(), this->E[0]->UpdateCoeffs());
    this->E[1]->FwdTrans(this->E[1]->GetPhys(), this->E[1]->UpdateCoeffs());
    this->E[2]->FwdTrans(this->E[2]->GetPhys(), this->E[2]->UpdateCoeffs());
}

void ReducedBraginskii::CalcVelocities(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    const Array<OneD, NekDouble> &ne,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &adv_vel)
{
    int npts = inarray[0].size();
    for (int f = 0; f < adv_vel.size(); ++f)
    {
        for (int d = 0; d < adv_vel[f].size(); ++d)
        {
            Vmath::Zero(npts, adv_vel[f][d], 1);
        }
    }

    // Zero Electron velocity
    Vmath::Zero(npts, m_fields[1]->UpdatePhys(), 1);
    for (const auto &[s, v] : GetIons())
    {
        int ni_idx = v.fields.at(field_to_index["n"]);
        int vi_idx = v.fields.at(field_to_index["v"]);
        int pi_idx = v.fields.at(field_to_index["p"]);
        // Calculate Ion parallel velocities
        Vmath::Smul(npts, 1.0 / v.mass, inarray[vi_idx], 1, this->v_i_par[s],
                    1);
        Vmath::Svtvp(npts, v.charge, this->v_i_par[s], 1,
                     m_fields[1]->GetPhys(), 1, m_fields[1]->UpdatePhys(), 1);

        Vmath::Vdiv(npts, this->v_i_par[s], 1, inarray[ni_idx], 1,
                    this->v_i_par[s], 1);

        for (int d = 0; d < m_spacedim; ++d)
        {
            Vmath::Vmul(npts, this->b_unit[d], 1, this->v_i_par[s], 1,
                        adv_vel[ni_idx][d], 1);
            Vmath::Vcopy(npts, adv_vel[ni_idx][d], 1, adv_vel[vi_idx][d], 1);
            Vmath::Vcopy(npts, adv_vel[ni_idx][d], 1, adv_vel[pi_idx][d], 1);
        }
    }
    Vmath::Vdiv(npts, m_fields[1]->GetPhys(), 1, ne, 1, this->v_e_par, 1);

    for (int d = 0; d < m_spacedim; ++d)
    {
        Vmath::Vmul(npts, this->b_unit[d], 1, this->v_e_par, 1,
                    adv_vel[pe_idx][d], 1);
    }

    for (const auto &[s, v] : GetNeutrals())
    {
        int nn_idx = v.fields.at(field_to_index["n"]);
        int vn_idx = v.fields.at(field_to_index["v"]);
        int pn_idx = v.fields.at(field_to_index["p"]);
        // Calculate Ion parallel velocities
        Vmath::Smul(npts, 1.0 / v.mass, inarray[vn_idx], 1, this->v_i_par[s],
                    1);
        Vmath::Vdiv(npts, this->v_i_par[s], 1, inarray[nn_idx], 1,
                    this->v_i_par[s], 1);

        for (int d = 0; d < m_spacedim; ++d)
        {
            Vmath::Vmul(npts, this->b_unit[d], 1, this->v_i_par[s], 1,
                        adv_vel[nn_idx][d], 1);
            Vmath::Vcopy(npts, adv_vel[nn_idx][d], 1, adv_vel[vn_idx][d], 1);
            Vmath::Vcopy(npts, adv_vel[nn_idx][d], 1, adv_vel[pn_idx][d], 1);
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
            Vmath::Vvtvp(num_trace_pts, normals[d], 1,
                         this->adv_vel_trace[j][d], 1, this->trace_vel_norm[j],
                         1, this->trace_vel_norm[j], 1);
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
    int npts       = field_vals[0].size();
    int nVariables = field_vals.size();

    for (int d = 0; d < m_spacedim; ++d)
    {
        // Electron Energy Flux
        Vmath::Vmul(npts, field_vals[pe_idx], 1, this->adv_vel[pe_idx][d], 1,
                    fluxes[pe_idx][d], 1);
    }

    for (const auto &[s, v] : this->GetSpecies())
    {
        int ni_idx = v.fields.at(field_to_index["n"]);
        int vi_idx = v.fields.at(field_to_index["v"]);
        int pi_idx = v.fields.at(field_to_index["p"]);

        for (int d = 0; d < m_spacedim; ++d)
        {
            // Ion Density Flux
            Vmath::Vmul(npts, field_vals[ni_idx], 1, this->adv_vel[ni_idx][d],
                        1, fluxes[ni_idx][d], 1);
            // Ion Momentum Flux
            Vmath::Vmul(npts, field_vals[vi_idx], 1, this->adv_vel[vi_idx][d],
                        1, fluxes[vi_idx][d], 1);
            // Ion Energy Flux
            Vmath::Vmul(npts, field_vals[pi_idx], 1, this->adv_vel[pi_idx][d],
                        1, fluxes[pi_idx][d], 1);
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
    int npointsIn  = GetNpoints();
    int npointsOut = npointsIn;
    int nTracePts  = GetTraceTotPoints();

    // this should be preallocated
    Array<OneD, Array<OneD, NekDouble>> outarrayDiff(nvariables);
    for (int i = 0; i < nvariables; ++i)
    {
        outarrayDiff[i] = Array<OneD, NekDouble>(npointsOut, 0.0);
    }

    Array<OneD, Array<OneD, NekDouble>> inarrayDiff(nvariables);
    Array<OneD, Array<OneD, NekDouble>> inFwd(nvariables);
    Array<OneD, Array<OneD, NekDouble>> inBwd(nvariables);

    for (int i = 0; i < nvariables; ++i)
    {
        inarrayDiff[i] = Array<OneD, NekDouble>(npointsIn, 0.0);
        inFwd[i]       = Array<OneD, NekDouble>(nTracePts, 0.0);
        inBwd[i]       = Array<OneD, NekDouble>(nTracePts, 0.0);
    }

    // Extract temperature
    m_varConv->GetElectronTemperature(inarray, inarrayDiff[pe_idx]);

    for (const auto &[s, v] : this->GetIons())
    {
        int pi_idx = v.fields.at(field_to_index["p"]);
        m_varConv->GetIonTemperature(s, v.mass, inarray, inarrayDiff[pi_idx]);
    }

    // Repeat calculation for trace space
    if (pFwd == NullNekDoubleArrayOfArray || pBwd == NullNekDoubleArrayOfArray)
    {
        inFwd = NullNekDoubleArrayOfArray;
        inBwd = NullNekDoubleArrayOfArray;
    }
    else
    {
        m_varConv->GetElectronTemperature(pFwd, inFwd[pe_idx]);
        m_varConv->GetElectronTemperature(pBwd, inBwd[pe_idx]);
        for (const auto &[s, v] : this->GetIons())
        {
            int pi_idx = v.fields.at(field_to_index["p"]);
            m_varConv->GetIonTemperature(s, v.mass, pFwd, inFwd[pi_idx]);
            m_varConv->GetIonTemperature(s, v.mass, pBwd, inBwd[pi_idx]);
        }
    }

    m_diffusion->Diffuse(nvariables, m_difffields, inarrayDiff, outarrayDiff,
                         inFwd, inBwd);

    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Vadd(npointsOut, outarrayDiff[i], 1, outarray[i], 1, outarray[i],
                    1);
    }
}

void ReducedBraginskii::CalcDiffTensor()
{
    int npoints = m_fields[0]->GetNpoints();

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Array<OneD, NekDouble> d(npoints, 0.0);
            for (int k = 0; k < npoints; k++)
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

void ReducedBraginskii::CalcK(const Array<OneD, Array<OneD, NekDouble>> &in_arr,
                              int f)
{
    int npoints = m_fields[0]->GetNpoints();
    double Z    = this->m_ions[f].charge;
    double A    = this->m_ions[f].mass;
    int ni_idx  = this->m_ions[f].fields[field_to_index["n"]];
    auto ne     = this->m_fields[0]->GetPhys();

    for (int p = 0; p < npoints; ++p)
    {
        m_kpar[p] = this->k_ci * this->k_par * pow(in_arr[pe_idx][p], 2.5) /
                    (Z * Z * in_arr[ni_idx][p]);
        m_kperp[p] = this->k_perp * Z * Z * std::sqrt(A) * in_arr[ni_idx][p] /
                     (sqrt(in_arr[pe_idx][p]) * this->mag_B[p]);
        m_kcross[p] =
            this->k_cross * ne[p] * in_arr[pe_idx][p] / (sqrt(this->mag_B[p]));
    }
}

// Ion thermal conductivity
void ReducedBraginskii::CalcKappa(
    const Array<OneD, Array<OneD, NekDouble>> &in_arr, int f)
{
    int npoints = m_fields[0]->GetNpoints();
    double Z    = this->m_ions[f].charge;
    double A    = this->m_ions[f].mass;
    int pi_idx  = this->m_ions[f].fields[field_to_index["p"]];
    int ni_idx  = this->m_ions[f].fields[field_to_index["n"]];

    Array<OneD, NekDouble> tmp(npoints, 0.0);
    for (const auto &[s2, v2] : this->GetIons())
    {
        double Z2   = this->m_ions[s2].charge;
        double A2   = this->m_ions[s2].mass;
        int ni_idx2 = this->m_ions[s2].fields.at(field_to_index["n"]);

        for (int p = 0; p < npoints; ++p)
        {
            tmp[p] += Z2 * Z2 * sqrt(A2 / (A + A2)) * in_arr[ni_idx2][p];
        }
    }
    for (int p = 0; p < npoints; ++p)
    {
        this->m_kpar[p] = this->kappa_i_par * in_arr[ni_idx][p] *
                          (in_arr[pi_idx][p], 2.5) / (sqrt(A) * Z * Z * tmp[p]);
        this->m_kperp[p] = this->kappa_i_perp * sqrt(A) * tmp[p] *
                           in_arr[ni_idx][p] /
                           (this->mag_B[p] * sqrt(in_arr[pi_idx][p]));
        this->m_kcross[p] = this->kappa_i_cross * in_arr[ni_idx][p] *
                            in_arr[pi_idx][p] / (Z * sqrt(this->mag_B[p]));
    }
}

// Electron thermal conductivity
void ReducedBraginskii::CalcKappa(
    const Array<OneD, Array<OneD, NekDouble>> &in_arr)
{
    int npoints = m_fields[0]->GetNpoints();
    auto ne     = this->m_fields[0]->GetPhys();
    for (int p = 0; p < npoints; ++p)
    {
        this->m_kpar[p]  = this->kappa_e_par * pow(in_arr[pe_idx][p], 2.5);
        this->m_kperp[p] = this->kappa_e_perp * ne[p] * ne[p] /
                           (this->mag_B[p] * sqrt(in_arr[pe_idx][p]));
        this->m_kcross[p] = this->kappa_e_cross * ne[p] * in_arr[pe_idx][p] /
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
    unsigned int nPts = qfield[0][0].size();

    for (const auto &[s, v] : this->GetIons())
    {
        int ni_idx = v.fields.at(field_to_index["n"]);
        int pi_idx = v.fields.at(field_to_index["p"]);

        CalcK(in_arr, s);
        CalcDiffTensor();

        for (unsigned int j = 0; j < nDim; ++j)
        {
            Vmath::Vmul(nPts, m_D[vc[j][0]].GetValue(), 1, qfield[0][ni_idx], 1,
                        fluxes[j][ni_idx], 1);
            for (unsigned int k = 1; k < nDim; ++k)
            {
                Vmath::Vvtvp(nPts, m_D[vc[j][k]].GetValue(), 1,
                             qfield[k][ni_idx], 1, fluxes[j][ni_idx], 1,
                             fluxes[j][ni_idx], 1);
            }
        }

        if (nDim == 3)
        {
            Vmath::Vvtvvtm(nPts, b_unit[1], 1, qfield[2][pi_idx], 1, b_unit[2],
                           1, qfield[1][pi_idx], 1, fluxes[0][pi_idx], 1);
            Vmath::Vvtvvtm(nPts, b_unit[2], 1, qfield[0][pi_idx], 1, b_unit[0],
                           1, qfield[2][pi_idx], 1, fluxes[1][pi_idx], 1);
            Vmath::Vvtvvtm(nPts, b_unit[0], 1, qfield[1][pi_idx], 1, b_unit[1],
                           1, qfield[0][pi_idx], 1, fluxes[2][pi_idx], 1);
        }
        else
        {
            Vmath::Vmul(nPts, b_unit[2], 1, qfield[1][pi_idx], 1,
                        fluxes[0][pi_idx], 1);
            Vmath::Neg(nPts, fluxes[0][pi_idx], 1);
            Vmath::Vmul(nPts, b_unit[2], 1, qfield[0][pi_idx], 1,
                        fluxes[1][pi_idx], 1);
        }

        CalcKappa(in_arr, s);
        CalcDiffTensor();
        for (unsigned int j = 0; j < nDim; ++j)
        {
            Vmath::Vmul(nPts, m_kcross, 1, fluxes[j][pi_idx], 1,
                        fluxes[j][pi_idx], 1);
            // Calc diffusion of n with D tensor and n field
            for (unsigned int k = 0; k < nDim; ++k)
            {
                Vmath::Vvtvp(nPts, m_D[vc[j][k]].GetValue(), 1,
                             qfield[k][pi_idx], 1, fluxes[j][pi_idx], 1,
                             fluxes[j][pi_idx], 1);
            }
        }
    }
}

void ReducedBraginskii::CalcNeutralSources_nvp(
    const double m, int pe_idx, int ni_idx, int vi_idx, int pi_idx, int nn_idx,
    int vn_idx, int pn_idx, const Array<OneD, Array<OneD, NekDouble>> &inarray,
    const Array<OneD, NekDouble> &ne,
    Array<OneD, Array<OneD, NekDouble>> &outarray, Array<OneD, NekDouble> &Spe)
{
    unsigned int nPts = ne.size();

    const Array<OneD, NekDouble> &pe = inarray[pe_idx];
    const Array<OneD, NekDouble> &nn = inarray[nn_idx];
    const Array<OneD, NekDouble> &vn = inarray[vn_idx];
    const Array<OneD, NekDouble> &pn = inarray[pn_idx];
    const Array<OneD, NekDouble> &ni = inarray[ni_idx];
    const Array<OneD, NekDouble> &vi = inarray[vi_idx];
    const Array<OneD, NekDouble> &pi = inarray[pi_idx];

    for (int p = 0; p < nPts; ++p)
    {
        double exponent = 13.6 / pe[p];
        double krec     = 0.7e-19 * std::sqrt(exponent);
        double kIZ      = (2e-13 / (6 + 1.0 / exponent)) *
                     std::sqrt(1.0 / exponent) * std::exp(-exponent);
        double kCX = 3.2e-15 * std::sqrt(pi[p] / 0.026);

        double SN = -kIZ * ne[p] * nn[p] + krec * ne[p] * ni[p];

        double SGN = -kIZ * ne[p] * vn[p] + krec * ne[p] * vi[p] +
                     kCX * nn[p] * vi[p] - kCX * ni[p] * vn[p];

        double vdiffsq =
            (vi[p] / ni[p] - vn[p] / nn[p]) * (vi[p] / ni[p] - vn[p] / nn[p]);
        double SP = -kIZ * ne[p] * pn[p] + krec * ne[p] * pi[p] +
                    kCX * nn[p] * pi[p] - kCX * ni[p] * pn[p] +
                    (m * ni[p] / 3) * (krec * ne[p] + kCX * nn[p]) * vdiffsq;

        double SPI = ni[p] * nn[p] * (kIZ + kCX) *
                     ((pn[p] / nn[p] - pi[p] / ni[p]) + (m / 3) * vdiffsq);

        double Wiz  = 1.0; // TODO get from AMJUEL
        double Wrec = 1.0;

        Spe[p] -= ne[p] * pe[p] * (kIZ * nn[p] - krec * ni[p]) +
                  ne[p] * (2.0 / 3.0) * (Wiz * nn[p] + Wrec * ni[p]);

        outarray[nn_idx][p] += SN;
        outarray[ni_idx][p] -= SN;
        outarray[vn_idx][p] += SGN;
        outarray[vi_idx][p] -= SGN + SN * vi[p] / ni[p];
        outarray[pn_idx][p] += SP;
        outarray[pi_idx][p] += SPI;
    }
}

void ReducedBraginskii::CalcNeutralSources_nv(
    const double m, int pe_idx, int ni_idx, int vi_idx, int nn_idx, int vn_idx,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    const Array<OneD, NekDouble> &ne,
    Array<OneD, Array<OneD, NekDouble>> &outarray, Array<OneD, NekDouble> &Spe)
{
    unsigned int nPts = ne.size();

    const Array<OneD, NekDouble> &pe = inarray[pe_idx];
    const Array<OneD, NekDouble> &nn = inarray[nn_idx];
    const Array<OneD, NekDouble> &vn = inarray[vn_idx];
    const Array<OneD, NekDouble> &ni = inarray[ni_idx];
    const Array<OneD, NekDouble> &vi = inarray[vi_idx];

    for (int p = 0; p < nPts; ++p)
    {
        double exponent = 13.6 / pe[p];
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

        Spe[p] -= ne[p] * pe[p] * (kIZ * nn[p] - krec * ni[p]) +
                  ne[p] * (2.0 / 3.0) * (Wiz * nn[p] + Wrec * ni[p]);

        outarray[nn_idx][p] += SN;
        outarray[ni_idx][p] -= SN;
        outarray[vn_idx][p] += SGN;
        outarray[vi_idx][p] -= SGN + SN * vi[p] / ni[p];
    }
}

void ReducedBraginskii::CalcNeutralSources_nv(
    const double m, int pe_idx, int ni_idx, int vi_idx, int pi_idx, int nn_idx,
    int vn_idx, const Array<OneD, Array<OneD, NekDouble>> &inarray,
    const Array<OneD, NekDouble> &ne,
    Array<OneD, Array<OneD, NekDouble>> &outarray, Array<OneD, NekDouble> &Spe)
{
    unsigned int nPts = ne.size();

    const Array<OneD, NekDouble> &pe = inarray[pe_idx];
    const Array<OneD, NekDouble> &nn = inarray[nn_idx];
    const Array<OneD, NekDouble> &vn = inarray[vn_idx];
    const Array<OneD, NekDouble> &ni = inarray[ni_idx];
    const Array<OneD, NekDouble> &vi = inarray[vi_idx];
    const Array<OneD, NekDouble> &pi = inarray[pi_idx];

    for (int p = 0; p < nPts; ++p)
    {
        double exponent = 13.6 / pe[p];
        double krec     = 0.7e-19 * std::sqrt(exponent);
        double kIZ      = (2e-13 / (6 + 1.0 / exponent)) *
                     std::sqrt(1.0 / exponent) * std::exp(-exponent);
        double kCX = 3.2e-15 * std::sqrt(pi[p] / 0.026);

        double SN = -kIZ * ne[p] * nn[p] + krec * ne[p] * ni[p];

        double SGN = -kIZ * ne[p] * vn[p] + krec * ne[p] * vi[p] +
                     kCX * nn[p] * vi[p] - kCX * ni[p] * vn[p];

        double vdiffsq =
            (vi[p] / ni[p] - vn[p] / nn[p]) * (vi[p] / ni[p] - vn[p] / nn[p]);

        double SPI = ni[p] * nn[p] * (kIZ + kCX) * (m / 3) * vdiffsq;

        double Wiz  = 1.0; // TODO get from AMJUEL
        double Wrec = 1.0;

        Spe[p] -= ne[p] * pe[p] * (kIZ * nn[p] - krec * ni[p]) +
                  ne[p] * (2.0 / 3.0) * (Wiz * nn[p] + Wrec * ni[p]);

        outarray[nn_idx][p] += SN;
        outarray[ni_idx][p] -= SN;
        outarray[vn_idx][p] += SGN;
        outarray[vi_idx][p] -= SGN + SN * vi[p] / ni[p];
        outarray[pi_idx][p] += SPI;
    }
}

void ReducedBraginskii::CalcNeutralSources_np(
    const double m, int pe_idx, int ni_idx, int pi_idx, int nn_idx, int pn_idx,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    const Array<OneD, NekDouble> &ne,
    Array<OneD, Array<OneD, NekDouble>> &outarray, Array<OneD, NekDouble> &Spe)
{
    unsigned int nPts = ne.size();

    const Array<OneD, NekDouble> &pe = inarray[pe_idx];
    const Array<OneD, NekDouble> &nn = inarray[nn_idx];
    const Array<OneD, NekDouble> &pn = inarray[pn_idx];
    const Array<OneD, NekDouble> &ni = inarray[ni_idx];
    const Array<OneD, NekDouble> &pi = inarray[pi_idx];

    for (int p = 0; p < nPts; ++p)
    {
        double exponent = 13.6 / pe[p];
        double krec     = 0.7e-19 * std::sqrt(exponent);
        double kIZ      = (2e-13 / (6 + 1.0 / exponent)) *
                     std::sqrt(1.0 / exponent) * std::exp(-exponent);
        double kCX = 3.2e-15 * std::sqrt(pi[p] / 0.026);

        double SN = -kIZ * ne[p] * nn[p] + krec * ne[p] * ni[p];

        double SP = -kIZ * ne[p] * pn[p] + krec * ne[p] * pi[p] +
                    kCX * nn[p] * pi[p] - kCX * ni[p] * pn[p];

        double SPI =
            ni[p] * nn[p] * (kIZ + kCX) * (pn[p] / nn[p] - pi[p] / ni[p]);

        double Wiz  = 1.0; // TODO get from AMJUEL
        double Wrec = 1.0;

        Spe[p] -= ne[p] * pe[p] * (kIZ * nn[p] - krec * ni[p]) +
                  ne[p] * (2.0 / 3.0) * (Wiz * nn[p] + Wrec * ni[p]);

        outarray[nn_idx][p] += SN;
        outarray[ni_idx][p] -= SN;
        outarray[pn_idx][p] += SP;
        outarray[pi_idx][p] += SPI;
    }
}

void ReducedBraginskii::CalcNeutralSources_n(
    const double m, int pe_idx, int ni_idx, int nn_idx,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    const Array<OneD, NekDouble> &ne,
    Array<OneD, Array<OneD, NekDouble>> &outarray, Array<OneD, NekDouble> &Spe)
{
    unsigned int nPts = ne.size();

    const Array<OneD, NekDouble> &pe = inarray[pe_idx];
    const Array<OneD, NekDouble> &nn = inarray[nn_idx];
    const Array<OneD, NekDouble> &ni = inarray[ni_idx];

    for (int p = 0; p < nPts; ++p)
    {
        double exponent = 13.6 / pe[p];
        double krec     = 0.7e-19 * std::sqrt(exponent);
        double kIZ      = (2e-13 / (6 + 1.0 / exponent)) *
                     std::sqrt(1.0 / exponent) * std::exp(-exponent);

        double SN = -kIZ * ne[p] * nn[p] + krec * ne[p] * ni[p];

        double Wiz  = 1.0; // TODO get from AMJUEL
        double Wrec = 1.0;

        Spe[p] -= ne[p] * pe[p] * (kIZ * nn[p] - krec * ni[p]) +
                  ne[p] * (2.0 / 3.0) * (Wiz * nn[p] + Wrec * ni[p]);

        outarray[nn_idx][p] += SN;
        outarray[ni_idx][p] -= SN;
    }
}

void ReducedBraginskii::AddNeutralSources(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    unsigned int nPts = inarray[0].size();
    auto ne           = m_fields[0]->GetPhys();
    Array<OneD, NekDouble> tmp(nPts, 0.0);

    for (const auto &[s, v] : this->GetNeutrals())
    {
        int nn_idx = v.fields.at(field_to_index["n"]);
        int ni_idx = this->m_ions[v.ion].fields.at(field_to_index["n"]);

        if (v.fields.find(field_to_index["p"]) != v.fields.end())
        {
            int pn_idx = v.fields.at(field_to_index["p"]);
            int pi_idx = this->m_ions[v.ion].fields.at(field_to_index["p"]);

            if (v.fields.find(field_to_index["v"]) != v.fields.end())
            {
                // Neutral n,v,p; Ion n,v,p
                int vn_idx = v.fields.at(field_to_index["v"]);
                int vi_idx = this->m_ions[v.ion].fields.at(field_to_index["v"]);
                CalcNeutralSources_nvp(v.mass, pe_idx, ni_idx, vi_idx, pi_idx,
                                       nn_idx, vn_idx, pn_idx, inarray, ne,
                                       outarray, tmp);
            }
            else
            {
                // Neutral n,p; Ion n,v,p; Ion n,p;
                CalcNeutralSources_np(v.mass, pe_idx, ni_idx, pi_idx, nn_idx,
                                      pn_idx, inarray, ne, outarray, tmp);
            }
        }

        else if (v.fields.find(field_to_index["v"]) != v.fields.end())
        {
            int vn_idx = v.fields.at(field_to_index["v"]);
            int vi_idx = this->m_ions[v.ion].fields.at(field_to_index["v"]);
            if (this->m_ions[v.ion].fields.find(field_to_index["p"]) !=
                this->m_ions[v.ion].fields.end())
            {
                // Neutral n,v; Ion n,v,p
                int pi_idx = this->m_ions[v.ion].fields.at(field_to_index["p"]);
                CalcNeutralSources_nv(v.mass, pe_idx, ni_idx, vi_idx, pi_idx,
                                      nn_idx, vn_idx, inarray, ne, outarray,
                                      tmp);
            }
            else
            {
                // Neutral n,v; Ion n,v
                CalcNeutralSources_nv(v.mass, pe_idx, ni_idx, vi_idx, nn_idx,
                                      vn_idx, inarray, ne, outarray, tmp);
            }
        }
        else
        {
            // Neutral n; Ion n; Ion n,v; Ion n,v,p
            CalcNeutralSources_n(v.mass, pe_idx, ni_idx, nn_idx, inarray, ne,
                                 outarray, tmp);
        }
    }
    Vmath::Vadd(nPts, outarray[pe_idx], 1, tmp, 1, outarray[pe_idx], 1);
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
    size_t npoints    = GetNpoints();
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
    Array<OneD, Array<OneD, NekDouble>> inarrayDiff{nvariables};
    Array<OneD, Array<OneD, NekDouble>> inFwd{nvariables};
    Array<OneD, Array<OneD, NekDouble>> inBwd{nvariables};

    for (int i = 0; i < nvariables; ++i)
    {
        inarrayDiff[i] = Array<OneD, NekDouble>{npoints};
        inFwd[i]       = Array<OneD, NekDouble>{nTracePts};
        inBwd[i]       = Array<OneD, NekDouble>{nTracePts};
    }

    // Extract temperature
    m_varConv->GetElectronTemperature(inarray, inarrayDiff[pe_idx]);
    for (const auto &[s, v] : this->GetIons())
    {
        int pi_idx = v.fields.at(field_to_index["p"]);
        m_varConv->GetIonTemperature(s, v.mass, inarray, inarrayDiff[pi_idx]);
    }

    // Repeat calculation for trace space
    if (pFwd == NullNekDoubleArrayOfArray || pBwd == NullNekDoubleArrayOfArray)
    {
        inFwd = NullNekDoubleArrayOfArray;
        inBwd = NullNekDoubleArrayOfArray;
    }
    else
    {
        m_varConv->GetElectronTemperature(pFwd, inFwd[pe_idx]);
        for (const auto &[s, v] : this->GetIons())
        {
            int pi_idx = v.fields.at(field_to_index["p"]);

            m_varConv->GetIonTemperature(s, v.mass, pFwd, inFwd[pi_idx]);
            m_varConv->GetIonTemperature(s, v.mass, pBwd, inBwd[pi_idx]);
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
    int npts   = GetNpoints();
    int ncoeff = GetNcoeffs();
    Array<OneD, NekDouble> tmp(ncoeff, 0.0);

    // Add contribution to electron energy
    m_indfields[pe_idx]->FwdTrans(this->src_fields[0]->GetPhys(), tmp);
    Vmath::Vadd(npts, outarray[pe_idx], 1, tmp, 1, outarray[pe_idx], 1);

    for (const auto &[s, v] : this->GetIons())
    {
        int ni_idx = v.fields.at(field_to_index["n"]);
        int vi_idx = v.fields.at(field_to_index["v"]);
        int pi_idx = v.fields.at(field_to_index["p"]);
        //  Add contribution to ion density
        m_indfields[ni_idx]->FwdTrans(
            this->src_fields[ni_src_idx[s]]->GetPhys(), tmp);
        Vmath::Vadd(npts, outarray[ni_idx], 1, tmp, 1, outarray[ni_idx], 1);

        // Add contribution to ion energy
        m_indfields[pi_idx]->FwdTrans(
            this->src_fields[pi_src_idx[s]]->GetPhys(), tmp);
        Vmath::Vadd(npts, outarray[pi_idx], 1, tmp, 1, outarray[pi_idx], 1);

        // Add number density source contribution to ion energy
        Array<OneD, NekDouble> dynamic_energy(npts);
        m_varConv->GetIonDynamicEnergy(s, v.mass, inarray, dynamic_energy);
        Vmath::Vmul(npts, dynamic_energy, 1,
                    this->src_fields[ni_src_idx[s]]->GetPhys(), 1,
                    dynamic_energy, 1);
        m_indfields[pi_idx]->FwdTrans(dynamic_energy, tmp);
        Vmath::Vadd(npts, outarray[pi_idx], 1, tmp, 1, outarray[pi_idx], 1);

        Vmath::Zero(npts, dynamic_energy, 1);
        for (int d = 0; d < m_spacedim; ++d)
        {
            Vmath::Vvtvp(npts, this->b_unit[d], 1,
                         this->src_fields[vi_src_idx[s] + d]->GetPhys(), 1,
                         dynamic_energy, 1, dynamic_energy, 1);
        }
        m_indfields[vi_idx]->FwdTrans(dynamic_energy, tmp);
        Vmath::Vadd(npts, outarray[vi_idx], 1, tmp, 1, outarray[vi_idx], 1);
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
    // Type of advection to use. Default is DG.
    m_session->LoadSolverInfo("AdvectionType", this->adv_type, "WeakDG");
    // Type of Riemann solver to use. Default = "Upwind"
    m_session->LoadSolverInfo("UpwindType", this->riemann_solver_type,
                              "MultiFieldUpwind");
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
    const int nPhys   = m_fields[0]->GetNpoints();
    const int nCoeffs = m_fields[0]->GetNcoeffs();

    if (this->particles_enabled)
    {
        int i = 0;
        for (auto &[k, v] : this->particle_sys->get_species())
        {
            variables.push_back(k + "_SOURCE_DENSITY");
            Array<OneD, NekDouble> SrcFwd1(nCoeffs);
            m_fields[0]->FwdTransLocalElmt(this->src_fields[i]->GetPhys(),
                                           SrcFwd1);
            fieldcoeffs.push_back(SrcFwd1);

            variables.push_back(k + "_SOURCE_ENERGY");
            Array<OneD, NekDouble> SrcFwd2(nCoeffs);
            m_fields[0]->FwdTransLocalElmt(this->src_fields[i + 1]->GetPhys(),
                                           SrcFwd2);
            fieldcoeffs.push_back(SrcFwd2);
            i += (2 + m_spacedim);
        }
    }
}
} // namespace PENKNIFE