#include "ElectrostaticTurbulence.hpp"
#include "../RiemannSolvers/PlasmaSolver.hpp"
#include <SolverUtils/Advection/AdvectionNonConservative.h>
#include <SolverUtils/Advection/AdvectionWeakDG.h>

namespace PENKNIFE
{

/// Name of class
static std::string class_name;
std::string ElectrostaticTurbulence::class_name =
    SU::GetEquationSystemFactory().RegisterCreatorFunction(
        "ElectrostaticTurbulence", ElectrostaticTurbulence::create,
        "Solves electrostatic turbulence with anisotropic diffusion");
/**
 * @brief Creates an instance of this class.
 */
static SU::EquationSystemSharedPtr create(
    const LU::SessionReaderSharedPtr &session,
    const SD::MeshGraphSharedPtr &graph)
{
    SU::EquationSystemSharedPtr p =
        MemoryManager<ElectrostaticTurbulence>::AllocateSharedPtr(session,
                                                                  graph);
    p->InitObject();
    return p;
}

ElectrostaticTurbulence::ElectrostaticTurbulence(
    const LU::SessionReaderSharedPtr &session,
    const SD::MeshGraphSharedPtr &graph)
    : PlasmaSystem(session, graph)
{
    this->n_indep_fields = 3; // p_e, w, phi
}

void ElectrostaticTurbulence::v_InitObject(bool DeclareFields)
{
    PlasmaSystem::v_InitObject(DeclareFields);

    std::string diffName;
    m_session->LoadSolverInfo("DiffusionType", diffName, "LDG");
    m_diffusion =
        SolverUtils::GetDiffusionFactory().CreateInstance(diffName, diffName);
    m_diffusion->SetFluxVector(&ElectrostaticTurbulence::GetFluxVectorDiff,
                               this);

    // workaround for bug in DiffusionLDG
    m_difffields = Array<OneD, MR::ExpListSharedPtr>(m_indfields.size() - 1);
    for (int f = 0; f < m_difffields.size(); ++f)
    {
        m_difffields[f] = m_indfields[f];
    }

    m_diffusion->InitObject(m_session, m_difffields);

    // Create storage for velocities
    int npts  = GetNpoints();
    pe_idx    = m_indfields.size() - this->n_indep_fields;
    omega_idx = m_indfields.size() - this->n_indep_fields + 1;
    phi_idx   = m_indfields.size() - this->n_indep_fields + 2;

    // ExB velocity
    this->v_ExB = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);

    // Drift velocities
    this->v_de = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
    this->v_di = std::vector<Array<OneD, Array<OneD, NekDouble>>>(n_species);

    // Per-field advection velocities (phi not advected, omega is calculates
    // separately)
    this->adv_vel = Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(
        m_indfields.size() - 1);

    for (int i = 0; i < this->adv_vel.size(); ++i)
    {
        this->adv_vel[i] = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
    }
    this->omega_flux = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);

    for (int d = 0; d < m_spacedim; ++d)
    {
        this->v_ExB[d]              = Array<OneD, NekDouble>(npts, 0.0);
        this->v_de[d]               = Array<OneD, NekDouble>(npts, 0.0);
        this->adv_vel[pe_idx][d]    = Array<OneD, NekDouble>(npts, 0.0);
        this->adv_vel[omega_idx][d] = Array<OneD, NekDouble>(npts, 0.0);
        this->omega_flux[d]         = Array<OneD, NekDouble>(npts, 0.0);
    }

    for (const auto &[s, v] : this->GetSpecies())
    {
        if (v.fields.find(field_to_index["v"]) != v.fields.end())
        {
            int ni_idx = v.fields.at(field_to_index["n"]);
            int vi_idx = v.fields.at(field_to_index["v"]);
            this->advected_fields.push_back(ni_idx);
            this->advected_fields.push_back(vi_idx);

            this->v_di[s] = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);

            for (int d = 0; d < m_graph->GetSpaceDimension(); ++d)
            {
                this->v_di[s][d] = Array<OneD, NekDouble>(npts, 0.0);
            }

            if (v.fields.find(field_to_index["p"]) != v.fields.end())
            {
                int pi_idx = v.fields.at(field_to_index["p"]);
                this->advected_fields.push_back(pi_idx);
            }
        }
    }

    this->advected_fields.push_back(pe_idx);
    this->advected_fields.push_back(omega_idx);

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
            this->adv_vel[i][d] = Array<OneD, NekDouble>(npts, 0.0);
        }
    }

    // Since we are starting from a setup where each field is defined to be a
    // discontinuous field (and thus support DG), the first thing we do is to
    // recreate the phi field so that it is continuous, in order to support the
    // Poisson solve. Note that you can still perform a Poisson solve using a
    // discontinuous field, which is done via the hybridisable discontinuous
    // Galerkin (HDG) approach.

    this->phi   = m_fields[4];
    m_fields[4] = MemoryManager<MR::ContField>::AllocateSharedPtr(
        m_session, m_graph, "phi", true, false);
    m_indfields[phi_idx] = m_fields[4];

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
        this->omega_flux_norm = Array<OneD, NekDouble>(nTrace, 0.0);
        this->omega_flux_trace =
            Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
        for (int d = 0; d < m_spacedim; ++d)
        {
            this->omega_flux_trace[d] = Array<OneD, NekDouble>(nTrace, 0.0);
        }
    }

    // Create Riemann solver and set normal velocity
    // callback functions
    this->riemann_solver = SU::GetRiemannSolverFactory().CreateInstance(
        this->riemann_solver_type, m_session);
    auto t = std::dynamic_pointer_cast<PlasmaSolver>(this->riemann_solver);
    t->omega_idx = omega_idx;

    this->riemann_solver->SetVector(
        "Vn", &ElectrostaticTurbulence::GetAdvVelNorm, this);
    this->riemann_solver->SetScalar(
        "wf", &ElectrostaticTurbulence::GetOmegaFlux, this);

    // Setup advection object
    m_session->LoadSolverInfo("AdvectionType", this->adv_type, "WeakDG");
    m_advection = SU::GetAdvectionFactory().CreateInstance(this->adv_type,
                                                           this->adv_type);
    m_advection->SetFluxVector(&ElectrostaticTurbulence::GetFluxVector, this);
    m_advection->SetRiemannSolver(this->riemann_solver);
    m_advection->InitObject(m_session, m_indfields);

    m_ode.DefineOdeRhs(&ElectrostaticTurbulence::DoOdeRhs, this);

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
            if (v.fields.find(field_to_index["p"]) != v.fields.end())
            {
                this->src_fields.emplace_back(
                    MemoryManager<MR::DisContField>::AllocateSharedPtr(
                        *std::dynamic_pointer_cast<MR::DisContField>(
                            m_fields[0])));

                src_syms.push_back(Sym<REAL>(v.name + "_SOURCE_ENERGY"));
                src_components.push_back(0);
                pi_src_idx.push_back(cnt++);
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

bool ElectrostaticTurbulence::v_PostIntegrate(int step)
{
    m_fields[0]->FwdTrans(m_fields[0]->GetPhys(), m_fields[0]->UpdateCoeffs());
    m_fields[1]->FwdTrans(m_fields[1]->GetPhys(), m_fields[1]->UpdateCoeffs());

    // Writes a step of the particle trajectory.

    return PlasmaSystem::v_PostIntegrate(step);
}

/**
 * @brief Populate rhs array ( @p outarray )
 *
 * @param inarray physical values of all fields
 * @param[out] outarray output array (RHSs of time integration equations)
 */
void ElectrostaticTurbulence::DoOdeRhs(
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
    Array<OneD, Array<OneD, NekDouble>> Fwd(nvariables - 1);
    Array<OneD, Array<OneD, NekDouble>> Bwd(nvariables - 1);

    for (int i = 0; i < nvariables - 1; ++i)
    {
        Fwd[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
        Bwd[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
        m_indfields[i]->GetFwdBwdTracePhys(inarray[i], Fwd[i], Bwd[i]);
    }

    // Helmholtz Solve for electrostatic potential
    SolvePhi(inarray, m_fields[0]->GetPhys());

    // Calculate E
    ComputeE();
    // Calculate ExB, parallel and diamagnetic velocities
    ComputevExB();
    CalcVelocities(inarray);
    // AddDriftVelocities(inarray, m_fields[0]->GetPhys(), this->adv_vel);

    // CalcOmegaFlux(inarray, this->omega_flux);

    // Perform advection
    DoAdvection(inarray, outarray, time, Fwd, Bwd);

    m_bndConds->Update(inarray, time);

    // DoExtra(inarray, outarray);

    for (int i = 0; i < nvariables - 1; ++i)
    {
        Vmath::Neg(npts, outarray[i], 1);
    }

    CalcKappaTensor();

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

void ElectrostaticTurbulence::DoExtra(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    int npts = GetNpoints();
    Array<OneD, NekDouble> grad(npts, 0.0);
    Array<OneD, NekDouble> extra(npts, 0.0);
    Array<OneD, NekDouble> vel(npts, 0.0);
    // Add extra terms to electron energy

    for (int d = 0; d < m_spacedim; ++d)
    {
        Vmath::Vvtvp(npts, this->b_unit[d], 1, this->v_e_par, 1, this->v_ExB[d],
                     1, vel, 1);
        m_fields[0]->PhysDeriv(d, vel, grad);
        Vmath::Vadd(npts, grad, 1, extra, 1, extra, 1);
    }
    Vmath::Vmul(npts, inarray[pe_idx], 1, extra, 1, extra, 1);
    Vmath::Smul(npts, 2.0 / 3.0, extra, 1, extra, 1);
    Vmath::Vadd(npts, extra, 1, outarray[pe_idx], 1, outarray[pe_idx], 1);
    Vmath::Zero(npts, extra, 1);

    for (const auto &[s, v] : this->GetIons())
    {
        int ni_idx    = v.fields.at(field_to_index["n"]);
        int vi_idx    = v.fields.at(field_to_index["v"]);
        int pi_idx    = v.fields.at(field_to_index["e"]);
        double charge = v.charge;
        // Add extra terms to ion momentum

        for (int d = 0; d < m_spacedim; ++d)
        {
            m_fields[0]->PhysDeriv(d, inarray[pi_idx], grad);
            Vmath::Vvtvp(npts, this->b_unit[d], 1, grad, 1, extra, 1, extra, 1);
        }
        Vmath::Smul(npts, 2.0 / 3.0, extra, 1, extra, 1);
        Vmath::Vadd(npts, extra, 1, outarray[vi_idx], 1, outarray[vi_idx], 1);
        Vmath::Zero(npts, extra, 1);
        for (int d = 0; d < m_spacedim; ++d)
        {
            Vmath::Vvtvp(npts, this->b_unit[d], 1, this->E[d]->GetPhys(), 1,
                         extra, 1, extra, 1);
        }
        Vmath::Vmul(npts, inarray[ni_idx], 1, extra, 1, extra, 1);
        Vmath::Smul(npts, charge, extra, 1, extra, 1);
        Vmath::Vadd(npts, extra, 1, outarray[vi_idx], 1, outarray[vi_idx], 1);
        Vmath::Zero(npts, extra, 1);

        // Add extra terms to ion energy
        for (int d = 0; d < m_spacedim; ++d)
        {
            Vmath::Vvtvp(npts, this->b_unit[d], 1, this->v_i_par[s], 1,
                         this->v_ExB[d], 1, vel, 1);
            m_fields[0]->PhysDeriv(d, vel, grad);
            Vmath::Vadd(npts, grad, 1, extra, 1, extra, 1);
        }
        Vmath::Smul(npts, 2.0 / 3.0, extra, 1, extra, 1);
        Vmath::Vadd(npts, extra, 1, outarray[pi_idx], 1, outarray[pi_idx], 1);
    }
}

/**
 * @brief Compute the advection terms for the right-hand side
 */
void ElectrostaticTurbulence::DoAdvection(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
    const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    const Array<OneD, Array<OneD, NekDouble>> &pBwd)
{
    int nvariables = this->adv_vel.size();
    int npointsIn  = GetNpoints();
    int npointsOut = npointsIn;
    int nTracePts  = GetTraceTotPoints();

    Array<OneD, Array<OneD, NekDouble>> advVel(m_spacedim);

    Array<OneD, Array<OneD, NekDouble>> outarrayAdv(nvariables);
    for (int i = 0; i < nvariables; ++i)
    {
        outarrayAdv[i] = Array<OneD, NekDouble>(npointsOut, 0.0);
    }

    Array<OneD, Array<OneD, NekDouble>> inarrayAdv(nvariables);
    Array<OneD, Array<OneD, NekDouble>> inFwd(nvariables);
    Array<OneD, Array<OneD, NekDouble>> inBwd(nvariables);

    for (int i = 0; i < nvariables; ++i)
    {
        inarrayAdv[i] = Array<OneD, NekDouble>(npointsIn, 0.0);
        inFwd[i]      = Array<OneD, NekDouble>(nTracePts, 0.0);
        inBwd[i]      = Array<OneD, NekDouble>(nTracePts, 0.0);
    }

    m_advection->Advect(nvariables, m_advfields, advVel, inarrayAdv,
                        outarrayAdv, time, pFwd, pBwd);

    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Vadd(npointsOut, outarrayAdv[i], 1,
                    outarray[this->advected_fields[i]], 1,
                    outarray[this->advected_fields[i]], 1);
    }
}

/**
 * @brief Compute the advection terms for the right-hand side
 */
void ElectrostaticTurbulence::DoParticles(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    int npts = GetNpoints();

    // Add contribution to electron energy
    Vmath::Vadd(npts, outarray[pe_idx], 1, this->src_fields[0]->GetPhys(), 1,
                outarray[pe_idx], 1);

    for (const auto &[s, v] : this->GetIons())
    {
        int ni_idx  = v.fields.at(field_to_index["n"]);
        int vi_idx  = v.fields.at(field_to_index["v"]);
        int pi_idx  = v.fields.at(field_to_index["e"]);
        double mass = v.mass;
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
        m_varConv->GetIonDynamicEnergy(s, mass, inarray, dynamic_energy);
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
 * @brief Calls HelmSolve to solve for the electric potential
 *
 * @param inarray Array of physical field values
 */
void ElectrostaticTurbulence::SolvePhi(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    [[maybe_unused]] const Array<OneD, NekDouble> &ne)
{
    int npts = GetNpoints();
    StdRegions::ConstFactorMap factors;
    // Helmholtz => Poisson (lambda = 0)
    factors[StdRegions::eFactorLambda] = 0.0;

    Array<OneD, NekDouble> tmp(npts, 0.0);
    Array<OneD, NekDouble> tmp2(npts, 0.0);
    for (const auto &[s, v] : this->GetIons())
    {
        int ni_idx = v.fields.at(field_to_index["n"]);
        int pi_idx = v.fields.at(field_to_index["e"]);
        Vmath::Svtvp(npts, v.mass, inarray[ni_idx], 1, tmp, 1, tmp, 1);
        Vmath::Vdiv(npts, inarray[pi_idx], 1, inarray[ni_idx], 1, tmp2, 1);
        Vmath::Smul(npts, 1.0 / v.charge, tmp2, 1, tmp2, 1);
    }
    Vmath::Vdiv(npts, inarray[omega_idx], 1, tmp, 1, tmp, 1);
    Vmath::Vmul(npts, tmp, 1, this->mag_B, 1, tmp, 1);

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Array<OneD, NekDouble> D(npts, 0.0);
            for (int k = 0; k < npts; k++)
            {
                D[k] = -b_unit[i][k] * b_unit[j][k];
                if (i == j)
                {
                    D[k] += 1;
                }
            }

            m_phi_varcoeff[vc[i][j]] = D;
        }
    }

    // Solve for phi. Output of this routine is in coefficient (spectral)
    // space, so backwards transform to physical space since we'll need that
    // for the advection step & computing drift velocity.

    m_indfields[phi_idx]->HelmSolve(tmp, this->phi->UpdateCoeffs(), factors,
                                    m_phi_varcoeff);

    m_indfields[phi_idx]->BwdTrans(this->phi->GetCoeffs(),
                                   this->phi->UpdatePhys());

    Vmath::Vsub(npts, this->phi->GetPhys(), 1, tmp2, 1, this->phi->UpdatePhys(),
                1);
    m_indfields[phi_idx]->FwdTrans(m_indfields[phi_idx]->GetPhys(),
                                   m_indfields[phi_idx]->UpdateCoeffs());
}

/**
 * @brief Calculates initial potential and gradient
 */
void ElectrostaticTurbulence::CalcInitPhi()
{
    int npts = GetNpoints();

    Array<OneD, Array<OneD, NekDouble>> inarray(m_indfields.size());
    Array<OneD, NekDouble> ne = m_fields[0]->UpdatePhys();
    for (int i = 0; i < m_indfields.size(); i++)
    {
        inarray[i] = m_indfields[i]->GetPhys();
    }
    m_varConv->GetElectronDensity(inarray, ne);
    m_fields[0]->FwdTrans(ne, m_fields[0]->UpdateCoeffs());
    CalcVelocities(inarray);
    CalcInitOmega();

    SolvePhi(inarray, ne);
    ComputeE();
}

/**
 * @brief Compute the gradient of phi for evaluation at the particle positions.
 */
void ElectrostaticTurbulence::ComputeE()
{
    m_indfields[phi_idx]->PhysDeriv(
        this->phi->GetPhys(), this->E[0]->UpdatePhys(),
        this->E[1]->UpdatePhys(), this->E[2]->UpdatePhys());
    int npts = GetNpoints();

    Vmath::Neg(npts, this->E[0]->UpdatePhys(), 1);
    Vmath::Neg(npts, this->E[1]->UpdatePhys(), 1);
    Vmath::Neg(npts, this->E[2]->UpdatePhys(), 1);

    this->E[0]->FwdTrans(this->E[0]->GetPhys(), this->E[0]->UpdateCoeffs());
    this->E[1]->FwdTrans(this->E[1]->GetPhys(), this->E[1]->UpdateCoeffs());
    this->E[2]->FwdTrans(this->E[2]->GetPhys(), this->E[2]->UpdateCoeffs());
}

void ElectrostaticTurbulence::ComputevExB()
{
    int npts = GetNpoints();
    // Calculate ExB velocity

    const Array<OneD, NekDouble> &Ex = this->E[0]->GetPhys();
    const Array<OneD, NekDouble> &Ey = this->E[1]->GetPhys();
    const Array<OneD, NekDouble> &Ez = this->E[2]->GetPhys();
    const Array<OneD, NekDouble> &Bx = this->B[0]->GetPhys();
    const Array<OneD, NekDouble> &By = this->B[1]->GetPhys();
    const Array<OneD, NekDouble> &Bz = this->B[2]->GetPhys();
    if (m_spacedim == 3)
    {
        for (int p = 0; p < npts; ++p)
        {
            this->v_ExB[0][p] =
                (Ey[p] * Bz[p] - Ez[p] * By[p]) / this->mag_B[p];
            this->v_ExB[1][p] =
                (Ez[p] * Bx[p] - Ex[p] * Bz[p]) / this->mag_B[p];
            this->v_ExB[2][p] =
                (Ex[p] * By[p] - Ey[p] * Bx[p]) / this->mag_B[p];
        }
    }
    else
    {
        for (int p = 0; p < npts; ++p)
        {
            this->v_ExB[0][p] = Ey[p] * Bz[p] / this->mag_B[p];
            this->v_ExB[1][p] = -Ex[p] * Bz[p] / this->mag_B[p];
        }
    }
}

void ElectrostaticTurbulence::CalcVelocities(
    const Array<OneD, Array<OneD, NekDouble>> &inarray)
{
    int npts = inarray[0].size();
    for (int f = 0; f < this->adv_vel.size(); ++f)
    {
        for (int d = 0; d < this->adv_vel[f].size(); ++d)
        {
            Vmath::Zero(npts, adv_vel[f][d], 1);
        }
    }
    const Array<OneD, NekDouble> &ne = m_fields[0]->GetPhys();

    // Zero Electron velocity
    Array<OneD, NekDouble> &j_i = m_fields[1]->UpdatePhys();
    Vmath::Zero(npts, j_i, 1);

    for (const auto &[s, v] : this->GetIons())
    {
        int ni_idx = v.fields.at(field_to_index["n"]);
        int vi_idx = v.fields.at(field_to_index["v"]);
        int pi_idx = v.fields.at(field_to_index["e"]);
        // Calculate Ion parallel velocities

        for (int p = 0; p < npts; ++p)
        {
            j_i[p] += v.charge * inarray[vi_idx][p] / v.mass;
            double v_i_par = inarray[vi_idx][p] / (v.mass * inarray[ni_idx][p]);
            for (int d = 0; d < m_spacedim; ++d)
            {
                this->adv_vel[ni_idx][d][p] =
                    this->v_ExB[d][p] + v_i_par * this->b_unit[d][p];
                this->adv_vel[vi_idx][d][p] =
                    this->v_ExB[d][p] + v_i_par * this->b_unit[d][p];
                this->adv_vel[pi_idx][d][p] =
                    this->v_ExB[d][p] + v_i_par * this->b_unit[d][p];
            }
        }
    }
    Array<OneD, NekDouble> j_par(npts, 0.0);
    // TODO calculate conductivity
    double sigma = 0;
    for (int p = 0; p < npts; ++p)
    {
        for (int d = 0; d < m_spacedim; ++d)
        {
            j_par[p] += sigma * this->E[d]->GetPhys()[p] * this->b_unit[d][p];
        }
    }
    for (int p = 0; p < npts; ++p)
    {
        for (int d = 0; d < m_spacedim; ++d)
        {
            this->adv_vel[pe_idx][d][p] =
                this->v_ExB[d][p] +
                this->b_unit[d][p] * (j_i[p] - j_par[p]) / ne[p];
        }
    }

    for (const auto &[s, v] : GetNeutrals())
    {
        int nn_idx = v.fields.at(field_to_index["n"]);
        int vn_idx = v.fields.at(field_to_index["v"]);
        int pn_idx = v.fields.at(field_to_index["e"]);
        // Calculate Neutral parallel velocities
        for (int p = 0; p < npts; ++p)
        {
            double v_i_par = inarray[vn_idx][p] / (v.mass * inarray[nn_idx][p]);
            for (int d = 0; d < m_spacedim; ++d)
            {
                this->adv_vel[nn_idx][d][p] = v_i_par * this->b_unit[d][p];
                this->adv_vel[vn_idx][d][p] = v_i_par * this->b_unit[d][p];
                this->adv_vel[pn_idx][d][p] = v_i_par * this->b_unit[d][p];
            }
        }
    }
}

void ElectrostaticTurbulence::AddDriftVelocities(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    const Array<OneD, NekDouble> &ne,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &adv_vel)
{
    int npts = inarray[0].size();

    for (const auto &[s, v] : this->GetIons())
    {
        int ni_idx = v.fields.at(field_to_index["n"]);
        int pi_idx = v.fields.at(field_to_index["e"]);
        Array<OneD, NekDouble> gradp[3];

        for (int d = 0; d < 3; ++d)
        {

            gradp[d] = Array<OneD, NekDouble>(npts, 0.0);
        }

        if (m_spacedim == 2)
            m_indfields[pi_idx]->PhysDeriv(inarray[pi_idx], gradp[0], gradp[1]);
        else if (m_spacedim == 3)
            m_indfields[pi_idx]->PhysDeriv(inarray[pi_idx], gradp[0], gradp[1],
                                           gradp[2]);

        Vmath::Vvtvvtm(npts, gradp[1], 1, this->B[2]->GetPhys(), 1, gradp[2], 1,
                       this->B[1]->GetPhys(), 1, this->v_di[s][0], 1);
        Vmath::Vdiv(npts, this->v_di[s][0], 1, this->mag_B, 1, this->v_di[s][0],
                    1);
        Vmath::Vdiv(npts, this->v_di[s][0], 1, inarray[ni_idx], 1,
                    this->v_di[s][0], 1);
        Vmath::Smul(npts, -2.0 / (3.0 * v.charge), this->v_di[s][0], 1,
                    this->v_di[s][0], 1);
        Vmath::Vvtvvtm(npts, gradp[2], 1, this->B[0]->GetPhys(), 1, gradp[0], 1,
                       this->B[2]->GetPhys(), 1, this->v_di[s][1], 1);
        Vmath::Vdiv(npts, this->v_di[s][1], 1, this->mag_B, 1, this->v_di[s][1],
                    1);
        Vmath::Vdiv(npts, this->v_di[s][1], 1, inarray[ni_idx], 1,
                    this->v_di[s][1], 1);
        Vmath::Smul(npts, -2.0 / (3.0 * v.charge), this->v_di[s][1], 1,
                    this->v_di[s][1], 1);

        if (m_spacedim == 3)
        {
            Vmath::Vvtvvtm(npts, gradp[0], 1, this->B[1]->GetPhys(), 1,
                           gradp[1], 1, this->B[0]->GetPhys(), 1,
                           this->v_di[s][2], 1);
            Vmath::Vdiv(npts, this->v_di[s][2], 1, this->mag_B, 1,
                        this->v_di[s][2], 1);
            Vmath::Vdiv(npts, this->v_di[s][2], 1, inarray[ni_idx], 1,
                        this->v_di[s][2], 1);
            Vmath::Smul(npts, -2.0 / (3.0 * v.charge), this->v_di[s][2], 1,
                        this->v_di[s][2], 1);
        }
        for (int d = 0; d < m_spacedim; ++d)
        {
            Vmath::Vadd(npts, this->v_di[s][d], 1, this->adv_vel[ni_idx][d], 1,
                        this->adv_vel[ni_idx][d], 1);
            Vmath::Svtvp(npts, 5.0 / 3.0, this->v_di[s][d], 1,
                         adv_vel[pi_idx][d], 1, adv_vel[pi_idx][d], 1);
        }
    }
    // Calculate electron diagmagnetic velocity

    Array<OneD, NekDouble> gradp[3];

    for (int d = 0; d < 3; ++d)
    {

        gradp[d] = Array<OneD, NekDouble>(npts, 0.0);
    }
    if (m_spacedim == 2)
    {
        m_indfields[pe_idx]->PhysDeriv(inarray[pe_idx], gradp[0], gradp[1],
                                       NullNekDouble1DArray);
    }

    else if (m_spacedim == 3)
    {
        m_indfields[pe_idx]->PhysDeriv(inarray[pe_idx], gradp[0], gradp[1],
                                       gradp[2]);
    }

    Vmath::Vvtvvtm(npts, gradp[1], 1, this->B[2]->GetPhys(), 1, gradp[2], 1,
                   this->B[1]->GetPhys(), 1, this->v_de[0], 1);
    Vmath::Vdiv(npts, this->v_de[0], 1, this->mag_B, 1, this->v_de[0], 1);
    Vmath::Vdiv(npts, this->v_de[0], 1, ne, 1, this->v_de[0], 1);
    Vmath::Smul(npts, 2.0 / 3.0, this->v_de[0], 1, this->v_de[0], 1);

    Vmath::Vvtvvtm(npts, gradp[2], 1, this->B[0]->GetPhys(), 1, gradp[0], 1,
                   this->B[2]->GetPhys(), 1, this->v_de[1], 1);
    Vmath::Vdiv(npts, this->v_de[1], 1, this->mag_B, 1, this->v_de[1], 1);
    Vmath::Vdiv(npts, this->v_de[1], 1, ne, 1, this->v_de[1], 1);
    Vmath::Smul(npts, 2.0 / 3.0, this->v_de[1], 1, this->v_de[1], 1);

    if (m_spacedim == 3)
    {
        Vmath::Vvtvvtm(npts, gradp[0], 1, this->B[1]->GetPhys(), 1, gradp[1], 1,
                       this->B[0]->GetPhys(), 1, this->v_de[2], 1);
        Vmath::Vdiv(npts, this->v_de[2], 1, this->mag_B, 1, this->v_de[2], 1);
        Vmath::Vdiv(npts, this->v_de[2], 1, ne, 1, this->v_de[2], 1);
        Vmath::Smul(npts, 2.0 / 3.0, this->v_de[2], 1, this->v_de[2], 1);
    }

    for (int d = 0; d < m_spacedim; ++d)
    {
        Vmath::Svtvp(npts, 5.0 / 3.0, this->v_de[d], 1, adv_vel[pe_idx][d], 1,
                     adv_vel[pe_idx][d], 1);
    }
}

void ElectrostaticTurbulence::CalcInitOmega()
{
    int npts = GetNpoints();

    const Array<OneD, NekDouble> &Bx = this->B[0]->GetPhys();
    const Array<OneD, NekDouble> &By = this->B[1]->GetPhys();
    const Array<OneD, NekDouble> &Bz = this->B[2]->GetPhys();

    Array<OneD, Array<OneD, NekDouble>> w(m_spacedim);
    for (int d = 0; d < m_spacedim; ++d)
    {
        w[d] = Array<OneD, NekDouble>(npts, 0.0);
    }

    for (const auto &[s, v] : this->GetIons())
    {
        int ni_idx                       = v.fields.at(field_to_index["n"]);
        const Array<OneD, NekDouble> &ni = m_indfields[ni_idx]->GetPhys();
        // Calculate sum of w = (n_i) v0 X b/|B|
        if (m_spacedim == 3)
        {
            for (int p = 0; p < npts; ++p)
            {
                w[0][p] += v.mass * ni[p] *
                           (this->adv_vel[ni_idx][1][p] * Bz[p] -
                            this->adv_vel[ni_idx][2][p] * By[p]) /
                           this->mag_B[p];
                w[1][p] += v.mass * ni[p] *
                           (this->adv_vel[ni_idx][2][p] * Bx[p] -
                            this->adv_vel[ni_idx][0][p] * Bz[p]) /
                           this->mag_B[p];
                w[2][p] += v.mass * ni[p] *
                           (this->adv_vel[ni_idx][0][p] * By[p] -
                            this->adv_vel[ni_idx][1][p] * Bx[p]) /
                           this->mag_B[p];
            }
        }
        else if (m_spacedim == 2)
        {
            for (int p = 0; p < npts; ++p)
            {
                w[0][p] += v.mass * ni[p] * this->adv_vel[ni_idx][1][p] *
                           Bz[p] / this->mag_B[p];
                w[1][p] += -v.mass * ni[p] * this->adv_vel[ni_idx][0][p] *
                           Bz[p] / this->mag_B[p];
            }
        }
    }
    for (int d = 0; d < m_spacedim; ++d)
    {
        m_indfields[omega_idx]->PhysDeriv(d, w[d], w[d]);

        Vmath::Vadd(npts, w[d], 1, m_indfields[omega_idx]->GetPhys(), 1,
                    m_indfields[omega_idx]->UpdatePhys(), 1);
    }
    Vmath::Smul(npts, 1.0 / this->omega_c, m_indfields[omega_idx]->GetPhys(), 1,
                m_indfields[omega_idx]->UpdatePhys(), 1);

    m_indfields[omega_idx]->FwdTransLocalElmt(
        m_indfields[omega_idx]->GetPhys(),
        m_indfields[omega_idx]->UpdateCoeffs());
}

void ElectrostaticTurbulence::CalcOmegaFlux(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &omega_flux)
{
    int npts                         = inarray[0].size();
    const Array<OneD, NekDouble> &Bx = this->B[0]->GetPhys();
    const Array<OneD, NekDouble> &By = this->B[1]->GetPhys();
    const Array<OneD, NekDouble> &Bz = this->B[2]->GetPhys();
    for (int d = 0; d < m_spacedim; ++d)
    {
        Vmath::Zero(npts, omega_flux[d], 1);
    }
    for (const auto &[s, v] : this->GetIons())
    {
        int ni_idx = v.fields.at(field_to_index["n"]);

        // Calculate w = (n_i) v0 X b/|B|
        Array<OneD, Array<OneD, NekDouble>> w(m_spacedim);
        if (m_spacedim == 3)
        {
            for (int p = 0; p < npts; ++p)
            {
                w[0][p] = inarray[ni_idx][p] *
                          (this->adv_vel[ni_idx][1][p] * Bz[p] -
                           this->adv_vel[ni_idx][2][p] * By[p]) /
                          this->mag_B[p];
                w[1][p] = inarray[ni_idx][p] *
                          (this->adv_vel[ni_idx][2][p] * Bx[p] -
                           this->adv_vel[ni_idx][0][p] * Bz[p]) /
                          this->mag_B[p];
                w[2][p] = inarray[ni_idx][p] *
                          (this->adv_vel[ni_idx][0][p] * By[p] -
                           this->adv_vel[ni_idx][1][p] * Bx[p]) /
                          this->mag_B[p];
            }
        }
        else if (m_spacedim == 2)
        {
            for (int p = 0; p < npts; ++p)
            {
                w[0][p] = inarray[ni_idx][p] * this->adv_vel[ni_idx][1][p] *
                          Bz[p] / this->mag_B[p];
                w[1][p] = -inarray[ni_idx][p] * this->adv_vel[ni_idx][0][p] *
                          Bz[p] / this->mag_B[p];
            }
        }

        for (int d = 0; d < m_spacedim; ++d)
        {
            Array<OneD, NekDouble> tmp(npts, 0.0);

            for (int d2 = 0; d2 < m_spacedim; ++d2)
            {
                Array<OneD, NekDouble> tmp2(npts, 0.0);

                // Calculate ∇⋅(w⊗v0)
                Vmath::Vmul(npts, this->adv_vel[ni_idx][d], 1, w[d2], 1, tmp2,
                            1);
                m_indfields[omega_idx]->PhysDeriv(d2, tmp2, tmp2);
                Vmath::Vadd(npts, tmp2, 1, tmp, 1, tmp, 1);
            }
            // Vorticity Flux
            Vmath::Smul(npts, v.mass / this->omega_c, tmp, 1, tmp, 1);
            Vmath::Vadd(npts, tmp, 1, omega_flux[d], 1, omega_flux[d], 1);
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
Array<OneD, Array<OneD, NekDouble>> &ElectrostaticTurbulence::GetAdvVelNorm()
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

Array<OneD, NekDouble> &ElectrostaticTurbulence::GetOmegaFlux()
{
    int num_trace_pts = GetTraceNpoints();
    Array<OneD, Array<OneD, NekDouble>> normals(m_spacedim);

    for (int d = 0; d < m_spacedim; ++d)
    {
        normals[d] = Array<OneD, NekDouble>(num_trace_pts);
    }

    m_indfields[omega_idx]->GetTrace()->GetNormals(normals);
    Vmath::Zero(num_trace_pts, this->omega_flux_norm, 1);
    for (int d = 0; d < m_spacedim; ++d)
    {
        m_indfields[omega_idx]->ExtractTracePhys(this->omega_flux[d],
                                                 this->omega_flux_trace[d]);
        for (int p = 0; p < num_trace_pts; ++p)
        {
            this->omega_flux_norm[p] +=
                normals[d][p] * this->omega_flux_trace[d][p];
        }
    }

    return this->omega_flux_norm;
}

/**
 *  @brief Construct flux array.
 *
 * @param  field_vals Physical values for each advection field
 * @param[out] flux       Flux array
 */
void ElectrostaticTurbulence::GetFluxVector(
    const Array<OneD, Array<OneD, NekDouble>> &field_vals,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fluxes)
{
    int npts = field_vals[0].size();

    for (int d = 0; d < m_spacedim; ++d)
    {
        for (int p = 0; p < npts; ++p)
        {
            fluxes[pe_idx][d][p] =
                this->adv_vel[pe_idx][d][p] * field_vals[pe_idx][p];
        }
    }

    // Omega flux
    for (int d = 0; d < m_spacedim; ++d)
    {
        // Vmath::Vcopy(npts, this->omega_flux[d], 1, fluxes[omega_idx][d], 1);
        Vmath::Vmul(npts, field_vals[omega_idx], 1, this->adv_vel[omega_idx][d],
                    1, fluxes[omega_idx][d], 1);
    }

    for (const auto &[s, v] : this->GetSpecies())
    {
        int ni_idx = v.fields.at(field_to_index["n"]);
        int vi_idx = v.fields.at(field_to_index["v"]);
        int pi_idx = v.fields.at(field_to_index["p"]);

        for (int d = 0; d < m_spacedim; ++d)
        {
            for (int p = 0; p < npts; ++p)
            {
                fluxes[ni_idx][d][p] =
                    this->adv_vel[ni_idx][d][p] * field_vals[ni_idx][p];
                fluxes[vi_idx][d][p] =
                    this->adv_vel[vi_idx][d][p] * field_vals[vi_idx][p];
                fluxes[pi_idx][d][p] =
                    this->adv_vel[pi_idx][d][p] * field_vals[pi_idx][p];
            }
        }
    }
}

void ElectrostaticTurbulence::DoDiffusion(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray,
    const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    const Array<OneD, Array<OneD, NekDouble>> &pBwd)
{
    int nvariables = inarray.size() - 1;
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

    for (const auto &[s, v] : this->GetSpecies())
    {
        int pi_idx = v.fields.at(field_to_index["e"]);
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
        for (const auto &[s, v] : this->GetSpecies())
        {
            int pi_idx = v.fields.at(field_to_index["e"]);
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

void ElectrostaticTurbulence::CalcKPar()
{
    // Change to fn of fields
    int npoints = m_fields[0]->GetNpoints();
    NekDouble k_par;
    m_session->LoadParameter("k_par", k_par, 100.0);
    m_kpar = Array<OneD, NekDouble>(npoints, k_par);
}

void ElectrostaticTurbulence::CalcKPerp()
{
    // Change to fn of fields
    int npoints = m_fields[0]->GetNpoints();
    NekDouble k_perp;
    m_session->LoadParameter("k_perp", k_perp, 1.0);
    m_kperp = Array<OneD, NekDouble>(npoints, k_perp);
}

void ElectrostaticTurbulence::CalcDiffTensor()
{
    int npoints = m_fields[0]->GetNpoints();
    CalcKPar();
    CalcKPerp();
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

void ElectrostaticTurbulence::CalcKappaPar()
{
    // Change to fn of T
    int npoints = m_fields[0]->GetNpoints();
    NekDouble kappa_par;
    m_session->LoadParameter("kappa_par", kappa_par, 0.1);
    m_kappapar = Array<OneD, NekDouble>(npoints, kappa_par);
}

void ElectrostaticTurbulence::CalcKappaPerp()
{
    // Change to fn of T
    int npoints = m_fields[0]->GetNpoints();
    NekDouble kappa_perp;
    m_session->LoadParameter("kappa_perp", kappa_perp, 0.1);
    m_kappaperp = Array<OneD, NekDouble>(npoints, kappa_perp);
}

void ElectrostaticTurbulence::CalcKappaTensor()
{
    int npoints = m_fields[0]->GetNpoints();
    CalcKappaPar();
    CalcKappaPerp();
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Array<OneD, NekDouble> kappa(npoints, 0.0);
            for (int k = 0; k < npoints; k++)
            {
                kappa[k] = (m_kappapar[k] - m_kappaperp[k]) * b_unit[i][k] *
                           b_unit[j][k];
                if (i == j)
                {
                    kappa[k] += m_kappaperp[k];
                }
            }
            m_kappa[vc[i][j]] = kappa;
        }
    }
}

/**
 * @brief Construct the flux vector for the anisotropic diffusion problem.
 */
void ElectrostaticTurbulence::GetFluxVectorDiff(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &qfield,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fluxes)
{
    unsigned int nPts = inarray[0].size();

    for (int j = 0; j < m_spacedim; ++j)
    {
        // Calc diffusion of n with D tensor
        Vmath::Vmul(nPts, m_kappa[vc[j][0]].GetValue(), 1, qfield[0][pe_idx], 1,
                    fluxes[j][pe_idx], 1);
        Vmath::Vmul(nPts, m_kappa[vc[j][0]].GetValue(), 1, qfield[0][omega_idx],
                    1, fluxes[j][omega_idx], 1);
        for (int k = 1; k < m_spacedim; ++k)
        {
            Vmath::Vvtvp(nPts, m_kappa[vc[j][k]].GetValue(), 1,
                         qfield[k][pe_idx], 1, fluxes[j][pe_idx], 1,
                         fluxes[j][pe_idx], 1);
            Vmath::Vvtvp(nPts, m_kappa[vc[j][0]].GetValue(), 1,
                         qfield[k][omega_idx], 1, fluxes[j][omega_idx], 1,
                         fluxes[j][omega_idx], 1);
        }
    }
    for (const auto &[s, v] : this->GetSpecies())
    {
        int pi_idx = v.fields.at(field_to_index["e"]);
        for (int j = 0; j < m_spacedim; ++j)
        {
            Vmath::Vmul(nPts, m_kappa[vc[j][0]].GetValue(), 1,
                        qfield[0][pi_idx], 1, fluxes[j][pi_idx], 1);
            for (int k = 1; k < m_spacedim; ++k)
            {
                Vmath::Vvtvp(nPts, m_kappa[vc[j][k]].GetValue(), 1,
                             qfield[k][pi_idx], 1, fluxes[j][pi_idx], 1,
                             fluxes[j][pi_idx], 1);
            }
        }
    }
}

void ElectrostaticTurbulence::AddNeutralSources(
    const Array<OneD, Array<OneD, NekDouble>> &in_arr,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
}

/**
 * @brief Populate rhs array ( @p outarray )
 *
 * @param inarray physical values of all fields
 * @param[out] outarray output array (RHSs of time integration equations)
 */
void ElectrostaticTurbulence::DoOdeImplicitRhs(
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
void ElectrostaticTurbulence::DoOdeRhsCoeff(
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
void ElectrostaticTurbulence::DoAdvectionCoeff(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
    const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    const Array<OneD, Array<OneD, NekDouble>> &pBwd)
{
    int nvariables = inarray.size() - 1;

    std::dynamic_pointer_cast<SU::AdvectionWeakDG>(m_advection)
        ->AdvectCoeffs(nvariables, m_indfields, this->v_ExB, inarray, outarray,
                       time, pFwd, pBwd);
}

void ElectrostaticTurbulence::DoDiffusionCoeff(
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
    Array<OneD, Array<OneD, NekDouble>> inarrayDiff{nvariables - 1};
    Array<OneD, Array<OneD, NekDouble>> inFwd{nvariables - 1};
    Array<OneD, Array<OneD, NekDouble>> inBwd{nvariables - 1};

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
        int pi_idx = v.fields.at(field_to_index["e"]);
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
            int pi_idx = v.fields.at(field_to_index["e"]);
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
void ElectrostaticTurbulence::DoParticlesCoeff(
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
        int pi_idx = v.fields.at(field_to_index["e"]);
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
        m_fields[pi_idx]->FwdTrans(dynamic_energy, tmp);
        Vmath::Vadd(npts, outarray[pi_idx], 1, tmp, 1, outarray[pi_idx], 1);

        Vmath::Zero(npts, dynamic_energy, 1);
        for (int d = 0; d < m_spacedim; ++d)
        {
            Vmath::Vvtvp(npts, this->b_unit[d], 1,
                         this->src_fields[vi_src_idx[s] + d]->GetPhys(), 1,
                         dynamic_energy, 1, dynamic_energy, 1);
        }
        m_fields[vi_idx]->FwdTrans(dynamic_energy, tmp);
        Vmath::Vadd(npts, outarray[vi_idx], 1, tmp, 1, outarray[vi_idx], 1);
    }
}

/**
 * @brief After reading ICs, calculate phi and grad(phi)
 */
void ElectrostaticTurbulence::v_SetInitialConditions(NekDouble init_time,
                                                     bool dump_ICs,
                                                     const int domain)
{
    PlasmaSystem::v_SetInitialConditions(init_time, dump_ICs, domain);
    CalcInitPhi();
    Checkpoint_Output(0);
}

// void ElectrostaticTurbulence::v_SetBoundaryConditions(
//     Array<OneD, Array<OneD, NekDouble>> &physarray, NekDouble time)
// {
//     size_t nTracePts  = GetTraceTotPoints();
//     size_t nvariables = physarray.size() - 1;
//     Array<OneD, Array<OneD, NekDouble>> Fwd(nvariables);
//     for (size_t i = 0; i < nvariables; ++i)
//     {
//         Fwd[i] = Array<OneD, NekDouble>(nTracePts);
//         int p  = i < n_species * n_fields_per_species
//                      ? i % n_fields_per_species
//                      : i - n_species * n_fields_per_species;
//         m_fields[p]->ExtractTracePhys(physarray[i], Fwd[i]);
//     }
//     for (auto &bc : m_bndConds)
//     {
//         bc->Apply(Fwd, physarray, time);
//     }
// }

void ElectrostaticTurbulence::load_params()
{
    PlasmaSystem::load_params();
    // Type of advection to use. Default is DG.
    m_session->LoadSolverInfo("AdvectionType", this->adv_type, "WeakDG");
    // Type of Riemann solver to use. Default = "Upwind"
    m_session->LoadSolverInfo("UpwindType", this->riemann_solver_type,
                              "VorticityUpwind");
    std::string boussinesq_str;
    m_session->LoadSolverInfo("Boussinesq Approximation", boussinesq_str, "On");
    this->m_boussinesq = (boussinesq_str == "On");
}

void ElectrostaticTurbulence::v_ExtraFldOutput(
    std::vector<Array<OneD, NekDouble>> &fieldcoeffs,
    std::vector<std::string> &variables)
{
    PlasmaSystem::v_ExtraFldOutput(fieldcoeffs, variables);
    const int nPhys   = m_fields[0]->GetNpoints();
    const int nCoeffs = m_fields[0]->GetNcoeffs();

    m_fields[0]->FwdTransLocalElmt(this->phi->GetPhys(), fieldcoeffs[4]);

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