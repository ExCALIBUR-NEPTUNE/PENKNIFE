#include "DoubleDiffusiveField.hpp"

namespace PENKNIFE
{

/// Name of class
static std::string class_name;
std::string DoubleDiffusiveField::class_name =
    SU::GetEquationSystemFactory().RegisterCreatorFunction(
        "DoubleDiffusiveField", DoubleDiffusiveField::create,
        "Solves for two diffusive fields (n, e) with anisotropy");
/**
 * @brief Creates an instance of this class.
 */
static SU::EquationSystemSharedPtr create(
    const LU::SessionReaderSharedPtr &session,
    const SD::MeshGraphSharedPtr &graph)
{
    SU::EquationSystemSharedPtr p =
        MemoryManager<DoubleDiffusiveField>::AllocateSharedPtr(session, graph);
    p->InitObject();
    return p;
}

DoubleDiffusiveField::DoubleDiffusiveField(
    const LU::SessionReaderSharedPtr &session,
    const SD::MeshGraphSharedPtr &graph)
    : PlasmaSystem(session, graph)
{
    this->n_indep_fields = 1;
}

void DoubleDiffusiveField::v_InitObject(bool DeclareFields)
{
    PlasmaSystem::v_InitObject(DeclareFields);

    std::string diffName;
    m_session->LoadSolverInfo("DiffusionType", diffName, "LDG");
    m_diffusion =
        SolverUtils::GetDiffusionFactory().CreateInstance(diffName, diffName);
    m_diffusion->SetFluxVector(&DoubleDiffusiveField::GetFluxVectorDiff, this);

    // workaround for bug in DiffusionLDG
    m_difffields = Array<OneD, MR::ExpListSharedPtr>(m_indfields.size());
    for (int f = 0; f < m_difffields.size(); ++f)
    {
        m_difffields[f] = m_indfields[f];
    }

    m_diffusion->InitObject(m_session, m_difffields);
    int npts       = GetNpoints();
    this->m_kpar   = Array<OneD, NekDouble>(npts, 0.0);
    this->m_kcross = Array<OneD, NekDouble>(npts, 0.0);
    this->m_kperp  = Array<OneD, NekDouble>(npts, 0.0);

    // this->ne = std::dynamic_pointer_cast<MR::DisContField>(m_fields[0]);
    // this->Te = MemoryManager<MR::DisContField>::AllocateSharedPtr(*ne);

    pe_idx = m_indfields.size() - this->n_indep_fields;

    m_ode.DefineOdeRhs(&DoubleDiffusiveField::DoOdeRhs, this);

    if (this->particles_enabled)
    {
        std::vector<Sym<REAL>> src_syms;
        std::vector<int> src_components;
        this->src_fields.emplace_back(
            MemoryManager<MR::DisContField>::AllocateSharedPtr(
                *std::dynamic_pointer_cast<MR::DisContField>(m_fields[0])));
        src_syms.push_back(Sym<REAL>("ELECTRON_SOURCE_DENSITY"));
        src_components.push_back(0);

        this->src_fields.emplace_back(
            MemoryManager<MR::DisContField>::AllocateSharedPtr(
                *std::dynamic_pointer_cast<MR::DisContField>(m_fields[0])));
        src_syms.push_back(Sym<REAL>("ELECTRON_SOURCE_ENERGYY"));
        src_components.push_back(0);

        for (auto &[k, v] : this->particle_sys->get_species())
        {
            this->src_fields.emplace_back(
                MemoryManager<MR::DisContField>::AllocateSharedPtr(
                    *std::dynamic_pointer_cast<MR::DisContField>(m_fields[0])));
            this->src_fields.emplace_back(
                MemoryManager<MR::DisContField>::AllocateSharedPtr(
                    *std::dynamic_pointer_cast<MR::DisContField>(m_fields[0])));
            src_syms.push_back(Sym<REAL>(k + "_SOURCE_DENSITY"));
            src_components.push_back(0);
            src_syms.push_back(Sym<REAL>(k + "_SOURCE_ENERGY"));
            src_components.push_back(0);
        }

        this->src_fields.emplace_back(
            MemoryManager<MR::DisContField>::AllocateSharedPtr(
                *std::dynamic_pointer_cast<MR::DisContField>(m_fields[0])));
        src_syms.push_back(Sym<REAL>("ELECTRON_SOURCE_ENERGY"));
        src_components.push_back(0);

        this->particle_sys->setup_evaluate_fields(this->E, this->B, this->ne,
                                                  this->Te, this->ve);

        this->particle_sys->finish_setup(this->src_fields, src_syms,
                                         src_components);
    }
}

void DoubleDiffusiveField::ImplicitTimeIntCG(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray,
    [[maybe_unused]] const NekDouble time, const NekDouble lambda)
{

    int nvariables                       = inarray.size();
    int npoints                          = m_indfields[0]->GetNpoints();
    m_factors[StdRegions::eFactorLambda] = 1.0 / lambda / m_epsilon;

    // We solve ( \nabla^2 - HHlambda ) Y[i] = rhs [i]
    // inarray = input: \hat{rhs} -> output: \hat{Y}
    // outarray = output: nabla^2 \hat{Y}
    // where \hat = modal coeffs
    if (m_intScheme->GetIntegrationSchemeType() == LibUtilities::eImplicit)
    {
        for (int i = 0; i < nvariables; ++i)
        {
            Vmath::Zero(npoints, outarray[i], 1);
        }

        if (this->particles_enabled)
        {
            for (int i = 0; i < this->src_fields.size(); ++i)
            {
                Vmath::Vadd(outarray[i].size(), outarray[i], 1,
                            this->src_fields[i]->GetPhys(), 1, outarray[i], 1);
            }
        }

        for (auto &x : m_forcing)
        {
            x->Apply(m_indfields, inarray, outarray, time);
        }
    }
    CalcDiffTensor();
    Vmath::Zero(npoints, m_fields[0]->UpdatePhys(), 1);

    for (int i = 0; i < nvariables; ++i)
    {
        if (m_intScheme->GetIntegrationSchemeType() == LibUtilities::eImplicit)
        {
            // Multiply forcing term by -1 for definition of HelmSolve function
            Vmath::Smul(npoints, -1.0, outarray[i], 1, outarray[i], 1);

            // Multiply 1.0/timestep/lambda
            Vmath::Svtvp(npoints, -m_factors[StdRegions::eFactorLambda],
                         inarray[i], 1, outarray[i], 1, outarray[i], 1);
        }
        else
        {
            // Multiply 1.0/timestep/lambda
            Vmath::Smul(npoints, -m_factors[StdRegions::eFactorLambda],
                        inarray[i], 1, outarray[i], 1);
        }

        if (i % 2 == 1)
        {
            Vmath::Vadd(npoints, inarray[i - 1], 1, m_fields[0]->UpdatePhys(),
                        1, m_fields[0]->UpdatePhys(), 1);
            Vmath::Vdiv(npoints, outarray[i], 1, inarray[i - 1], 1, outarray[i],
                        1);
        }
        else if (i == nvariables - 1)
        {
            Vmath::Vdiv(npoints, outarray[i], 1, m_fields[0]->GetPhys(), 1,
                        outarray[i], 1);
        }
        m_indfields[i]->HelmSolve(outarray[i], m_indfields[i]->UpdateCoeffs(),
                                  m_factors, m_D);

        m_indfields[i]->BwdTrans(m_indfields[i]->GetCoeffs(), outarray[i]);
    }

    Vmath::Zero(npoints, m_fields[0]->UpdatePhys(), 1);
    for (int i = 0; i < nvariables; ++i)
    {
        if (i % 2 == 1)
        {
            Vmath::Vadd(npoints, outarray[i - 1], 1, m_fields[0]->UpdatePhys(),
                        1, m_fields[0]->UpdatePhys(), 1);
            Vmath::Vmul(npoints, outarray[i], 1, outarray[i - 1], 1,
                        outarray[i], 1);
        }
        else if (i == nvariables - 1)
        {
            Vmath::Vmul(npoints, outarray[i], 1, m_fields[0]->GetPhys(), 1,
                        outarray[i], 1);
        }
        m_indfields[i]->SetPhysState(false);
    }
}

void DoubleDiffusiveField::CalcK(
    const Array<OneD, Array<OneD, NekDouble>> &in_arr, int f)
{
    int npoints = m_fields[0]->GetNpoints();
    double Z, A;
    // this->neso_config->load_species_parameter(f, "Charge", Z);
    // this->neso_config->load_species_parameter(f, "Mass", A);
    auto ne    = this->m_fields[0]->GetPhys();
    int ni_idx = m_ions[f].fields.at(field_to_index.at("n"));

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
void DoubleDiffusiveField::CalcKappa(
    const Array<OneD, Array<OneD, NekDouble>> &in_arr, int f)
{
    int npoints = m_fields[0]->GetNpoints();
    double Z    = this->m_ions[f].charge;
    double A    = this->m_ions[f].mass;
    int ei_idx  = this->m_ions[f].fields[field_to_index["e"]];
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
                          (in_arr[ei_idx][p], 2.5) / (sqrt(A) * Z * Z * tmp[p]);
        this->m_kperp[p] = this->kappa_i_perp * sqrt(A) * tmp[p] *
                           in_arr[ni_idx][p] /
                           (this->mag_B[p] * sqrt(in_arr[ei_idx][p]));
        this->m_kcross[p] = this->kappa_i_cross * in_arr[ni_idx][p] *
                            in_arr[ei_idx][p] / (Z * sqrt(this->mag_B[p]));
    }
}

// Electron thermal conductivity
void DoubleDiffusiveField::CalcKappa(
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

void DoubleDiffusiveField::CalcDiffTensor()
{
    int npoints    = m_fields[0]->GetNpoints();
    int nvariables = m_indfields.size();

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

/**
 * @brief Populate rhs array ( @p out_arr ) with diffused quantitites
 * @param in_arr physical values of all fields
 * @param[out] out_arr output array (RHSs of time integration equations)
 */
void DoubleDiffusiveField::DoOdeRhs(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
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

    DoDiffusion(inarray, outarray, Fwd, Bwd);

    if (this->particles_enabled)
    {
        for (int i = 0; i < this->src_fields.size(); ++i)
        {
            Vmath::Vadd(outarray[i].size(), outarray[i], 1,
                        this->src_fields[i]->GetPhys(), 1, outarray[i], 1);
        }
    }

    // Add forcing terms
    for (auto &x : m_forcing)
    {
        x->Apply(m_fields, inarray, outarray, time);
    }
}

void DoubleDiffusiveField::DoDiffusion(
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

    for (const auto &[s, v] : this->GetSpecies())
    {
        int ni_idx = v.fields.at(field_to_index["n"]);
        int ei_idx = v.fields.at(field_to_index["e"]);
        Vmath::Vcopy(npointsIn, inarray[ni_idx], 1, inarrayDiff[ni_idx], 1);
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
        m_varConv->GetElectronTemperature(pFwd, inFwd[pe_idx]);
        m_varConv->GetElectronTemperature(pBwd, inBwd[pe_idx]);

        for (const auto &[s, v] : this->GetSpecies())
        {
            int ni_idx = v.fields.at(field_to_index["n"]);
            int ei_idx = v.fields.at(field_to_index["e"]);
            Vmath::Vcopy(nTracePts, pFwd[ni_idx], 1, inFwd[ni_idx], 1);
            Vmath::Vcopy(nTracePts, pBwd[ni_idx], 1, inBwd[ni_idx], 1);

            m_varConv->GetIonTemperature(s, v.mass, pFwd, inFwd[ei_idx]);
            m_varConv->GetIonTemperature(s, v.mass, pBwd, inBwd[ei_idx]);
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
/**
 * @brief Construct the flux vector for the anisotropic diffusion problem.
 */
void DoubleDiffusiveField::GetFluxVectorDiff(
    const Array<OneD, Array<OneD, NekDouble>> &in_arr,
    const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &qfield,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fluxes)
{
    unsigned int nDim = qfield.size();
    unsigned int nFld = qfield[0].size();
    unsigned int nPts = qfield[0][0].size();

    for (const auto &[s, v] : this->GetSpecies())
    {
        CalcK(in_arr, s);
        CalcDiffTensor();
        int ni_idx = v.fields.at(field_to_index["n"]);
        int ei_idx = v.fields.at(field_to_index["e"]);

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
            Vmath::Vvtvvtm(nPts, b_unit[1], 1, qfield[2][ei_idx], 1, b_unit[2],
                           1, qfield[1][ei_idx], 1, fluxes[0][ei_idx], 1);
            Vmath::Vvtvvtm(nPts, b_unit[2], 1, qfield[0][ei_idx], 1, b_unit[0],
                           1, qfield[2][ei_idx], 1, fluxes[1][ei_idx], 1);
            Vmath::Vvtvvtm(nPts, b_unit[0], 1, qfield[1][ei_idx], 1, b_unit[1],
                           1, qfield[0][ei_idx], 1, fluxes[2][ei_idx], 1);
        }
        else
        {
            Vmath::Vmul(nPts, b_unit[2], 1, qfield[1][ei_idx], 1,
                        fluxes[0][ei_idx], 1);
            Vmath::Neg(nPts, fluxes[0][ei_idx], 1);
            Vmath::Vmul(nPts, b_unit[2], 1, qfield[0][ei_idx], 1,
                        fluxes[1][ei_idx], 1);
        }

        CalcKappa(in_arr, s);
        CalcDiffTensor();
        for (unsigned int j = 0; j < nDim; ++j)
        {
            Vmath::Vmul(nPts, m_kcross, 1, fluxes[j][ei_idx], 1,
                        fluxes[j][ei_idx], 1);
            // Calc diffusion of n with D tensor and n field
            for (unsigned int k = 0; k < nDim; ++k)
            {
                Vmath::Vvtvp(nPts, m_D[vc[j][k]].GetValue(), 1,
                             qfield[k][ei_idx], 1, fluxes[j][ei_idx], 1,
                             fluxes[j][ei_idx], 1);
            }
        }
    }

    if (nDim == 3)
    {
        Vmath::Vvtvvtm(nPts, b_unit[1], 1, qfield[2][pe_idx], 1, b_unit[2], 1,
                       qfield[1][pe_idx], 1, fluxes[0][pe_idx], 1);
        Vmath::Vvtvvtm(nPts, b_unit[2], 1, qfield[0][pe_idx], 1, b_unit[0], 1,
                       qfield[2][pe_idx], 1, fluxes[1][pe_idx], 1);
        Vmath::Vvtvvtm(nPts, b_unit[0], 1, qfield[1][pe_idx], 1, b_unit[1], 1,
                       qfield[0][pe_idx], 1, fluxes[2][pe_idx], 1);
    }
    else
    {
        Vmath::Vmul(nPts, b_unit[2], 1, qfield[1][pe_idx], 1, fluxes[0][pe_idx],
                    1);
        Vmath::Neg(nPts, fluxes[0][pe_idx], 1);
        Vmath::Vmul(nPts, b_unit[2], 1, qfield[0][pe_idx], 1, fluxes[1][pe_idx],
                    1);
    }

    CalcKappa(in_arr);
    CalcDiffTensor();
    for (unsigned int j = 0; j < nDim; ++j)
    {
        Vmath::Vmul(nPts, m_kcross, 1, fluxes[j][pe_idx], 1, fluxes[j][pe_idx],
                    1);
        for (unsigned int k = 0; k < nDim; ++k)
        {
            Vmath::Vvtvp(nPts, m_D[vc[j][k]].GetValue(), 1, qfield[k][pe_idx],
                         1, fluxes[j][pe_idx], 1, fluxes[j][pe_idx], 1);
        }
    }
}

void DoubleDiffusiveField::load_params()
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

bool DoubleDiffusiveField::v_PostIntegrate(int step)
{
    m_fields[0]->FwdTrans(m_fields[0]->GetPhys(), m_fields[0]->UpdateCoeffs());
    m_fields[1]->FwdTrans(m_fields[1]->GetPhys(), m_fields[1]->UpdateCoeffs());

    // Writes a step of the particle trajectory.

    return PlasmaSystem::v_PostIntegrate(step);
}

void DoubleDiffusiveField::v_ExtraFldOutput(
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
            i += 2;
        }
    }
}
} // namespace PENKNIFE