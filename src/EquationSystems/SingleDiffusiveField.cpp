#include "SingleDiffusiveField.hpp"

namespace PENKNIFE
{

/// Name of class
static std::string class_name;
std::string SingleDiffusiveField::class_name =
    SU::GetEquationSystemFactory().RegisterCreatorFunction(
        "SingleDiffusiveField", SingleDiffusiveField::create,
        "Solves for a single diffusive field (n) with anisotropy");
/**
 * @brief Creates an instance of this class.
 */
static SU::EquationSystemSharedPtr create(
    const LU::SessionReaderSharedPtr &session,
    const SD::MeshGraphSharedPtr &graph)
{
    SU::EquationSystemSharedPtr p =
        MemoryManager<SingleDiffusiveField>::AllocateSharedPtr(session, graph);
    p->InitObject();
    return p;
}

SingleDiffusiveField::SingleDiffusiveField(
    const LU::SessionReaderSharedPtr &session,
    const SD::MeshGraphSharedPtr &graph)
    : PlasmaSystem(session, graph)
{
    this->n_indep_fields = 0;
}

void SingleDiffusiveField::v_InitObject(bool DeclareFields)
{
    PlasmaSystem::v_InitObject(DeclareFields);
    this->ne = std::dynamic_pointer_cast<MR::DisContField>(m_fields[0]);

    int npoints = m_indfields[0]->GetNpoints();
    m_kpar      = Array<OneD, NekDouble>(npoints);
    m_kperp     = Array<OneD, NekDouble>(npoints);
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            m_D[i][j] = Array<OneD, NekDouble>(npoints);
        }
    }

    m_session->MatchSolverInfo("SpectralVanishingViscosity", "True",
                               m_useSpecVanVisc, false);
    m_session->LoadParameter("epsilon", m_epsilon, 1.0);
    if (m_useSpecVanVisc)
    {
        m_session->LoadParameter("SVVCutoffRatio", m_sVVCutoffRatio, 0.75);
        m_session->LoadParameter("SVVDiffCoeff", m_sVVDiffCoeff, 0.1);
    }

    // Setup diffusion object
    m_ode.DefineOdeRhs(&SingleDiffusiveField::DoOdeRhs, this);

    switch (m_projectionType)
    {
        case MultiRegions::eDiscontinuous:
        {
            std::string diffName;
            m_session->LoadSolverInfo("DiffusionType", diffName, "LDG");
            m_diffusion = SolverUtils::GetDiffusionFactory().CreateInstance(
                diffName, diffName);
            m_diffusion->SetFluxVector(&SingleDiffusiveField::GetFluxVectorDiff,
                                       this);
            m_diffusion->InitObject(m_session, m_indfields);
            break;
        }
        case MultiRegions::eGalerkin:
        {
            // Enable implicit solver for CG.
            m_ode.DefineImplicitSolve(&SingleDiffusiveField::ImplicitTimeIntCG,
                                      this);

            break;
        }
        default:
        {
            ASSERTL0(false, "Unknown projection scheme");
            break;
        }
    }
    if (this->particles_enabled)
    {
        std::vector<Sym<REAL>> src_syms;
        std::vector<int> src_components;
        for (auto &[k, v] : this->particle_sys->get_species())
        {
            this->src_fields.emplace_back(
                MemoryManager<MR::DisContField>::AllocateSharedPtr(
                    *std::dynamic_pointer_cast<MR::DisContField>(m_fields[0])));
            src_syms.push_back(Sym<REAL>(k + "_SOURCE_DENSITY"));
            src_components.push_back(0);
        }

        this->particle_sys->setup_evaluate_fields(this->E, this->B, this->ne,
                                                  this->Te, this->ve);

        this->particle_sys->finish_setup(this->src_fields, src_syms,
                                         src_components);
        this->diag_field = MemoryManager<MR::DisContField>::AllocateSharedPtr(
            *std::dynamic_pointer_cast<MR::DisContField>(m_fields[0]));
        this->particle_sys->diag_setup(this->diag_field);
    }
}

void SingleDiffusiveField::ImplicitTimeIntCG(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray,
    [[maybe_unused]] const NekDouble time, const NekDouble lambda)
{

    int npoints = m_indfields[0]->GetNpoints();
    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorLambda] = 1.0 / lambda / m_epsilon;
    if (m_useSpecVanVisc)
    {
        factors[StdRegions::eFactorSVVCutoffRatio] = m_sVVCutoffRatio;
        factors[StdRegions::eFactorSVVDiffCoeff]   = m_sVVDiffCoeff / m_epsilon;
    }

    // We solve ( \nabla^2 - HHlambda ) Y[i] = rhs [i]
    // inarray = input: \hat{rhs} -> output: \hat{Y}
    // outarray = output: nabla^2 \hat{Y}
    // where \hat = modal coeffs
    if (m_intScheme->GetIntegrationSchemeType() == LibUtilities::eImplicit)
    {
        for (const auto &[s, v] : GetIons())
        {
            int ni_idx = v.fields.at(field_to_index.at("n"));
            Vmath::Zero(npoints, outarray[ni_idx], 1);
        }

        for (auto &x : m_forcing)
        {
            x->Apply(m_indfields, inarray, outarray, time);
        }
    }

    for (const auto &[s, v] : GetIons())
    {
        int ni_idx = v.fields.at(field_to_index.at("n"));
        if (m_intScheme->GetIntegrationSchemeType() == LibUtilities::eImplicit)
        {
            // Multiply forcing term by -1 for definition of HelmSolve function
            Vmath::Smul(npoints, -1.0, outarray[ni_idx], 1, outarray[ni_idx],
                        1);

            // Multiply 1.0/timestep/lambda
            Vmath::Svtvp(npoints, -factors[StdRegions::eFactorLambda],
                         inarray[ni_idx], 1, outarray[ni_idx], 1,
                         outarray[ni_idx], 1);
        }
        else
        {
            // Multiply 1.0/timestep/lambda
            Vmath::Smul(npoints, -factors[StdRegions::eFactorLambda],
                        inarray[ni_idx], 1, outarray[ni_idx], 1);
        }
        CalcDiffTensor(s);
        StdRegions::VarCoeffMap varcoeffs;
        MultiRegions::VarFactorsMap varfactors =
            MultiRegions::NullVarFactorsMap;

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                varcoeffs[vc[i][j]] = m_D[i][j];
            }
        }

        // Solve a system of equations with Helmholtz solver
        auto key = m_indfields[ni_idx]->HelmSolve(
            outarray[ni_idx], m_indfields[ni_idx]->UpdateCoeffs(), factors,
            varcoeffs, varfactors);

        if (key.GetMatrixType() == StdRegions::eHelmholtz ||
            key.GetMatrixType() == StdRegions::eHelmholtzGJP)
        {
            m_indfields[ni_idx]->UnsetGlobalLinSys(key, true);
        }

        m_indfields[ni_idx]->BwdTrans(m_indfields[ni_idx]->GetCoeffs(),
                                      outarray[ni_idx]);

        m_indfields[ni_idx]->SetPhysState(false);
    }
}

void SingleDiffusiveField::CalcKPar(int f)
{
    int npoints = m_fields[0]->GetNpoints();
    if (m_session->DefinesParameter("k_par"))
    {
        double k = m_session->GetParameter("k_par");
        Vmath::Fill(npoints, k, m_kpar, 1);
    }
    else
    {
        double Z   = m_ions[f].charge;
        int ni_idx = m_ions[f].fields.at(field_to_index.at("n"));

        Vmath::Fill(npoints, this->k_par / (Z * Z), m_kpar, 1);
        Vmath::Vdiv(npoints, m_kpar, 1, m_indfields[ni_idx]->GetPhys(), 1,
                    m_kpar, 1);
    }
}

void SingleDiffusiveField::CalcKPerp(int f)
{
    int npoints = m_fields[0]->GetNpoints();
    if (m_session->DefinesParameter("k_perp"))
    {
        double k = m_session->GetParameter("k_perp");
        Vmath::Fill(npoints, k, m_kperp, 1);
    }
    else
    {
        double Z   = m_ions[f].charge;
        double A   = m_ions[f].mass;
        int ni_idx = m_ions[f].fields.at(field_to_index.at("n"));

        Vmath::Fill(npoints, this->k_perp * Z * Z * std::sqrt(A), m_kperp, 1);
        Vmath::Vmul(npoints, m_kperp, 1, m_indfields[ni_idx]->GetPhys(), 1,
                    m_kperp, 1);
        Vmath::Vdiv(npoints, m_kperp, 1, this->mag_B, 1, m_kperp, 1);
    }
}

void SingleDiffusiveField::CalcKPerpAnomalous(int f)
{
    int npoints = m_fields[0]->GetNpoints();

    for (int p = 0; p < npoints; ++p)
    {
        m_kperp[p] = this->k_perp / std::sqrt(this->mag_B[p]);
    }
}

void SingleDiffusiveField::CalcDiffTensor(int f)
{
    int npoints = m_fields[0]->GetNpoints();

    CalcKPar(f);
    CalcKPerp(f);

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < npoints; k++)
            {
                m_D[i][j][k] =
                    (m_kpar[k] - m_kperp[k]) * b_unit[i][k] * b_unit[j][k];
                if (i == j)
                {
                    m_D[i][j][k] += m_kperp[k];
                }
            }
        }
    }
}

/**
 * @brief Populate rhs array ( @p out_arr ) with diffused quantitites
 * @param in_arr physical values of all fields
 * @param[out] out_arr output array (RHSs of time integration equations)
 */
void SingleDiffusiveField::DoOdeRhs(
    const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
    Array<OneD, Array<OneD, NekDouble>> &out_arr, const NekDouble time)
{
    if (m_explicitDiffusion || m_projectionType == MR::eDiscontinuous)
    {
        size_t nvariables = in_arr.size();
        m_diffusion->Diffuse(nvariables, m_indfields, in_arr, out_arr);
    }
    else
    {
        // RHS should be set to zero.
        for (int i = 0; i < out_arr.size(); ++i)
        {
            Vmath::Zero(out_arr[i].size(), &out_arr[i][0], 1);
        }
    }

    if (this->particles_enabled)
    {
        for (int i = 0; i < this->particle_sys->get_species().size(); ++i)
        {
            Vmath::Vadd(out_arr[i].size(), out_arr[i], 1,
                        this->src_fields[i]->GetPhys(), 1, out_arr[i], 1);
        }
    }

    // Add forcing terms
    for (auto &x : m_forcing)
    {
        x->Apply(m_fields, in_arr, out_arr, time);
    }
}

/**
 * @brief Construct the flux vector for the anisotropic diffusion problem.
 */
void SingleDiffusiveField::GetFluxVectorDiff(
    const Array<OneD, Array<OneD, NekDouble>> &in_arr,
    const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &qfield,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fluxes)
{
    unsigned int nDim = qfield.size();
    unsigned int nFld = qfield[0].size();
    unsigned int nPts = qfield[0][0].size();

    for (auto &[k, v] : this->GetIons())
    {
        int f = v.fields[field_to_index["n"]];
        CalcDiffTensor(k);
        for (unsigned int j = 0; j < nDim; ++j)
        {
            // Calc diffusion of n with D tensor
            Vmath::Vmul(nPts, m_D[j][0], 1, qfield[0][f], 1, fluxes[j][f], 1);
            for (unsigned int k = 1; k < nDim; ++k)
            {
                Vmath::Vvtvp(nPts, m_D[j][k], 1, qfield[k][f], 1, fluxes[j][f],
                             1, fluxes[j][f], 1);
            }
        }
    }
}

void SingleDiffusiveField::load_params()
{
    PlasmaSystem::load_params();
    NekDouble k_c, lambda, T_bg;

    m_session->LoadParameter("k_c", k_c);
    m_session->LoadParameter("lambda", lambda);
    m_session->LoadParameter("T_bg", T_bg);

    k_par = 6.0 * k_c * (sqrt(2.0 * pow(M_PI, 3)) / lambda) *
            constants::epsilon_0 * constants::epsilon_0 * constants::c *
            pow(this->Tnorm * T_bg, 2.5) / sqrt(constants::m_e);
    // Correct for microns in epsilon_0 and density scale
    k_par *= 1e12 / this->Nnorm;
    // Convert to solver length and time scale
    k_par /= (this->omega_c * this->mesh_length * this->mesh_length);
    // multiply k_par by Z^-2 n^-1 in solver

    k_perp = lambda / (6.0 * sqrt(this->Tnorm * T_bg * pow(M_PI, 3.0))) *
             (sqrt(constants::m_p) / constants::c) *
             pow((constants::e / constants::epsilon_0_si), 2);
    k_perp *= this->Nnorm;
    k_perp /= (this->omega_c * this->mesh_length * this->mesh_length);
    // multiply k_perp by A^0.5 Z^2 n B^-1 in solver
}

bool SingleDiffusiveField::v_PostIntegrate(int step)
{
    m_fields[0]->FwdTrans(m_fields[0]->GetPhys(), m_fields[0]->UpdateCoeffs());

    if (this->particles_enabled)
        this->particle_sys->diag_project();
    // Writes a step of the particle trajectory.

    return PlasmaSystem::v_PostIntegrate(step);
}

void SingleDiffusiveField::v_ExtraFldOutput(
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
            Array<OneD, NekDouble> SrcFwd(nCoeffs);
            m_fields[0]->FwdTransLocalElmt(this->src_fields[i++]->GetPhys(),
                                           SrcFwd);
            fieldcoeffs.push_back(SrcFwd);
        }
        variables.push_back("DIAGNOSTIC");
        Array<OneD, NekDouble> SrcFwd(nCoeffs);
        m_fields[0]->FwdTransLocalElmt(this->diag_field->GetPhys(), SrcFwd);
        fieldcoeffs.push_back(SrcFwd);
    }
}
} // namespace PENKNIFE