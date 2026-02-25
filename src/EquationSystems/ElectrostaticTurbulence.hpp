#ifndef ELECTROSTATICTURBULENCE_HPP
#define ELECTROSTATICTURBULENCE_HPP
#include "../Misc/VariableConverter.hpp"
#include "PlasmaSystem.hpp"

namespace PENKNIFE
{
// [omega, {n, mnv, 1.5p}_e, {n, mnv, 1.5p}_i, ...]
class ElectrostaticTurbulence : public PlasmaSystem
{
public:
    friend class MemoryManager<ElectrostaticTurbulence>;

    /// Name of class
    static std::string class_name;

    /**
     * @brief Create an instance of this class and initialise it.
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
    ~ElectrostaticTurbulence() override = default;

protected:
    ElectrostaticTurbulence(const LU::SessionReaderSharedPtr &session,
                            const SD::MeshGraphSharedPtr &graph);
    void v_InitObject(bool DeclareFields = true) override;
    void v_SetInitialConditions(NekDouble init_time, bool dump_ICs,
                                const int domain) override;

    bool v_PostIntegrate(int step) override;
    void DoOdeRhs(const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                  Array<OneD, Array<OneD, NekDouble>> &outarray,
                  const NekDouble time);

    /// Advection functions
    void DoAdvection(const Array<OneD, Array<OneD, NekDouble>> &inarray,
                     Array<OneD, Array<OneD, NekDouble>> &outarray,
                     const NekDouble time,
                     const Array<OneD, Array<OneD, NekDouble>> &pFwd,
                     const Array<OneD, Array<OneD, NekDouble>> &pBwd);

    void DoExtra(const Array<OneD, Array<OneD, NekDouble>> &inarray,
                 Array<OneD, Array<OneD, NekDouble>> &outarray);

    void DoParticles(const Array<OneD, Array<OneD, NekDouble>> &inarray,
                     Array<OneD, Array<OneD, NekDouble>> &outarray);

    void CalcInitPhi();
    void CalcInitOmega();
    void SolvePhi(const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                  [[maybe_unused]] const Array<OneD, NekDouble> &ne);
    void ComputeE();
    void ComputevExB();
    void AddDriftVelocities(
        const Array<OneD, Array<OneD, NekDouble>> &inarray,
        const Array<OneD, NekDouble> &ne,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &adv_vel);
    void CalcOmegaFlux(const Array<OneD, Array<OneD, NekDouble>> &inarray,
                       Array<OneD, Array<OneD, NekDouble>> &omega_flux);

    void CalcVelocities(const Array<OneD, Array<OneD, NekDouble>> &inarray);

    Array<OneD, Array<OneD, NekDouble>> &GetAdvVelNorm();
    Array<OneD, NekDouble> &GetOmegaFlux();

    // Advective Flux vector
    void GetFluxVector(
        const Array<OneD, Array<OneD, NekDouble>> &field_vals,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fluxes);

    /// Diffusion functions
    void DoDiffusion(const Array<OneD, Array<OneD, NekDouble>> &inarray,
                     Array<OneD, Array<OneD, NekDouble>> &outarray,
                     const Array<OneD, Array<OneD, NekDouble>> &pFwd,
                     const Array<OneD, Array<OneD, NekDouble>> &pBwd);

    void CalcKPar();
    void CalcKPerp();
    void CalcDiffTensor();
    void CalcKappaPar();
    void CalcKappaPerp();
    void CalcKappaTensor();
    // Diffusive Flux vector
    void GetFluxVectorDiff(
        const Array<OneD, Array<OneD, NekDouble>> &in_arr,
        const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &q_field,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fluxes);

    void AddNeutralSources(const Array<OneD, Array<OneD, NekDouble>> &in_arr,
                           Array<OneD, Array<OneD, NekDouble>> &outarray);

    // Functions for the Implicit Solve
    void DoOdeImplicitRhs(
        const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
        Array<OneD, Array<OneD, NekDouble>> &out_arr, const NekDouble time);

    void DoOdeRhsCoeff(const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                       Array<OneD, Array<OneD, NekDouble>> &outarray,
                       const NekDouble time);
    void DoAdvectionCoeff(const Array<OneD, Array<OneD, NekDouble>> &inarray,
                          Array<OneD, Array<OneD, NekDouble>> &outarray,
                          const NekDouble time,
                          const Array<OneD, Array<OneD, NekDouble>> &pFwd,
                          const Array<OneD, Array<OneD, NekDouble>> &pBwd);
    void DoParticlesCoeff(const Array<OneD, Array<OneD, NekDouble>> &inarray,
                          Array<OneD, Array<OneD, NekDouble>> &out_arr);
    void DoDiffusionCoeff(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray,
        const Array<OneD, const Array<OneD, NekDouble>> &pFwd,
        const Array<OneD, const Array<OneD, NekDouble>> &pBwd);

    void load_params() override;

    void v_ExtraFldOutput(std::vector<Array<OneD, NekDouble>> &fieldcoeffs,
                          std::vector<std::string> &variables) override;

private:
    int ee_idx;
    int omega_idx;
    int phi_idx;

    std::vector<int> ni_src_idx;
    std::vector<int> vi_src_idx;
    std::vector<int> ei_src_idx;

    /// Hasegawa-Wakatani α
    NekDouble alpha;
    /// Hasegawa-Wakatani κ
    NekDouble kappa;

    /// Velocities
    /// Storage for ExB drift velocity
    Array<OneD, Array<OneD, NekDouble>> v_ExB;
    // Electron parallel velocity
    Array<OneD, NekDouble> v_e_par;
    // Ion parallel velocities
    std::vector<Array<OneD, NekDouble>> v_i_par;
    // Electron diamagnetic drift velocity
    Array<OneD, Array<OneD, NekDouble>> v_de;
    // Ion diamagnetic drift velocities
    std::vector<Array<OneD, Array<OneD, NekDouble>>> v_di;
    // Per field advection velocities
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> adv_vel;
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> adv_vel_trace;
    // Per field advection velocities normal to trace elements
    Array<OneD, Array<OneD, NekDouble>> trace_vel_norm;

    // For Advection
    std::vector<int> advected_fields;
    Array<OneD, MR::ExpListSharedPtr> m_advfields;

    // Vorticity Flux
    Array<OneD, Array<OneD, NekDouble>> omega_flux;
    Array<OneD, Array<OneD, NekDouble>> omega_flux_trace;
    // Vorticity flux normal to trace elements
    Array<OneD, NekDouble> omega_flux_norm;
    Array<OneD, NekDouble> trace_b_norm;

    MR::ExpListSharedPtr phi;

    StdRegions::VarCoeffMap m_phi_varcoeff;

    /// Whether the Boussinesq approximation is used for the vorticity
    bool m_boussinesq;
    /// Riemann solver type (used for all advection terms)
    std::string riemann_solver_type;
    /// Riemann solver object used in electron advection
    SU::RiemannSolverSharedPtr riemann_solver;
    /// Advection object used in the electron density equation
    SU::AdvectionSharedPtr m_advection;
    /// Advection type
    std::string adv_type;

    // For Diffusion
    // workaround for bug in DiffusionLDG
    Array<OneD, MR::ExpListSharedPtr> m_difffields;
    //
    StdRegions::ConstFactorMap m_factors;

    Array<OneD, NekDouble> m_kperp;
    Array<OneD, NekDouble> m_kpar;
    StdRegions::VarCoeffMap m_D;

    Array<OneD, NekDouble> m_kappaperp;
    Array<OneD, NekDouble> m_kappapar;
    StdRegions::VarCoeffMap m_kappa;
};

} // namespace PENKNIFE
#endif