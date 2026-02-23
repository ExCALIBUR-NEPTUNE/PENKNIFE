#ifndef REDUCEDBRAGINSKII_HPP
#define REDUCEDBRAGINSKII_HPP
#include "../Misc/VariableConverter.hpp"
#include "PlasmaSystem.hpp"

namespace PENKNIFE
{
// [omega, {n, mnv, 1.5p}_e, {n, mnv, 1.5p}_i, ...]
class ReducedBraginskii : public PlasmaSystem
{
public:
    friend class MemoryManager<ReducedBraginskii>;

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
            MemoryManager<ReducedBraginskii>::AllocateSharedPtr(session, graph);
        p->InitObject();
        return p;
    }
    ~ReducedBraginskii() override = default;

protected:
    ReducedBraginskii(const LU::SessionReaderSharedPtr &session,
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

    void DoParticles(const Array<OneD, Array<OneD, NekDouble>> &inarray,
                     Array<OneD, Array<OneD, NekDouble>> &outarray);

    void ComputeE();

    void CalcVelocities(const Array<OneD, Array<OneD, NekDouble>> &inarray,
                        const Array<OneD, NekDouble> &ne,
                        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &);

    Array<OneD, Array<OneD, NekDouble>> &GetAdvVelNorm();

    // Advective Flux vector
    void GetFluxVector(
        const Array<OneD, Array<OneD, NekDouble>> &field_vals,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fluxes);

    /// Diffusion functions
    void DoDiffusion(const Array<OneD, Array<OneD, NekDouble>> &inarray,
                     Array<OneD, Array<OneD, NekDouble>> &outarray,
                     const Array<OneD, Array<OneD, NekDouble>> &pFwd,
                     const Array<OneD, Array<OneD, NekDouble>> &pBwd);

    void CalcK(const Array<OneD, Array<OneD, NekDouble>> &in_arr, int f);
    void CalcKappa(const Array<OneD, Array<OneD, NekDouble>> &in_arr, int f);
    void CalcKappa(const Array<OneD, Array<OneD, NekDouble>> &in_arr);
    void CalcDiffTensor();
    // Diffusive Flux vector
    void GetFluxVectorDiff(
        const Array<OneD, Array<OneD, NekDouble>> &in_arr,
        const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &q_field,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fluxes);

    void CalcNeutralSources(const double m, int pe_idx, int ni_idx, int vi_idx,
                            int pi_idx, int nn_idx, int vn_idx, int pn_idx,
                            const Array<OneD, Array<OneD, NekDouble>> &inarray,
                            const Array<OneD, NekDouble> &ne,
                            Array<OneD, Array<OneD, NekDouble>> &outarray,
                            Array<OneD, NekDouble> &Spe);
    void CalcNeutralSources(const double m, int pe_idx, int ni_idx, int vi_idx,
                            int nn_idx, int vn_idx,
                            const Array<OneD, Array<OneD, NekDouble>> &inarray,
                            const Array<OneD, NekDouble> &ne,
                            Array<OneD, Array<OneD, NekDouble>> &outarray,
                            Array<OneD, NekDouble> &Spe);
    void CalcNeutralSources(const double m, int pe_idx, int ni_idx, int vi_idx,
                            int pi_idx, int nn_idx, int vn_idx,
                            const Array<OneD, Array<OneD, NekDouble>> &inarray,
                            const Array<OneD, NekDouble> &ne,
                            Array<OneD, Array<OneD, NekDouble>> &outarray,
                            Array<OneD, NekDouble> &Spe);
    void CalcNeutralSources(const double m, int pe_idx, int ni_idx, int pi_idx,
                            int nn_idx, int pn_idx,
                            const Array<OneD, Array<OneD, NekDouble>> &inarray,
                            const Array<OneD, NekDouble> &ne,
                            Array<OneD, Array<OneD, NekDouble>> &outarray,
                            Array<OneD, NekDouble> &Spe);
    void CalcNeutralSources(const double m, int pe_idx, int ni_idx, int nn_idx,
                            const Array<OneD, Array<OneD, NekDouble>> &inarray,
                            const Array<OneD, NekDouble> &ne,
                            Array<OneD, Array<OneD, NekDouble>> &outarray,
                            Array<OneD, NekDouble> &Spe);

    void AddNeutralSources(const Array<OneD, Array<OneD, NekDouble>> &inarray,
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
    int pe_idx;

    std::vector<int> ni_src_idx;
    std::vector<int> vi_src_idx;
    std::vector<int> pi_src_idx;

    /// Velocities

    // Electron parallel velocity
    Array<OneD, NekDouble> v_e_par;
    // Ion parallel velocities
    std::vector<Array<OneD, NekDouble>> v_i_par;

    // Per field advection velocities
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> adv_vel;
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> adv_vel_trace;
    // Per field advection velocities normal to trace elements
    Array<OneD, Array<OneD, NekDouble>> trace_vel_norm;

    Array<OneD, NekDouble> trace_b_norm;
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

    NekDouble k_par;
    NekDouble k_perp;
    NekDouble k_cross;
    NekDouble kappa_i_par;
    NekDouble kappa_i_perp;
    NekDouble kappa_i_cross;
    NekDouble kappa_e_par;
    NekDouble kappa_e_perp;
    NekDouble kappa_e_cross;
    NekDouble k_ci;
    NekDouble k_ce;

    Array<OneD, NekDouble> m_kpar;
    Array<OneD, NekDouble> m_kperp;
    Array<OneD, NekDouble> m_kcross;

    StdRegions::VarCoeffMap m_D;

    Array<OneD, NekDouble> kIZ;
    Array<OneD, NekDouble> kCX;
    Array<OneD, NekDouble> krec;
};

} // namespace PENKNIFE
#endif