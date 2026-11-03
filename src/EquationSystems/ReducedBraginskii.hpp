#ifndef REDUCEDBRAGINSKII_HPP
#define REDUCEDBRAGINSKII_HPP

#include "../Misc/VariableConverter.hpp"
#include "PlasmaSystem.hpp"
#include <SolverCore/AdvDiffusion/AdvWeakDGDiffusionIP/AdvWeakDGDiffusionIPOp.hpp>
#include <SolverCore/TraceFlux/TraceFluxOp.hpp>

namespace PENKNIFE
{
// [{n, mnv, 1.5p}_i,..., 1.5p_e]
class ReducedBraginskii : public PlasmaSystem
{
public:
    friend class MemoryManager<ReducedBraginskii>;

    /// Name of class
    static std::string class_name;

    /**
     * @brief Create an instance of this class and initialise it.
     */
    static SC::EquationSystemSharedPtr create(
        const LU::SessionReaderSharedPtr &session,
        const SD::MeshGraphSharedPtr &graph)
    {
        SC::EquationSystemSharedPtr p =
            MemoryManager<ReducedBraginskii>::AllocateSharedPtr(session, graph);
        p->InitObject();
        return p;
    }
    ~ReducedBraginskii() override = default;

protected:
    ReducedBraginskii(const LU::SessionReaderSharedPtr &session,
                      const SD::MeshGraphSharedPtr &graph);

    void SetUpBoundaryConditions();
    virtual void v_InitObject(bool DeclareFields = true) override;
    virtual void v_InitialiseFields() override;
    virtual void v_InitialiseOperators() override;

    SC::TraceFluxOp<double>::ViscousTraceForm ReadViscousTraceForm();

    void v_SetInitialConditions(NekDouble init_time, bool dump_ICs,
                                const int domain) override;
    bool v_PreIntegrate(int step) override;
    bool v_PostIntegrate(int step) override;
    void DoOdeRhs(LU::Field<double, FieldState::Phys> &in,
                  LU::Field<double, FieldState::Phys> &out,
                  [[maybe_unused]] const double &time, const double &dt);

    void DoParticles(LU::Field<double, FieldState::Phys> &in,
                     LU::Field<double, FieldState::Phys> &out);

    void ComputeE();

    void CalcVelocities(LU::Field<double, FieldState::Phys> &in);

    void InitBraginskii();

    void CalcCollisionFrequencies(
        const Array<OneD, Array<OneD, NekDouble>> &in_arr);

    void CalcNeutralSources_nvp(
        const double m, int ee_idx, int ni_idx, int vi_idx, int pi_idx,
        int nn_idx, int vn_idx, int pn_idx,
        const Array<OneD, Array<OneD, NekDouble>> &inarray,
        const Array<OneD, NekDouble> &ne,
        Array<OneD, Array<OneD, NekDouble>> &outarray,
        Array<OneD, NekDouble> &Spe);
    void CalcNeutralSources_nv(
        const double m, int ee_idx, int ni_idx, int vi_idx, int nn_idx,
        int vn_idx, const Array<OneD, Array<OneD, NekDouble>> &inarray,
        const Array<OneD, NekDouble> &ne,
        Array<OneD, Array<OneD, NekDouble>> &outarray,
        Array<OneD, NekDouble> &Spe);
    void CalcNeutralSources_nv(
        const double m, int ee_idx, int ni_idx, int vi_idx, int pi_idx,
        int nn_idx, int vn_idx,
        const Array<OneD, Array<OneD, NekDouble>> &inarray,
        const Array<OneD, NekDouble> &ne,
        Array<OneD, Array<OneD, NekDouble>> &outarray,
        Array<OneD, NekDouble> &Spe);
    void CalcNeutralSources_np(
        const double m, int ee_idx, int ni_idx, int pi_idx, int nn_idx,
        int pn_idx, const Array<OneD, Array<OneD, NekDouble>> &inarray,
        const Array<OneD, NekDouble> &ne,
        Array<OneD, Array<OneD, NekDouble>> &outarray,
        Array<OneD, NekDouble> &Spe);
    void CalcNeutralSources_n(
        const double m, int ee_idx, int ni_idx, int nn_idx,
        const Array<OneD, Array<OneD, NekDouble>> &inarray,
        const Array<OneD, NekDouble> &ne,
        Array<OneD, Array<OneD, NekDouble>> &outarray,
        Array<OneD, NekDouble> &Spe);

    void AddNeutralSources(const Array<OneD, Array<OneD, NekDouble>> &inarray,
                           Array<OneD, Array<OneD, NekDouble>> &outarray);

    void v_ExtraFldOutput(std::vector<Array<OneD, NekDouble>> &fieldcoeffs,
                          std::vector<std::string> &variables) override;

private:
    /// Velocities

    // Setup workspaces
    LU::Field<double, FieldState::Phys> m_advectionVel;
    /// The advection velocity on the local trace, which the trace flux
    /// operator reads through a pointer; kept here so that it outlives it.
    LU::Field<double, FieldState::Phys> m_traceAdvVel;

    std::map<double, std::shared_ptr<PreconOp<double>>> m_preconOp;
    std::shared_ptr<AdvDiffusionOp<double>> m_advDiffusionOp;
    std::shared_ptr<PhysTraceExtractOp<double>> m_physTraceExtractOp;

    Array<OneD, NekDouble> trace_b_norm;
    /// Riemann solver type (used for all advection terms)
    std::string riemann_solver_type;
    /// Riemann solver object used in electron advection
    SU::RiemannSolverSharedPtr riemann_solver;
    /// Advection object used in the electron density equation
    SU::AdvectionSharedPtr m_advection;
    /// Advection type
    std::string adv_type;

    // For Advection
    std::vector<int> advected_fields;
    Array<OneD, MR::ExpListSharedPtr> m_advfields;

    // For Diffusion
    Array<OneD, MR::ExpListSharedPtr> m_difffields;
    //
    StdRegions::ConstFactorMap m_factors;

    ClosureSharedPtr m_closure;

    Array<OneD, Array<OneD, NekDouble>> friction;

    Array<OneD, NekDouble> kIZ;
    Array<OneD, NekDouble> kCX;
    Array<OneD, NekDouble> krec;
};

} // namespace PENKNIFE
#endif